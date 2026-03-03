use anyhow::Result;
use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt};
use tracing::info;

use crate::metadata::{
    decode_completion_msg, decode_handshake, LambdaMetadata, COMPLETION_MSG,
    COMPLETION_MSG_SIZE, HANDSHAKE_SIZE,
};

const BLOCK_HEADER_SIZE: usize = 24;

struct BlockHeader {
    #[allow(dead_code)]
    block_id: i64,
    block_size: u64,
    #[allow(dead_code)]
    is_last: i8,
    raw: [u8; BLOCK_HEADER_SIZE],
}

#[derive(Debug, Clone)]
pub enum ReceiverProgress {
    Chunk {
        success: bool,
        num_of_completed_chunks: u64,
        num_of_partial_bytes: u64,
    },
    Raw {
        success: bool,
        num_of_recv_bytes: u64,
    },
}

impl ReceiverProgress {
    pub fn success(&self) -> bool {
        match self {
            Self::Chunk { success, .. } | Self::Raw { success, .. } => *success,
        }
    }
}

pub async fn read_metadata<R>(reader: &mut R) -> Result<LambdaMetadata>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0u8; HANDSHAKE_SIZE];
    reader.read_exact(&mut buf).await?;
    let metadata = decode_handshake(&buf)?;
    info!(
        is_stateless = metadata.is_stateless,
        chunk_start_id = metadata.chunk_start_id,
        "[receiver.rs] Decoded metadata"
    );
    Ok(metadata)
}

pub struct RawReader {
    to_stdout: bool,
    fifo_name: String,
    file: Option<File>,
    stdout: tokio::io::Stdout,
    total_received: u64,
    carry: Vec<u8>,
}

impl RawReader {
    pub async fn new(fifo_name: &str) -> Result<Self> {
        let to_stdout = fifo_name == "-";
        let file = if to_stdout {
            None
        } else {
            Some(File::create(fifo_name).await?)
        };
        Ok(Self {
            to_stdout,
            fifo_name: fifo_name.to_string(),
            file,
            stdout: tokio::io::stdout(),
            total_received: 0,
            carry: Vec::new(),
        })
    }

    pub async fn read_from<R>(&mut self, reader: &mut R) -> Result<ReceiverProgress>
    where
        R: AsyncRead + Unpin,
    {
        info!("[receiver.rs] Starting to read raw bytes to {}", if self.to_stdout { "stdout" } else { self.fifo_name.as_str() });
        let mut buf = [0u8; 8192];
        let mut recv_bytes = 0u64;
        let mut skip_remaining = self.total_received;
        self.carry.clear();
        if skip_remaining > 0 {
            info!(
                skip_remaining,
                "[receiver.rs] RawReader skipping previously received bytes"
            );
        }

        loop {
            match reader.read(&mut buf).await {
                Ok(0) => {
                    self.total_received = self.total_received.saturating_add(recv_bytes);
                    return Ok(ReceiverProgress::Raw {
                        success: false,
                        num_of_recv_bytes: recv_bytes,
                    });
                }
                Ok(n) => {
                    if !self.carry.is_empty() {
                        let mut combined = Vec::with_capacity(self.carry.len() + n);
                        combined.extend_from_slice(&self.carry);
                        combined.extend_from_slice(&buf[..n]);

                        if let Some(pos) = find_completion_msg(&combined) {
                            info!(
                                pos,
                                "[receiver.rs] RawReader detected completion message (combined)"
                            );
                            let payload = &combined[..pos];
                            let (_discarded, written) = discard_and_write(
                                payload,
                                &mut skip_remaining,
                                self.to_stdout,
                                &mut self.file,
                                &mut self.stdout,
                            )
                            .await?;
                            recv_bytes += written as u64;
                            self.total_received = self.total_received.saturating_add(recv_bytes);
                            return Ok(ReceiverProgress::Raw {
                                success: true,
                                num_of_recv_bytes: recv_bytes,
                            });
                        }

                        if combined.len() >= COMPLETION_MSG_SIZE {
                            let keep = COMPLETION_MSG_SIZE - 1;
                            let write_len = combined.len() - keep;
                            if write_len > 0 {
                                let payload = &combined[..write_len];
                                let (_discarded, written) = discard_and_write(
                                    payload,
                                    &mut skip_remaining,
                                    self.to_stdout,
                                    &mut self.file,
                                    &mut self.stdout,
                                )
                                .await?;
                                recv_bytes += written as u64;
                            }
                            self.carry.clear();
                            self.carry.extend_from_slice(&combined[write_len..]);
                        } else {
                            self.carry = combined;
                        }
                    } else {
                        if let Some(pos) = find_completion_msg(&buf[..n]) {
                            info!(
                                pos,
                                "[receiver.rs] RawReader detected completion message (buffer)"
                            );
                            let payload = &buf[..pos];
                            let (_discarded, written) = discard_and_write(
                                payload,
                                &mut skip_remaining,
                                self.to_stdout,
                                &mut self.file,
                                &mut self.stdout,
                            )
                            .await?;
                            recv_bytes += written as u64;
                            self.total_received = self.total_received.saturating_add(recv_bytes);
                            return Ok(ReceiverProgress::Raw {
                                success: true,
                                num_of_recv_bytes: recv_bytes,
                            });
                        }

                        if n >= COMPLETION_MSG_SIZE {
                            let keep = COMPLETION_MSG_SIZE - 1;
                            let write_len = n - keep;
                            if write_len > 0 {
                                let payload = &buf[..write_len];
                                let (_discarded, written) = discard_and_write(
                                    payload,
                                    &mut skip_remaining,
                                    self.to_stdout,
                                    &mut self.file,
                                    &mut self.stdout,
                                )
                                .await?;
                                recv_bytes += written as u64;
                            }
                            self.carry.clear();
                            self.carry.extend_from_slice(&buf[write_len..n]);
                        } else {
                            self.carry.clear();
                            self.carry.extend_from_slice(&buf[..n]);
                        }
                    }
                }
                Err(_) => {
                    self.total_received = self.total_received.saturating_add(recv_bytes);
                    return Ok(ReceiverProgress::Raw {
                        success: false,
                        num_of_recv_bytes: recv_bytes,
                    });
                }
            }
        }
    }
}

pub struct ChunkReader {
    to_stdout: bool,
    file: Option<File>,
    stdout: tokio::io::Stdout,
    completed_chunks: u64,
    partial_bytes: u64,
    partial_header_written: bool,
}

impl ChunkReader {
    pub async fn new(fifo_name: &str) -> Result<Self> {
        let to_stdout = fifo_name == "-";
        let file = if to_stdout {
            None
        } else {
            Some(File::create(fifo_name).await?)
        };
        Ok(Self {
            to_stdout,
            file,
            stdout: tokio::io::stdout(),
            completed_chunks: 0,
            partial_bytes: 0,
            partial_header_written: false,
        })
    }

    pub async fn read_from<R>(&mut self, reader: &mut R) -> Result<ReceiverProgress>
    where
        R: AsyncRead + Unpin,
    {
        let mut buf = [0u8; 8192];
        let mut attempt_completed = 0u64;
        let mut skip_chunks = self.completed_chunks;
        let mut skip_chunk_bytes = self.partial_bytes;
        let mut header_already_written = self.partial_header_written;
        if skip_chunks > 0 || skip_chunk_bytes > 0 {
            info!(
                skip_chunks,
                skip_chunk_bytes,
                "[receiver.rs] ChunkReader skipping previously received data"
            );
        }

        loop {
            let header = match read_block_header_or_completion(reader).await {
                Ok(Some(header)) => header,
                Ok(None) => {
                    return Ok(ReceiverProgress::Chunk {
                        success: true,
                        num_of_completed_chunks: attempt_completed,
                        num_of_partial_bytes: 0,
                    });
                }
                Err(_) => {
                    return Ok(ReceiverProgress::Chunk {
                        success: false,
                        num_of_completed_chunks: attempt_completed,
                        num_of_partial_bytes: 0,
                    });
                }
            };

            info!(
                block_id = header.block_id,
                block_size = header.block_size,
                "[receiver.rs] ChunkReader received block header"
            );

            if skip_chunks == 0 && skip_chunk_bytes == 0 && !header_already_written {
                if self.to_stdout {
                    self.stdout.write_all(&header.raw).await?;
                } else if let Some(f) = &mut self.file {
                    f.write_all(&header.raw).await?;
                }
                header_already_written = true;
            }

            let mut remaining = header.block_size;
            let mut current_written: u64 = 0;
            let initial_skip = if skip_chunks == 0 { skip_chunk_bytes } else { 0 };

            while remaining > 0 {
                let to_read = (buf.len() as u64).min(remaining) as usize;
                let n = match reader.read(&mut buf[..to_read]).await {
                    Ok(0) => {
                        let partial = if skip_chunks == 0 {
                            initial_skip + current_written
                        } else {
                            0
                        };
                        if skip_chunks == 0 {
                            self.partial_bytes = partial;
                            self.partial_header_written = header_already_written;
                        }
                        return Ok(ReceiverProgress::Chunk {
                            success: false,
                            num_of_completed_chunks: attempt_completed,
                            num_of_partial_bytes: partial,
                        });
                    }
                    Ok(n) => n,
                    Err(_) => {
                        let partial = if skip_chunks == 0 {
                            initial_skip + current_written
                        } else {
                            0
                        };
                        if skip_chunks == 0 {
                            self.partial_bytes = partial;
                            self.partial_header_written = header_already_written;
                        }
                        return Ok(ReceiverProgress::Chunk {
                            success: false,
                            num_of_completed_chunks: attempt_completed,
                            num_of_partial_bytes: partial,
                        });
                    }
                };

                if skip_chunks == 0 {
                    let mut start = 0usize;
                    if skip_chunk_bytes > 0 {
                        let discard = (skip_chunk_bytes as usize).min(n);
                        skip_chunk_bytes -= discard as u64;
                        start = discard;
                    }
                    if start < n {
                        let data = &buf[start..n];
                        if self.to_stdout {
                            self.stdout.write_all(data).await?;
                        } else if let Some(f) = &mut self.file {
                            f.write_all(data).await?;
                        }
                        current_written += data.len() as u64;
                    }
                }
                remaining = remaining.saturating_sub(n as u64);
            }

            if skip_chunks > 0 {
                skip_chunks -= 1;
                header_already_written = false;
                info!("[receiver.rs] ChunkReader skipped full chunk");
                continue;
            }

            attempt_completed += 1;
            self.completed_chunks = self.completed_chunks.saturating_add(1);
            self.partial_bytes = 0;
            self.partial_header_written = false;
            skip_chunk_bytes = 0;
            header_already_written = false;
        }
    }
}

async fn read_block_header_or_completion<R>(reader: &mut R) -> Result<Option<BlockHeader>>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0u8; BLOCK_HEADER_SIZE];
    let mut read = 0usize;

    while read < BLOCK_HEADER_SIZE {
        let n = reader.read(&mut buf[read..]).await?;
        if n == 0 {
            break;
        }
        read += n;
    }

    if read == BLOCK_HEADER_SIZE {
        let block_id = i64::from_le_bytes(buf[0..8].try_into().unwrap());
        let block_size = u64::from_le_bytes(buf[8..16].try_into().unwrap());
        let is_last = i8::from_le_bytes([buf[16]]);
        return Ok(Some(BlockHeader {
            block_id,
            block_size,
            is_last,
            raw: buf,
        }));
    }

    if read == COMPLETION_MSG_SIZE {
        let mut completion_buf = [0u8; COMPLETION_MSG_SIZE];
        completion_buf.copy_from_slice(&buf[..COMPLETION_MSG_SIZE]);
        if decode_completion_msg(&completion_buf) {
            info!("[receiver.rs] Detected completion message");
            return Ok(None);
        }
    }

    anyhow::bail!("incomplete block header or missing completion message");
}

fn find_completion_msg(buf: &[u8]) -> Option<usize> {
    if buf.len() < COMPLETION_MSG_SIZE {
        return None;
    }
    for i in 0..=buf.len() - COMPLETION_MSG_SIZE {
        if buf[i..i + COMPLETION_MSG_SIZE] == COMPLETION_MSG {
            return Some(i);
        }
    }
    None
}

async fn discard_and_write(
    payload: &[u8],
    skip_bytes: &mut u64,
    to_stdout: bool,
    file: &mut Option<File>,
    stdout: &mut tokio::io::Stdout,
) -> Result<(usize, usize)> {
    if payload.is_empty() {
        return Ok((0, 0));
    }

    let mut discarded = 0usize;
    let mut written = 0usize;

    if *skip_bytes > 0 {
        let skip_now = (*skip_bytes as usize).min(payload.len());
        *skip_bytes -= skip_now as u64;
        discarded = skip_now;
    }

    let write_start = discarded;
    if write_start < payload.len() {
        let data = &payload[write_start..];
        if to_stdout {
            stdout.write_all(data).await?;
        } else if let Some(f) = file {
            f.write_all(data).await?;
        }
        written = data.len();
    }

    Ok((discarded, written))
}
