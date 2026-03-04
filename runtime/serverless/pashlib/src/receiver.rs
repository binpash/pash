use anyhow::Result;
use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt};
use tracing::info;

use crate::metadata::{
    decode_completion_msg, decode_handshake, LambdaMetadata, COMPLETION_MSG,
    COMPLETION_MSG_SIZE, HANDSHAKE_SIZE,
};

const BLOCK_HEADER_SIZE: usize = 24;

fn short_rdv_key(rdv_key: &str) -> &str {
    rdv_key.get(..6).unwrap_or(rdv_key)
}

struct BlockHeader {
    #[allow(dead_code)]
    block_id: i64,
    block_size: u64,
    is_last: i8,
    raw: [u8; BLOCK_HEADER_SIZE],
}

#[derive(Debug, Clone)]
pub enum ReceiverProgress {
    Chunk {
        success: bool,
        // Number of fully completed chunks in this attempt.
        // A chunk is complete when we see a block with is_last != 0.
        num_of_completed_chunks: u64,
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

pub async fn read_metadata<R>(reader: &mut R, rdv_key: &str) -> Result<LambdaMetadata>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0u8; HANDSHAKE_SIZE];
    reader.read_exact(&mut buf).await?;
    let metadata = decode_handshake(&buf)?;
    info!(
        is_stateless = metadata.is_stateless,
        chunk_start_idx = metadata.chunk_start_idx,
        "[receiver.rs][{}] Decoded metadata",
        short_rdv_key(rdv_key)
    );
    Ok(metadata)
}

pub struct RawReader {
    rdv_key: String,
    to_stdout: bool,
    fifo_name: String,
    file: Option<File>,
    stdout: tokio::io::Stdout,
    total_received: u64,
    carry: Vec<u8>,
}

impl RawReader {
    pub async fn new(fifo_name: &str, rdv_key: &str) -> Result<Self> {
        let to_stdout = fifo_name == "-";
        let file = if to_stdout {
            None
        } else {
            Some(File::create(fifo_name).await?)
        };
        Ok(Self {
            rdv_key: short_rdv_key(rdv_key).to_string(),
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
        info!(
            output = if self.to_stdout { "stdout" } else { self.fifo_name.as_str() },
            prev_forwarded_bytes = self.total_received,
            "[receiver.rs][{}] RawReader attempt start",
            self.rdv_key
        );
        let mut buf = [0u8; 8192];
        // Number of bytes forwarded to downstream in this attempt.
        let mut recv_bytes = 0u64;
        // Retry dedup: bytes already forwarded in previous attempts.
        let mut skip_remaining = self.total_received;
        let mut skipped_this_attempt = 0u64;
        self.carry.clear();
        loop {
            match reader.read(&mut buf).await {
                Ok(0) => {
                    self.total_received = self.total_received.saturating_add(recv_bytes);
                    info!(
                        skipped_this_attempt,
                        recv_bytes,
                        total_received = self.total_received,
                        "[receiver.rs][{}] RawReader attempt end (incomplete)",
                        self.rdv_key
                    );
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
                            let payload = &combined[..pos];
                            let (discarded, written) = discard_and_write(
                                payload,
                                &mut skip_remaining,
                                self.to_stdout,
                                &mut self.file,
                                &mut self.stdout,
                            )
                            .await?;
                            skipped_this_attempt = skipped_this_attempt.saturating_add(discarded as u64);
                            recv_bytes += written as u64;
                            self.total_received = self.total_received.saturating_add(recv_bytes);
                            info!(
                                skipped_this_attempt,
                                recv_bytes,
                                total_received = self.total_received,
                                "[receiver.rs][{}] RawReader attempt end (success)",
                                self.rdv_key
                            );
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
                                let (discarded, written) = discard_and_write(
                                    payload,
                                    &mut skip_remaining,
                                    self.to_stdout,
                                    &mut self.file,
                                    &mut self.stdout,
                                )
                                .await?;
                                skipped_this_attempt = skipped_this_attempt.saturating_add(discarded as u64);
                                recv_bytes += written as u64;
                            }
                            self.carry.clear();
                            self.carry.extend_from_slice(&combined[write_len..]);
                        } else {
                            self.carry = combined;
                        }
                    } else {
                        if let Some(pos) = find_completion_msg(&buf[..n]) {
                            let payload = &buf[..pos];
                            let (discarded, written) = discard_and_write(
                                payload,
                                &mut skip_remaining,
                                self.to_stdout,
                                &mut self.file,
                                &mut self.stdout,
                            )
                            .await?;
                            skipped_this_attempt = skipped_this_attempt.saturating_add(discarded as u64);
                            recv_bytes += written as u64;
                            self.total_received = self.total_received.saturating_add(recv_bytes);
                            info!(
                                skipped_this_attempt,
                                recv_bytes,
                                total_received = self.total_received,
                                "[receiver.rs][{}] RawReader attempt end (success)",
                                self.rdv_key
                            );
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
                                let (discarded, written) = discard_and_write(
                                    payload,
                                    &mut skip_remaining,
                                    self.to_stdout,
                                    &mut self.file,
                                    &mut self.stdout,
                                )
                                .await?;
                                skipped_this_attempt = skipped_this_attempt.saturating_add(discarded as u64);
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
                    info!(
                        skipped_this_attempt,
                        recv_bytes,
                        total_received = self.total_received,
                        "[receiver.rs][{}] RawReader attempt end (read error)",
                        self.rdv_key
                    );
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
    rdv_key: String,
    to_stdout: bool,
    file: Option<File>,
    stdout: tokio::io::Stdout,
    // Total chunks completed across attempts.
    completed_chunks: u64,
    // - "block" = a single r_merge header + payload.
    // - "chunk" = one or more blocks within the same logical S3 chunk,
    //             consisting of a series block with the same block_id, and
    //             terminated when a block has is_last != 0.
    //
    // For the current (incomplete) chunk, how many full blocks were forwarded.
    // On retry, we skip these blocks.
    partial_blocks_forwarded: u64,
    // For the current block within the current chunk, how many bytes were forwarded.
    // On retry, we discard these bytes from the payload.
    partial_block_bytes_forwarded: u64,
    // Whether the current block's header has already been forwarded.
    partial_header_written: bool,
}

impl ChunkReader {
    pub async fn new(fifo_name: &str, rdv_key: &str) -> Result<Self> {
        let to_stdout = fifo_name == "-";
        let file = if to_stdout {
            None
        } else {
            Some(File::create(fifo_name).await?)
        };
        Ok(Self {
            rdv_key: short_rdv_key(rdv_key).to_string(),
            to_stdout,
            file,
            stdout: tokio::io::stdout(),
            completed_chunks: 0,
            partial_blocks_forwarded: 0,
            partial_block_bytes_forwarded: 0,
            partial_header_written: false,
        })
    }

    pub async fn read_from<R>(&mut self, reader: &mut R) -> Result<ReceiverProgress>
    where
        R: AsyncRead + Unpin,
    {
        let mut buf = [0u8; 8192];
        // Retry dedup state:
        // - We do not skip already completed chunks here. Recovery restarts from
        //   the next chunk using chunk_start_idx.
        // - We only skip duplicated blocks/bytes in the current incomplete chunk.
        let mut attempt_completed = 0u64;
        let mut skip_blocks_remaining = self.partial_blocks_forwarded;
        let mut skip_bytes_remaining = self.partial_block_bytes_forwarded;
        let mut header_already_written = self.partial_header_written;

        info!(
            prev_completed_chunks = self.completed_chunks,
            skip_blocks_remaining,
            skip_bytes_remaining,
            "[receiver.rs][{}] ChunkReader attempt start",
            self.rdv_key
        );

        loop {
            let header = match read_block_header_or_completion(reader).await {
                Ok(Some(header)) => header,
                Ok(None) => {
                    info!(
                        attempt_completed,
                        total_completed_chunks = self.completed_chunks,
                        "[receiver.rs][{}] ChunkReader attempt end (success)",
                        self.rdv_key
                    );
                    return Ok(ReceiverProgress::Chunk {
                        success: true,
                        num_of_completed_chunks: attempt_completed,
                    });
                }
                Err(_) => {
                    info!(
                        attempt_completed,
                        total_completed_chunks = self.completed_chunks,
                        "[receiver.rs][{}] ChunkReader attempt end (header read error)",
                        self.rdv_key
                    );
                    return Ok(ReceiverProgress::Chunk {
                        success: false,
                        num_of_completed_chunks: attempt_completed,
                    });
                }
            };

            // Skip full blocks already forwarded in the current (incomplete) chunk.
            if skip_blocks_remaining > 0 {
                let mut remaining = header.block_size;
                while remaining > 0 {
                    let to_read = (buf.len() as u64).min(remaining) as usize;
                    let n = match reader.read(&mut buf[..to_read]).await {
                        Ok(0) => {
                            info!(
                                attempt_completed,
                                total_completed_chunks = self.completed_chunks,
                                "[receiver.rs][{}] ChunkReader attempt end (incomplete while skipping block)",
                                self.rdv_key
                            );
                            return Ok(ReceiverProgress::Chunk {
                                success: false,
                                num_of_completed_chunks: attempt_completed,
                            });
                        }
                        Ok(n) => n,
                        Err(_) => {
                            info!(
                                attempt_completed,
                                total_completed_chunks = self.completed_chunks,
                                "[receiver.rs][{}] ChunkReader attempt end (read error while skipping block)",
                                self.rdv_key
                            );
                            return Ok(ReceiverProgress::Chunk {
                                success: false,
                                num_of_completed_chunks: attempt_completed,
                            });
                        }
                    };
                    remaining = remaining.saturating_sub(n as u64);
                }
                skip_blocks_remaining -= 1;
                header_already_written = false;
                continue;
            }

            // Forward this block (and discard any previously forwarded bytes within it).
            if !header_already_written {
                if self.to_stdout {
                    self.stdout.write_all(&header.raw).await?;
                } else if let Some(f) = &mut self.file {
                    f.write_all(&header.raw).await?;
                }
                header_already_written = true;
            }

            let mut remaining = header.block_size;
            let mut current_written: u64 = 0;
            while remaining > 0 {
                let to_read = (buf.len() as u64).min(remaining) as usize;
                let n = match reader.read(&mut buf[..to_read]).await {
                    Ok(0) => {
                        self.partial_block_bytes_forwarded =
                            self.partial_block_bytes_forwarded.saturating_add(current_written);
                        self.partial_header_written = header_already_written;
                        info!(
                            attempt_completed,
                            total_completed_chunks = self.completed_chunks,
                            partial_blocks_forwarded = self.partial_blocks_forwarded,
                            partial_block_bytes_forwarded = self.partial_block_bytes_forwarded,
                            "[receiver.rs][{}] ChunkReader attempt end (incomplete while forwarding)",
                            self.rdv_key
                        );
                        return Ok(ReceiverProgress::Chunk {
                            success: false,
                            num_of_completed_chunks: attempt_completed,
                        });
                    }
                    Ok(n) => n,
                    Err(_) => {
                        self.partial_block_bytes_forwarded =
                            self.partial_block_bytes_forwarded.saturating_add(current_written);
                        self.partial_header_written = header_already_written;
                        info!(
                            attempt_completed,
                            total_completed_chunks = self.completed_chunks,
                            partial_blocks_forwarded = self.partial_blocks_forwarded,
                            partial_block_bytes_forwarded = self.partial_block_bytes_forwarded,
                            "[receiver.rs][{}] ChunkReader attempt end (read error while forwarding)",
                            self.rdv_key
                        );
                        return Ok(ReceiverProgress::Chunk {
                            success: false,
                            num_of_completed_chunks: attempt_completed,
                        });
                    }
                };

                let mut start = 0usize;
                if skip_bytes_remaining > 0 {
                    let discard = (skip_bytes_remaining as usize).min(n);
                    skip_bytes_remaining -= discard as u64;
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

                remaining = remaining.saturating_sub(n as u64);
            }

            self.partial_blocks_forwarded = self.partial_blocks_forwarded.saturating_add(1);
            self.partial_block_bytes_forwarded = 0;
            skip_bytes_remaining = 0;
            self.partial_header_written = false;
            header_already_written = false;

            if header.is_last != 0 {
                attempt_completed += 1;
                self.completed_chunks = self.completed_chunks.saturating_add(1);
                self.partial_blocks_forwarded = 0;
                self.partial_block_bytes_forwarded = 0;
            }
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
        let is_last = buf[16] as i8;
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
