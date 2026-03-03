use anyhow::Result;
use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt};

use crate::metadata::{
    decode_completion_msg, decode_handshake, LambdaMetadata, StreamMode, COMPLETION_MSG_SIZE,
    HANDSHAKE_SIZE,
};

const BLOCK_HEADER_SIZE: usize = 24;

struct BlockHeader {
    #[allow(dead_code)]
    block_id: i64,
    block_size: u64,
    #[allow(dead_code)]
    is_last: i8,
}

#[derive(Debug, Clone)]
pub enum ReceiverProgress {
    Chunk {
        success: bool,
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

pub async fn read_metadata<R>(reader: &mut R) -> Result<LambdaMetadata>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0u8; HANDSHAKE_SIZE];
    reader.read_exact(&mut buf).await?;
    decode_handshake(&buf)
}

pub async fn read_payload_to_fifo<R>(
    reader: &mut R,
    fifo_name: &str,
    mode: StreamMode,
) -> Result<ReceiverProgress>
where
    R: AsyncRead + Unpin,
{
    match mode {
        StreamMode::RawBytes => read_raw(reader, fifo_name).await,
        StreamMode::Chunked => read_chunked(reader, fifo_name).await,
    }
}

async fn read_raw<R>(reader: &mut R, fifo_name: &str) -> Result<ReceiverProgress>
where
    R: AsyncRead + Unpin,
{
    let to_stdout = fifo_name == "-";
    let mut file = if to_stdout {
        None
    } else {
        Some(File::open(fifo_name).await?)
    };
    let mut stdout = tokio::io::stdout();
    let mut buf = [0u8; 8192];
    let mut recv_bytes = 0u64;

    loop {
        match reader.read(&mut buf).await {
            Ok(0) => {
                return Ok(ReceiverProgress::Raw {
                    success: true,
                    num_of_recv_bytes: recv_bytes,
                });
            }
            Ok(n) => {
                if to_stdout {
                    stdout.write_all(&buf[..n]).await?;
                } else if let Some(f) = &mut file {
                    f.write_all(&buf[..n]).await?;
                }
                recv_bytes += n as u64;
            }
            Err(_) => {
                return Ok(ReceiverProgress::Raw {
                    success: false,
                    num_of_recv_bytes: recv_bytes,
                });
            }
        }
    }
}

async fn read_chunked<R>(reader: &mut R, fifo_name: &str) -> Result<ReceiverProgress>
where
    R: AsyncRead + Unpin,
{
    let to_stdout = fifo_name == "-";
    let mut file = if to_stdout {
        None
    } else {
        Some(File::create(fifo_name).await?)
    };
    let mut stdout = tokio::io::stdout();
    let mut buf = [0u8; 8192];
    let mut completed_chunks = 0u64;

    loop {
        let header = match read_block_header_or_completion(reader).await {
            Ok(Some(header)) => header,
            Ok(None) => {
                return Ok(ReceiverProgress::Chunk {
                    success: true,
                    num_of_completed_chunks: completed_chunks,
                });
            }
            Err(_) => {
                return Ok(ReceiverProgress::Chunk {
                    success: false,
                    num_of_completed_chunks: completed_chunks,
                });
            }
        };

        let mut remaining = header.block_size;
        while remaining > 0 {
            let to_read = (buf.len() as u64).min(remaining) as usize;
            let n = match reader.read(&mut buf[..to_read]).await {
                Ok(0) => {
                    return Ok(ReceiverProgress::Chunk {
                        success: false,
                        num_of_completed_chunks: completed_chunks,
                    });
                }
                Ok(n) => n,
                Err(_) => {
                    return Ok(ReceiverProgress::Chunk {
                        success: false,
                        num_of_completed_chunks: completed_chunks,
                    });
                }
            };

            if to_stdout {
                stdout.write_all(&buf[..n]).await?;
            } else if let Some(f) = &mut file {
                f.write_all(&buf[..n]).await?;
            }
            remaining = remaining.saturating_sub(n as u64);
        }

        completed_chunks += 1;
    }
}

async fn read_block_header_or_completion<R>(reader: &mut R) -> Result<Option<BlockHeader>>
where
    R: AsyncRead + Unpin,
{
    const BLOCK_HEADER_SIZE: usize = 24;
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
