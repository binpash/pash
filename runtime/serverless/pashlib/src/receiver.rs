use anyhow::Result;
use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt};

use crate::metadata::{decode_handshake, LambdaMetadata, StreamMode, HANDSHAKE_SIZE};

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
        match reader.read(&mut buf).await {
            Ok(0) => {
                return Ok(ReceiverProgress::Chunk {
                    success: true,
                    num_of_completed_chunks: completed_chunks,
                });
            }
            Ok(n) => {
                if to_stdout {
                    stdout.write_all(&buf[..n]).await?;
                } else if let Some(f) = &mut file {
                    f.write_all(&buf[..n]).await?;
                }
                completed_chunks += 1;
            }
            Err(_) => {
                return Ok(ReceiverProgress::Chunk {
                    success: false,
                    num_of_completed_chunks: completed_chunks,
                });
            }
        }
    }
}
