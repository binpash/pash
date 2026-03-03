use anyhow::Result;
use tokio::fs::File;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::metadata::{encode_handshake, LambdaMetadata};

pub async fn write_metadata<W>(writer: &mut W, metadata: &LambdaMetadata) -> Result<()>
where
    W: AsyncWrite + Unpin,
{
    let buf = encode_handshake(metadata);
    writer.write_all(&buf).await?;
    Ok(())
}

pub async fn write_fifo_payload<W>(writer: &mut W, fifo_name: &str) -> Result<u64>
where
    W: AsyncWrite + Unpin,
{
    if fifo_name == "-" {
        let mut stdin = io::stdin();
        let copied = io::copy(&mut stdin, writer).await?;
        return Ok(copied);
    }

    let mut file = File::open(fifo_name).await?;
    let copied = io::copy(&mut file, writer).await?;
    Ok(copied)
}
