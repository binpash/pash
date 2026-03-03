use anyhow::Result;
use tokio::fs::File;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};
use tracing::info;

use crate::metadata::{encode_handshake, LambdaMetadata, COMPLETION_MSG};

pub async fn write_metadata<W>(writer: &mut W, metadata: &LambdaMetadata) -> Result<()>
where
    W: AsyncWrite + Unpin,
{
    let buf = encode_handshake(metadata);
    info!(
        is_stateless = metadata.is_stateless,
        chunk_start_id = metadata.chunk_start_id,
        "[sender.rs] Sending metadata"
    );
    writer.write_all(&buf).await?;
    Ok(())
}

pub async fn write_fifo_payload<W>(writer: &mut W, fifo_name: &str) -> Result<u64>
where
    W: AsyncWrite + Unpin,
{
    let copied = if fifo_name == "-" {
        info!("[sender.rs] Writing payload from stdin");
        let mut stdin = io::stdin();
        let copied = io::copy(&mut stdin, writer).await?;
        info!("[sender.rs] Finished writing payload from stdin, bytes copied: {}", copied);
        copied
    } else {
        let mut file = File::open(fifo_name).await?;
        info!("[sender.rs] Writing payload from fifo {}", fifo_name);
        let copied = io::copy(&mut file, writer).await?;
        info!("[sender.rs] Finished writing payload from fifo {}, bytes copied: {}", fifo_name, copied);
        copied
    };

    // Send completion message
    info!("[sender.rs] Sending completion message");
    writer.write_all(&COMPLETION_MSG).await?;
    writer.flush().await?;
    info!("[sender.rs] Sent completion message");

    Ok(copied)
}
