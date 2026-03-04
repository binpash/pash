use anyhow::Result;
use tokio::fs::File;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};
use tracing::info;

use crate::metadata::{encode_handshake, LambdaMetadata, COMPLETION_MSG};

fn short_rdv_key(rdv_key: &str) -> &str {
    rdv_key.get(..6).unwrap_or(rdv_key)
}

pub async fn write_metadata<W>(writer: &mut W, metadata: &LambdaMetadata) -> Result<()>
where
    W: AsyncWrite + Unpin,
{
    let buf = encode_handshake(metadata);
    writer.write_all(&buf).await?;
    Ok(())
}

pub async fn write_metadata_with_rdv<W>(
    writer: &mut W,
    metadata: &LambdaMetadata,
    rdv_key: &str,
) -> Result<()>
where
    W: AsyncWrite + Unpin,
{
    let buf = encode_handshake(metadata);
    let rdv_key = short_rdv_key(rdv_key);
    info!(
        is_stateless = metadata.is_stateless,
        chunk_start_idx = metadata.chunk_start_idx,
        "[sender.rs][{}] Sending metadata",
        rdv_key
    );
    writer.write_all(&buf).await?;
    Ok(())
}

pub async fn write_fifo_payload<W>(writer: &mut W, fifo_name: &str, rdv_key: &str) -> Result<u64>
where
    W: AsyncWrite + Unpin,
{
    let rdv_key = short_rdv_key(rdv_key);
    let copied = if fifo_name == "-" {
        info!("[sender.rs][{}] Writing payload from stdin", rdv_key);
        let mut stdin = io::stdin();
        let copied = io::copy(&mut stdin, writer).await?;
        info!(
            "[sender.rs][{}] Finished writing payload from stdin, bytes copied: {}",
            rdv_key, copied
        );
        copied
    } else {
        let mut file = File::open(fifo_name).await?;
        info!("[sender.rs][{}] Writing payload from fifo {}", rdv_key, fifo_name);
        let copied = io::copy(&mut file, writer).await?;
        info!(
            "[sender.rs][{}] Finished writing payload from fifo {}, bytes copied: {}",
            rdv_key, fifo_name, copied
        );
        copied
    };

    // Send completion message
    info!("[sender.rs][{}] Sending completion message", rdv_key);
    writer.write_all(&COMPLETION_MSG).await?;
    writer.flush().await?;
    info!("[sender.rs][{}] Sent completion message", rdv_key);

    Ok(copied)
}
