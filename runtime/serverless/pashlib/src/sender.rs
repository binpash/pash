use anyhow::Result;
use tokio::fs::OpenOptions;
use tokio::fs::File;
use tokio::io::{self, AsyncReadExt, AsyncWrite, AsyncWriteExt};
use tracing::info;
use std::process::Command;
use std::os::unix::fs::FileTypeExt;

use crate::holepunch::PashCtx;
use crate::metadata::{
    encode_handshake, read_resumability_ack, resumability_disabled, write_resumability_request,
    LambdaMetadata, ResumeAck, ResumeRequest, COMPLETION_MSG,
};

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
    // Stream stdin or fifo into the data channel and append a completion marker.
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

pub async fn monitor_resumability_and_forward(
    me: &str,
    peer: &str,
    rdv_key: &str,
    metadata: &LambdaMetadata,
    mut shutdown: tokio::sync::watch::Receiver<bool>,
) -> Result<()> {
    // Forward local resumability requests to the receiver over the control plane.
    let safe = if metadata.script_id.is_empty() {
        "unknown"
    } else {
        metadata.script_id.as_str()
    };
    let req_path = format!("/tmp/pash_resume_{}.req", safe);
    let ack_path = format!("/tmp/pash_resume_{}.ack", safe);
    let ctrl_key = format!("ctrl::{}", rdv_key);
    let mut ctx = PashCtx::new(me, &ctrl_key).await;
    let mut ctrl_stream = ctx.connect(peer).await;
    info!(
        "[sender.rs][{}] Resumability monitor connected to control channel",
        short_rdv_key(rdv_key)
    );

    for p in [&req_path, &ack_path] {
        if !std::path::Path::new(p).exists() {
            let status = Command::new("mkfifo").arg(p).status()?;
            if !status.success() {
                // If another process created it in the meantime, accept it.
                if let Ok(meta) = std::fs::metadata(p) {
                    if meta.file_type().is_fifo() {
                        continue;
                    }
                }
                info!(error = ?status, path = p, "[sender.rs][{}] Failed to create fifo for resumability monitor", short_rdv_key(rdv_key));
                return Err(anyhow::anyhow!("mkfifo failed for {}", p));
            }
        }
    }

    // Open FIFO in read-only mode; s3-reader keeps the writer end open.
    let mut req_file = OpenOptions::new().read(true).open(&req_path).await?;
    info!(
        "[sender.rs][{}] Resumability monitor started, watching for local resumability requests at {}",
        short_rdv_key(rdv_key),
        req_path
    );
    let mut buf = [0u8; 4];
    let n = tokio::select! {
        _ = shutdown.changed() => {
            info!(
                "[sender.rs][{}] Resumability monitor received shutdown signal, exiting",
                short_rdv_key(rdv_key)
            );
            return Ok(());
        }
        res = req_file.read_exact(&mut buf) => {
            res?;
            4usize
        },
    };
    if n == 0 {
        return Ok(());
    }
    let next_chunk_start_idx = u32::from_be_bytes(buf);
    info!(
        next_chunk_start_idx,
        "[sender.rs][{}] Local resumability request received",
        short_rdv_key(rdv_key)
    );

    // Relay the resumability request to the receiver over the control channel.
    let req = ResumeRequest {
        next_chunk_start_idx,
    };
    write_resumability_request(&mut ctrl_stream, &req).await?;
    info!(
        next_chunk_start_idx = req.next_chunk_start_idx,
        "[sender.rs][{}] Forwarded resumability request to receiver",
        short_rdv_key(rdv_key)
    );
    let ack: ResumeAck = read_resumability_ack(&mut ctrl_stream).await?;
    info!(
        accepted = ack.accepted,
        "[sender.rs][{}] Received resumability ack from receiver",
        short_rdv_key(rdv_key)
    );

    if std::path::Path::new(&ack_path).exists() {
        // Write ack back so the local s3-reader can close its downstream and finish.
        let mut ack_file = OpenOptions::new().write(true).open(&ack_path).await?;
        let ack_payload = serde_json::to_string(&ack)?;
        ack_file.write_all(ack_payload.as_bytes()).await?;
        ack_file.write_all(b"\n").await?;
        ack_file.flush().await?;
    }
    info!(
        "[sender.rs][{}] Resumability monitor finished handling request",
        short_rdv_key(rdv_key)
    );
    Ok(())
}

pub async fn send(
    me: &str,
    peer: &str,
    rdv_key: &str,
    fifo_name: &str,
) -> Result<LambdaMetadata> {

    // Establish data connection, send metadata, and stream payload.
    let metadata = LambdaMetadata::from_env();
    let mut ctx = PashCtx::new(me, rdv_key).await;
    let mut stream = ctx.connect(peer).await;
    write_metadata_with_rdv(&mut stream, &metadata, rdv_key).await?;

    // Always start a control-plane monitor so the s3-reader can trigger resumability at any time.
    let enable_resumability = !resumability_disabled();
    let (mut shutdown_tx, monitor_handle) = if enable_resumability {
        info!(
            "[sender.rs][{}] Resumability enabled, starting control-plane monitor",
            short_rdv_key(rdv_key)
        );
        let me = me.to_string();
        let peer = peer.to_string();
        let rdv_key = rdv_key.to_string();
        let md = metadata.clone();
        let (shutdown_tx, shutdown_rx) = tokio::sync::watch::channel(false);
        let handle = tokio::spawn(async move {
            if let Err(err) =
                monitor_resumability_and_forward(&me, &peer, &rdv_key, &md, shutdown_rx).await
            {
                let suppress = err
                    .downcast_ref::<std::io::Error>()
                    .map(|e| e.kind() == std::io::ErrorKind::UnexpectedEof)
                    .unwrap_or(false);
                if !suppress {
                    info!(
                        error = %err,
                        "[sender.rs] resumability monitor failed"
                    );
                }
            }
        });
        (Some(shutdown_tx), Some(handle))
    } else {
        (None, None)
    };

    write_fifo_payload(&mut stream, fifo_name, rdv_key).await?;

    if let Some(h) = monitor_handle {
        // Signal the control thread to exit cleanly once data is sent.
        if let Some(tx) = shutdown_tx.as_mut() {
            let _ = tx.send(true);
        }
        let _ = h.await;
    }

    Ok(metadata)
}
