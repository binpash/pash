use anyhow::{anyhow, Result};
use tracing::{info};

use crate::db_helper::{create_rdv_table_if_not_exists, make_db_client};
use crate::holepunch::PashCtx;
use crate::lambda_helper::invoke_recovery_lambda;
use crate::metadata::{LambdaMetadata, StreamMode};
use crate::receiver::{read_metadata, ChunkReader, RawReader, ReceiverProgress};
use crate::sender::{write_fifo_payload, write_metadata};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndpointMode {
    Send,
    Recv,
}

#[derive(Debug, Clone)]
pub struct FtExecutor {
    mode: EndpointMode,
    rdv_key: String,
    me: String,
    peer: String,
    fifo_name: String,
}

impl FtExecutor {
    pub fn new(arg: &str) -> Result<Self> {
        let parts: Vec<&str> = arg.split('*').collect();
        if parts.len() < 5 {
            return Err(anyhow!(
                "invalid endpoint arg: expected at least 5 fields, got {}",
                parts.len()
            ));
        }

        let mode = if parts[0].starts_with("send") {
            EndpointMode::Send
        } else if parts[0].starts_with("recv") {
            EndpointMode::Recv
        } else {
            return Err(anyhow!(
                "invalid mode '{}': expected send*... or recv*...",
                parts[0]
            ));
        };

        Ok(Self {
            mode,
            rdv_key: parts[1].to_string(),
            me: parts[2].to_string(),
            peer: parts[3].to_string(),
            fifo_name: parts[4].to_string(),
        })
    }

    async fn connect_to_peer(&self) -> tokio::net::TcpStream {
        let mut ctx = PashCtx::new(&self.me, &self.rdv_key).await;
        ctx.connect(&self.peer).await
    }

    fn completed_chunks(progress: &ReceiverProgress) -> u64 {
        match progress {
            ReceiverProgress::Chunk {
                num_of_completed_chunks,
                ..
            } => *num_of_completed_chunks,
            // Raw mode does not use chunk-based resume.
            ReceiverProgress::Raw { .. } => 0,
        }
    }

    async fn send_attempt<W>(&self, writer: &mut W) -> Result<LambdaMetadata>
    where
        W: tokio::io::AsyncWrite + Unpin,
    {
        let metadata = LambdaMetadata::from_env();
        write_metadata(writer, &metadata).await?;
        // we do not handle recovery at the sender side
        write_fifo_payload(writer, &self.fifo_name).await?;
        Ok(metadata)
    }

    async fn recv_attempt<R>(
        &self,
        reader: &mut R,
        raw_reader: &mut Option<RawReader>,
        chunk_reader: &mut Option<ChunkReader>,
    ) -> Result<(ReceiverProgress, LambdaMetadata)>
    where
        R: tokio::io::AsyncRead + Unpin,
    {
        let metadata = read_metadata(reader).await?;
        let recv_mode = if metadata.is_stateless {
            StreamMode::Chunked
        } else {
            StreamMode::RawBytes
        };
        info!(
            is_stateless = metadata.is_stateless,
            recv_mode = ?recv_mode,
            "[recovery.rs] Selected recv mode"
        );
        let progress = match recv_mode {
            StreamMode::RawBytes => {
                if raw_reader.is_none() {
                    *raw_reader = Some(RawReader::new(&self.fifo_name).await?);
                }
                raw_reader.as_mut().unwrap().read_from(reader).await
            }
            StreamMode::Chunked => {
                if chunk_reader.is_none() {
                    *chunk_reader = Some(ChunkReader::new(&self.fifo_name).await?);
                }
                chunk_reader.as_mut().unwrap().read_from(reader).await
            }
        }?;
        Ok((progress, metadata))
    }

    pub async fn execute(&self) -> Result<Vec<ReceiverProgress>> {
        let mut progress_history = Vec::new();
        let mut attempt: u64 = 0;
        let mut raw_reader: Option<RawReader> = None;
        let mut chunk_reader: Option<ChunkReader> = None;
        let client = make_db_client().await;
        create_rdv_table_if_not_exists(&client).await;
        let lambda_client = {
            let cfg = aws_config::defaults(aws_config::BehaviorVersion::latest())
                .load()
                .await;
            Some(aws_sdk_lambda::Client::new(&cfg))
        };

        // Recovery loop: connect -> handshake -> read/write per attempt.
        loop {
            attempt += 1;
            info!(
                attempt,
                mode = ?self.mode,
                me = %self.me,
                peer = %self.peer,
                rdv_key = %self.rdv_key,
                "[recovery.rs] Attempt start"
            );
            let attempt_result: Result<Option<(ReceiverProgress, Option<LambdaMetadata>)>> = async {
                // 1) connect
                let stream = self.connect_to_peer().await;
                let (mut rd, mut wr) = stream.into_split();
                info!(attempt, me = %self.me, peer = %self.peer, "[recovery.rs] Connected to peer");

                // 2) send/recv
                match self.mode {
                    EndpointMode::Send => {
                        info!(attempt, "[recovery.rs] Sending payload");
                        let metadata = self.send_attempt(&mut wr).await?;
                        info!(attempt, %metadata, "[recovery.rs] Sent metadata");
                        info!(attempt, "[recovery.rs] Sent payload");
                        Ok(None)
                    }
                    EndpointMode::Recv => {
                        info!(attempt, "[recovery.rs] Receiving payload");
                        let (progress, metadata) =
                            self.recv_attempt(&mut rd, &mut raw_reader, &mut chunk_reader).await?;
                        info!(attempt, metadata = %metadata, "[recovery.rs] Received metadata");
                        let (success, num, class) = match &progress {
                            ReceiverProgress::Chunk {
                                success,
                                num_of_completed_chunks: num,
                                ..
                            } => (*success, *num, "chunk"),

                            ReceiverProgress::Raw {
                                success,
                                num_of_recv_bytes: num,
                            } => (*success, *num, "raw"),
                        };

                        info!(
                            attempt,
                            success,
                            num,
                            class,
                            "[recovery.rs] Received payload"
                        );
                        Ok(Some((progress, Some(metadata))))
                    }
                }
            }
            .await;

            match attempt_result {
                Ok(Some((progress, metadata))) => {
                    let is_success = progress.success();
                    progress_history.push(progress.clone());
                    if is_success {
                        info!(attempt, "[recovery.rs] Final recv succeeded");
                        return Ok(progress_history);
                    }
                    if let (Some(lambda_client), Some(metadata)) =
                        (lambda_client.as_ref(), metadata.as_ref())
                    {
                        let num_of_completed_chunks = Self::completed_chunks(&progress);
                        if let Err(err) = invoke_recovery_lambda(
                            lambda_client,
                            "lambda",
                            metadata,
                            num_of_completed_chunks,
                        )
                        .await
                        {
                            info!(attempt, error = %err, "[recovery.rs] Recovery lambda invocation failed");
                        } else {
                            info!(attempt, function = "lambda", "[recovery.rs] Recovery lambda invoked");
                        }
                    }
                    info!(attempt, "[recovery.rs] Recv incomplete, retrying");
                }
                Ok(None) => {
                    info!(attempt, "[recovery.rs] Final send succeeded");
                    return Ok(progress_history);
                }
                Err(err) => {
                    if self.mode == EndpointMode::Send {
                        info!(attempt, error = %err, "[recovery.rs] Sender attempt failed, not retrying");
                        return Err(err);
                    }
                    info!(attempt, error = %err, "[recovery.rs] Attempt failed, retrying");
                    continue;
                }
            }
        }
    }
}
