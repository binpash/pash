use anyhow::{anyhow, Result};
use tracing::{info, warn};

use crate::holepunch::PashCtx;
use crate::metadata::{LambdaMetadata, StreamMode};
use crate::receiver::{read_metadata, read_payload_to_fifo, ReceiverProgress};
use crate::sender::{write_fifo_payload, write_metadata};
use crate::db_helper::{make_db_client, create_rdv_table_if_not_exists};

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

    pub async fn execute(&self) -> Result<Vec<ReceiverProgress>> {
        let mut progress_history = Vec::new();
        let mut attempt: u64 = 0;
        let client = make_db_client().await;
        create_rdv_table_if_not_exists(&client).await;

        // Recovery loop: connect -> handshake -> read/write per attempt.
        loop {
            attempt += 1;
            info!(
                attempt,
                mode = ?self.mode,
                me = %self.me,
                peer = %self.peer,
                rdv_key = %self.rdv_key,
                "[recovery.rs] attempt start"
            );
            let attempt_result: Result<Option<ReceiverProgress>> = async {
                // 1) connect
                let stream = self.connect_to_peer().await;
                let (mut rd, mut wr) = stream.into_split();
                info!(attempt, me = %self.me, peer = %self.peer, "[recovery.rs] connected");

                // 2) handshake
                let metadata = match self.mode {
                    EndpointMode::Send => {
                        let metadata = LambdaMetadata::from_env();
                        write_metadata(&mut wr, &metadata).await?;
                        info!(
                            attempt,
                            job = %metadata.leash_job_id,
                            chunk_start_id = metadata.chunk_start_id,
                            is_stateless = metadata.is_stateless,
                            "[recovery.rs] sent metadata"
                        );
                        metadata
                    }
                    EndpointMode::Recv => {
                        let lambda_metadata = read_metadata(&mut rd).await?;
                        info!(
                            attempt,
                            job = %lambda_metadata.leash_job_id,
                            chunk_start_id = lambda_metadata.chunk_start_id,
                            is_stateless = lambda_metadata.is_stateless,
                            "[recovery.rs] received metadata"
                        );
                        lambda_metadata
                    }
                };
                
                let recv_mode = if metadata.is_stateless {
                    StreamMode::RawBytes
                } else {
                    StreamMode::Chunked
                };

                // 3) read or write based on mode
                match self.mode {
                    EndpointMode::Send => {
                        // we do not handle recovery at the sender side
                        write_fifo_payload(&mut wr, &self.fifo_name).await?;
                        info!(attempt, fifo = %self.fifo_name, "[recovery.rs] sent payload");
                        Ok(None)
                    }
                    EndpointMode::Recv => {
                        let progress = read_payload_to_fifo(
                            &mut rd,
                            &self.fifo_name,
                            recv_mode,
                        )
                        .await?;
                        let (success, num, class) = match &progress {
                            ReceiverProgress::Chunk {
                                success,
                                num_of_completed_chunks: num,
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
                            fifo = %self.fifo_name,
                            "[recovery.rs] received payload"
                        );
                        Ok(Some(progress))
                    }
                }
            }.await;

            match attempt_result {
                Ok(Some(progress)) => {
                    let is_success = progress.success();
                    progress_history.push(progress);
                    if is_success {
                        info!(attempt, "[recovery.rs] recv succeeded");
                        return Ok(progress_history);
                    }
                    warn!(attempt, "[recovery.rs] recv incomplete, retrying");
                }
                Ok(None) => {
                    info!(attempt, "[recovery.rs] send succeeded");
                    return Ok(progress_history);
                }
                Err(err) => {
                    warn!(attempt, error = %err, "[recovery.rs] attempt failed, retrying");
                    continue;
                }
            }
        }
    }
}
