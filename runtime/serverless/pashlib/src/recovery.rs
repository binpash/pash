use anyhow::{anyhow, Result};
use tokio::sync::mpsc::{unbounded_channel, UnboundedSender};
use tracing::info;

use crate::aggregator::{Aggregator, SortMergeAggregator};
use crate::db_helper::{create_rdv_table_if_not_exists, make_db_client};
use crate::events::{Event, JobSpec};
use crate::lambda_helper::invoke_lambda;
use crate::receiver::{recv, ReceiverProgress};
use crate::sender::send;
use std::fs::OpenOptions;

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
    // Return a short rdv key for concise logs.
    fn rdv_tag(&self) -> &str {
        self.rdv_key.get(..6).unwrap_or(&self.rdv_key)
    }

    // Parse executor args from the pash runtime.
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

    fn completed_chunks(progress: &Option<ReceiverProgress>) -> u64 {
        match progress {
            Some(ReceiverProgress::Chunk {
                total_completed_chunks,
                ..
            }) => *total_completed_chunks,
            _ => 0,
        }
    }

    async fn handle_job_event(
        &self,
        job: JobSpec,
        lambda_client: &Option<aws_sdk_lambda::Client>,
        tx: &UnboundedSender<Event>,
        next_part_id: &mut u64,
    ) {
        // Resumability creates a new part id; recovery reuses the previous part id.
        let part_id = if job.part_id != 0 {
            job.part_id
        } else if job.resume_chunk_start_idx.is_some() {
            *next_part_id = next_part_id.saturating_add(1);
            *next_part_id
        } else {
            if *next_part_id == 0 {
                *next_part_id = 1;
            }
            *next_part_id
        };
        // Executor owns lambda invocation for both recovery and resumability.
        let chunk_start_idx = match job.resume_chunk_start_idx {
            Some(idx) => idx,
            None => {
                if job.metadata.is_stateless {
                    job.metadata
                        .chunk_start_idx
                        .saturating_add(Self::completed_chunks(&job.recovery_progress) as u32)
                } else {
                    0
                }
            }
        };
        if let Some(lambda_client) = lambda_client.as_ref() {
            if let Err(err) = invoke_lambda(
                lambda_client,
                "lambda",
                &job.metadata,
                chunk_start_idx,
                self.rdv_tag(),
            )
            .await
            {
                info!(
                    error = %err,
                    "[recovery.rs][{}] lambda invocation failed",
                    self.rdv_tag()
                );
            }
        }
        // Reader owns data + control channels; executor just schedules it.
        let me = self.me.clone();
        let peer = self.peer.clone();
        let rdv_key = self.rdv_key.clone();
        let fifo_name = self.fifo_name.clone();
        let tx = tx.clone();
        // Spawn the receiver task; it pushes completion/job events back into the queue.
        tokio::spawn(async move {
            if let Err(err) =
                recv(me, peer, rdv_key, fifo_name, tx, part_id, job.recovery_progress)
            .await
            {
                info!(error = %err, "[recovery.rs] recv task failed");
            }
        });
    }

    fn open_fifo_keepalive(&self) -> Result<std::fs::File> {
        // Open read+write so it doesn't block if the downstream hasn't opened yet.
        Ok(OpenOptions::new()
            .read(true)
            .write(true)
            .open(&self.fifo_name)?)
    }

    // Run the executor for send/recv endpoints.
    pub async fn execute(&self) -> Result<Vec<ReceiverProgress>> {
        let mut progress_history = Vec::new();
        let client = make_db_client().await;
        create_rdv_table_if_not_exists(&client).await;
        let lambda_client = {
            let cfg = aws_config::defaults(aws_config::BehaviorVersion::latest())
                .load()
                .await;
            Some(aws_sdk_lambda::Client::new(&cfg))
        };

        if self.mode == EndpointMode::Recv {

            // Keep the FIFO open to avoid premature EOF between reader jobs
            let _fifo_keepalive = self.open_fifo_keepalive()?;

            // Event-driven executor: consume events from readers/control plane
            let (tx, mut rx) = unbounded_channel::<Event>();
            let mut parts: Vec<String> = Vec::new();
            let mut any_temp = false;
            let mut next_part_id: u64 = 1;
            let mut active_jobs: u64 = 0;

            // Initial receiver task for the first lambda attempt.
            active_jobs = active_jobs.saturating_add(1);
            let first_part_id = next_part_id;
            let me = self.me.clone();
            let peer = self.peer.clone();
            let rdv_key = self.rdv_key.clone();
            let fifo_name = self.fifo_name.clone();
            let tx_spawn = tx.clone();
            tokio::spawn(async move {
                if let Err(err) = recv(
                    me,
                    peer,
                    rdv_key,
                    fifo_name,
                    tx_spawn,
                    first_part_id,
                    None,
                )
                .await
                {
                    info!(error = %err, "[recovery.rs] recv task failed");
                }
            });

            loop {
                let Some(event) = rx.recv().await else {
                    return Ok(progress_history);
                };
                match event {
                    Event::Job(job) => {
                        active_jobs = active_jobs.saturating_add(1);
                        self.handle_job_event(
                            job,
                            &lambda_client,
                            &tx,
                            &mut next_part_id,
                        )
                        .await;
                    }
                    Event::Completion(done) => {
                        active_jobs = active_jobs.saturating_sub(1);
                        progress_history.push(done.progress.clone());
                        if done.progress.success() && done.used_temp {
                            any_temp = true;
                            parts.push(done.part_path);
                        }
                    }
                }

                if active_jobs == 0 {
                    if any_temp && !parts.is_empty() {
                        // Merge all successful parts into the downstream output.
                        let aggr = SortMergeAggregator;
                        aggr.aggregate(&parts, &self.fifo_name)?;
                        for p in &parts {
                            let _ = std::fs::remove_file(p);
                        }
                        info!("[recovery.rs][{}] Merged {} parts into final output", self.rdv_tag(), parts.len());
                    }
                    info!("[recovery.rs][{}] All jobs completed", self.rdv_tag());
                    return Ok(progress_history);
                }
            }
        }
        // Send path: single attempt only.
        send(&self.me, &self.peer, &self.rdv_key, &self.fifo_name).await?;
        Ok(progress_history)
    }
}
