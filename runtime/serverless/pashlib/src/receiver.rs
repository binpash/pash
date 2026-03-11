use anyhow::Result;
use tokio::sync::mpsc::UnboundedSender;
use tokio::sync::watch;
use tokio::fs::File;
use tokio::io::{AsyncRead, AsyncReadExt, AsyncWriteExt};
use tracing::info;

use crate::events::{CompletionSpec, Event, JobSpec};
use crate::holepunch::PashCtx;
use crate::metadata::{
    decode_completion_msg, decode_handshake, read_resumability_request, write_resumability_ack,
    resumability_disabled, LambdaMetadata, ResumeAck, StreamMode, COMPLETION_MSG,
    COMPLETION_MSG_SIZE,
    HANDSHAKE_SIZE,
};

const BLOCK_HEADER_SIZE: usize = 24;

// Short rdv key for compact logging.
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
// Track receiver progress for recovery and aggregation.
pub enum ReceiverProgress {
    Chunk {
        success: bool,
        // Chunks completed in the latest attempt.
        last_completed_chunks: u64,
        // Total chunks completed across attempts.
        total_completed_chunks: u64,
        // Blocks forwarded in the current (incomplete) chunk.
        partial_blocks_forwarded: u64,
        // Bytes forwarded in the current block.
        partial_block_bytes_forwarded: u64,
        // Whether the current block header was forwarded.
        partial_header_written: bool,
    },
    Raw {
        success: bool,
        // Bytes forwarded in this attempt.
        last_recv_bytes: u64,
        // Total bytes forwarded across attempts.
        total_recv_bytes: u64,
    },
}

impl ReceiverProgress {
    pub fn success(&self) -> bool {
        match self {
            Self::Chunk { success, .. } | Self::Raw { success, .. } => *success,
        }
    }
}

// Receiver task entrypoint: handles data and control channels and pushes events.
pub async fn recv(
    me: String,
    peer: String,
    rdv_key: String,
    fifo_name: String,
    tx: UnboundedSender<Event>,
    part_id: u64,
    recovery_progress: Option<ReceiverProgress>,
) -> Result<()> {
    info!(
        rdv_key = short_rdv_key(&rdv_key),
        part_id,
        recovery_progress = ?recovery_progress,
        "[receiver.rs][{}] Starting recv task",
        short_rdv_key(&rdv_key)
    );
    // Connect data plane first, then start control monitoring.
    let mut ctx = PashCtx::new(&me, &rdv_key).await;
    let stream = ctx.connect(&peer).await;
    let (mut rd, _wr) = stream.into_split();
    let metadata = match read_metadata(&mut rd, &rdv_key).await {
        Ok(md) => md,
        Err(err) => {
            info!(
                error = %err,
                "[receiver.rs][{}] metadata decode failed; aborting reader",
                short_rdv_key(&rdv_key)
            );
            return Ok(());
        }
    };
    let (shutdown_tx, shutdown_rx) = watch::channel(false);
    let (resume_switch_tx, resume_switch_rx) = watch::channel(false);
    let enable_resumability = !resumability_disabled();
    info!(
        resumability_enabled = enable_resumability,
        "[receiver.rs][{}] Resumability enabled",
        short_rdv_key(&rdv_key)
    );
    let control_handle = if enable_resumability {
        Some(spawn_resumability_listener(
            me.clone(),
            peer.clone(),
            rdv_key.clone(),
            metadata.clone(),
            tx.clone(),
            resume_switch_tx,
            shutdown_rx,
        ))
    } else {
        None
    };

    let recv_mode = if metadata.is_stateless {
        StreamMode::Chunked
    } else {
        StreamMode::RawBytes
    };
    // Compute temp path once; used only if resumability triggers or this is a resume part.
    let temp_path = format!(
        "/tmp/pash_resume_{}_part_{}.part",
        short_rdv_key(&rdv_key),
        part_id
    );
    let start_on_temp = part_id > 1;
    let progress_result = match recv_mode {
        StreamMode::RawBytes => {
            // Raw reader can switch output to a temp file when resumability is triggered.
            match RawReader::new(
                &fifo_name,
                &rdv_key,
                Some(temp_path.clone()),
                if enable_resumability {
                    Some(resume_switch_rx)
                } else {
                    None
                },
                start_on_temp,
            )
            .await
            {
                Ok(mut raw_reader) => {
                    if let Some(ReceiverProgress::Raw { total_recv_bytes, .. }) =
                        recovery_progress
                    {
                        raw_reader.total_received = total_recv_bytes;
                    }
                    let progress = raw_reader.read_from(&mut rd).await;
                    let used_temp = raw_reader.used_temp();
                    (progress, used_temp)
                }
                Err(err) => (Err(err), false),
            }
        }
        StreamMode::Chunked => {
            match ChunkReader::new(&fifo_name, &rdv_key).await {
                Ok(mut chunk_reader) => {
                    if let Some(ReceiverProgress::Chunk {
                        total_completed_chunks,
                        partial_blocks_forwarded,
                        partial_block_bytes_forwarded,
                        partial_header_written,
                        ..
                    }) = recovery_progress
                    {
                        chunk_reader.completed_chunks = total_completed_chunks;
                        chunk_reader.partial_blocks_forwarded = partial_blocks_forwarded;
                        chunk_reader.partial_block_bytes_forwarded = partial_block_bytes_forwarded;
                        chunk_reader.partial_header_written = partial_header_written;
                    }
                    let progress = chunk_reader.read_from(&mut rd).await;
                    (progress, false)
                }
                Err(err) => (Err(err), false),
            }
        }
    };

    // Report completion or enqueue recovery before shutting down control plane.
    match progress_result {
        (Ok(progress), used_temp) => {
            if progress.success() {
                let _ = tx.send(Event::Completion(CompletionSpec {
                    part_id,
                    part_path: temp_path,
                    used_temp,
                    progress: progress.clone(),
                    metadata,
                }));
                info!(
                    part_id,
                    progress = ?progress,
                    "[receiver.rs][{}] Reported successful completion",
                    short_rdv_key(&rdv_key)
                );
            } else {
                // Fault recovery: enqueue job before completion.
                let _ = tx.send(Event::Job(JobSpec {
                    metadata: metadata.clone(),
                    resume_chunk_start_idx: None,
                    recovery_progress: Some(progress.clone()),
                    part_id,
                }));
                info!(
                    part_id,
                    progress = ?progress,
                    "[receiver.rs][{}] Enqueued recovery job due to failure",
                    short_rdv_key(&rdv_key)
                );
                let _ = tx.send(Event::Completion(CompletionSpec {
                    part_id,
                    part_path: temp_path,
                    used_temp: false,
                    progress: progress.clone(),
                    metadata,
                }));
                info!(
                    part_id,
                    progress = ?progress,
                    "[receiver.rs][{}] Reported completion after enqueuing recovery",
                    short_rdv_key(&rdv_key)
                );
            }
        }
        (Err(err), _used_temp) => {
            panic!(
                "[receiver.rs][{}] unexpected recv error: {}",
                short_rdv_key(&rdv_key),
                err
            );
        }
    }
    let _ = shutdown_tx.send(true);
    if let Some(h) = control_handle {
        let _ = h.await;
    }
    Ok(())
}

// Decode the handshake header into LambdaMetadata.
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
    file: Option<File>,
    temp_path: Option<String>,
    temp_file: Option<File>,
    use_temp: bool,
    resume_switch_rx: Option<watch::Receiver<bool>>,
    stdout: tokio::io::Stdout,
    total_received: u64,
    carry: Vec<u8>,
}

impl RawReader {
    // RawReader forwards the byte stream to fifo/stdout and tracks completion markers.
    pub async fn new(
        fifo_name: &str,
        rdv_key: &str,
        temp_path: Option<String>,
        resume_switch_rx: Option<watch::Receiver<bool>>,
        start_on_temp: bool,
    ) -> Result<Self> {
        let to_stdout = fifo_name == "-";
        let file = if to_stdout {
            None
        } else {
            Some(File::create(fifo_name).await?)
        };
        let temp_file = if start_on_temp {
            if let Some(path) = temp_path.as_ref() {
                Some(File::create(path).await?)
            } else {
                None
            }
        } else {
            None
        };
        Ok(Self {
            rdv_key: short_rdv_key(rdv_key).to_string(),
            to_stdout,
            file,
            temp_path,
            temp_file,
            use_temp: start_on_temp,
            resume_switch_rx,
            stdout: tokio::io::stdout(),
            total_received: 0,
            carry: Vec::new(),
        })
    }

    // Switch output to temp file once resumability is triggered.
    async fn maybe_switch_output(&mut self) -> Result<()> {
        if self.use_temp {
            return Ok(());
        }
        let switch = match self.resume_switch_rx.as_ref() {
            Some(rx) => *rx.borrow(),
            None => false,
        };
        if switch {
            if let Some(path) = self.temp_path.as_ref() {
                self.temp_file = Some(File::create(path).await?);
                self.use_temp = true;
            }
        }
        Ok(())
    }

    // Expose whether this reader ever switched to temp output.
    pub fn used_temp(&self) -> bool {
        self.use_temp
    }

    pub async fn read_from<R>(&mut self, reader: &mut R) -> Result<ReceiverProgress>
    where
        R: AsyncRead + Unpin,
    {
        info!(
            "[receiver.rs][{}] RawReader Starting reading data",
            self.rdv_key
        );
        // Stream bytes until the completion marker; track forwarded/skipped bytes for retries.
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
                        last_recv_bytes: recv_bytes,
                        total_recv_bytes: self.total_received,
                    });
                }
                Ok(n) => {
                    self.maybe_switch_output().await?;
                    let (to_stdout, file_ref) = if self.use_temp {
                        (false, &mut self.temp_file)
                    } else {
                        (self.to_stdout, &mut self.file)
                    };
                    if !self.carry.is_empty() {
                        let mut combined = Vec::with_capacity(self.carry.len() + n);
                        combined.extend_from_slice(&self.carry);
                        combined.extend_from_slice(&buf[..n]);

                        if let Some(pos) = find_completion_msg(&combined) {
                            let payload = &combined[..pos];
                            let (discarded, written) = discard_and_write(
                                payload,
                                &mut skip_remaining,
                                to_stdout,
                                file_ref,
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
                                last_recv_bytes: recv_bytes,
                                total_recv_bytes: self.total_received,
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
                                    to_stdout,
                                    file_ref,
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
                                to_stdout,
                                file_ref,
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
                                last_recv_bytes: recv_bytes,
                                total_recv_bytes: self.total_received,
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
                                    to_stdout,
                                    file_ref,
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
                        last_recv_bytes: recv_bytes,
                        total_recv_bytes: self.total_received,
                    });
                }
            }
        }
    }
}

pub struct ChunkReader {
    _rdv_key: String,
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
            _rdv_key: short_rdv_key(rdv_key).to_string(),
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
        info!(
            "[receiver.rs][{}] ChunkReader Starting reading data",
            self._rdv_key
        );
        let mut buf = [0u8; 8192];
        // Retry dedup state:
        // - We do not skip already completed chunks here. Recovery restarts from
        //   the next chunk using chunk_start_idx.
        // - We only skip duplicated blocks/bytes in the current incomplete chunk.
        let mut attempt_completed = 0u64;
        let mut skip_blocks_remaining = self.partial_blocks_forwarded;
        let mut skip_bytes_remaining = self.partial_block_bytes_forwarded;
        let mut header_already_written = self.partial_header_written;

        loop {
            let header = match read_block_header_or_completion(reader).await {
                Ok(Some(header)) => header,
                Ok(None) => {
                    return Ok(ReceiverProgress::Chunk {
                        success: true,
                        last_completed_chunks: attempt_completed,
                        total_completed_chunks: self.completed_chunks,
                        partial_blocks_forwarded: self.partial_blocks_forwarded,
                        partial_block_bytes_forwarded: self.partial_block_bytes_forwarded,
                        partial_header_written: self.partial_header_written,
                    });
                }
                Err(_) => {
                    return Ok(ReceiverProgress::Chunk {
                        success: false,
                        last_completed_chunks: attempt_completed,
                        total_completed_chunks: self.completed_chunks,
                        partial_blocks_forwarded: self.partial_blocks_forwarded,
                        partial_block_bytes_forwarded: self.partial_block_bytes_forwarded,
                        partial_header_written: self.partial_header_written,
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
                            return Ok(ReceiverProgress::Chunk {
                                success: false,
                                last_completed_chunks: attempt_completed,
                                total_completed_chunks: self.completed_chunks,
                                partial_blocks_forwarded: self.partial_blocks_forwarded,
                                partial_block_bytes_forwarded: self.partial_block_bytes_forwarded,
                                partial_header_written: self.partial_header_written,
                            });
                        }
                        Ok(n) => n,
                        Err(_) => {
                            return Ok(ReceiverProgress::Chunk {
                                success: false,
                                last_completed_chunks: attempt_completed,
                                total_completed_chunks: self.completed_chunks,
                                partial_blocks_forwarded: self.partial_blocks_forwarded,
                                partial_block_bytes_forwarded: self.partial_block_bytes_forwarded,
                                partial_header_written: self.partial_header_written,
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
                        return Ok(ReceiverProgress::Chunk {
                            success: false,
                            last_completed_chunks: attempt_completed,
                            total_completed_chunks: self.completed_chunks,
                            partial_blocks_forwarded: self.partial_blocks_forwarded,
                            partial_block_bytes_forwarded: self.partial_block_bytes_forwarded,
                            partial_header_written: self.partial_header_written,
                        });
                    }
                    Ok(n) => n,
                    Err(_) => {
                        self.partial_block_bytes_forwarded =
                            self.partial_block_bytes_forwarded.saturating_add(current_written);
                        self.partial_header_written = header_already_written;
                        return Ok(ReceiverProgress::Chunk {
                            success: false,
                            last_completed_chunks: attempt_completed,
                            total_completed_chunks: self.completed_chunks,
                            partial_blocks_forwarded: self.partial_blocks_forwarded,
                            partial_block_bytes_forwarded: self.partial_block_bytes_forwarded,
                            partial_header_written: self.partial_header_written,
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

fn spawn_resumability_listener(
    me: String,
    peer: String,
    rdv_key: String,
    env_md: LambdaMetadata,
    tx: UnboundedSender<Event>,
    resume_switch_tx: watch::Sender<bool>,
    mut shutdown: watch::Receiver<bool>,
) -> tokio::task::JoinHandle<()> {
    // Control-plane listener: receive resumability requests and enqueue jobs for executor.
    let rdv_tag = short_rdv_key(&rdv_key).to_string();
    let control_key = format!("ctrl::{}", rdv_key);

    tokio::spawn(async move {
        info!(
            "[receiver.rs][{}] Resumability listener started",
            rdv_tag
        );
        if *shutdown.borrow() {
            return;
        }
        let mut ctx = PashCtx::new(&me, &control_key).await;
        let stream = ctx.connect(&peer).await;
        info!(
            control_key = control_key,
            "[receiver.rs][{}] Resumability control channel connected",
            rdv_tag
        );
        let (mut rd, mut wr) = stream.into_split();
        tokio::select! {
            _ = shutdown.changed() => {
                return;
            }
            res = read_resumability_request(&mut rd) => {
                match res {
                    Ok(req) => {
                        let ack = ResumeAck { accepted: true };
                        if let Err(err) = write_resumability_ack(&mut wr, &ack).await {
                            info!(
                                error = %err,
                                "[receiver.rs][{}] control ack write failed",
                                rdv_tag
                            );
                            return;
                        }
                        let _ = resume_switch_tx.send(true);
                        let _ = tx.send(Event::Job(JobSpec {
                            metadata: env_md.clone(),
                            resume_chunk_start_idx: Some(req.next_chunk_start_idx),
                            recovery_progress: None,
                            part_id: 0,
                        }));
                        info!(
                            req = ?req,
                            "[receiver.rs][{}] Resumability request handled and job enqueued",
                            rdv_tag
                        );
                        
                    }
                    Err(err) => {
                        let suppress = err
                            .downcast_ref::<std::io::Error>()
                            .map(|e| e.kind() == std::io::ErrorKind::UnexpectedEof)
                            .unwrap_or(false);
                        if !suppress {
                            info!(
                                error = %err,
                                "[receiver.rs][{}] control channel read failed",
                                rdv_tag
                            );
                        }
                    }
                }
            }
        }
    })
}
