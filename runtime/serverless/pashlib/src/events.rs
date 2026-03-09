use crate::metadata::LambdaMetadata;
use crate::receiver::ReceiverProgress;

// Job specification passed through the executor event loop (always invokes lambda).
#[derive(Clone)]
pub struct JobSpec {
    // Lambda metadata is required for both recovery and resume invocation.
    pub metadata: LambdaMetadata,
    // Presence means resumability (start from next chunk).
    pub resume_chunk_start_idx: Option<u32>,
    // Reader recovery progress used for dedup on retry.
    pub recovery_progress: Option<ReceiverProgress>,
    // 0 means "assign in executor".
    pub part_id: u64,
}

// Completion record emitted by a reader.
pub struct CompletionSpec {
    // Part id is used to reuse the same output file on recovery retries.
    pub part_id: u64,
    pub part_path: String,
    // Whether the reader switched to temp output for this part.
    pub used_temp: bool,
    pub progress: ReceiverProgress,
    pub metadata: LambdaMetadata,
}

// Executor event stream.
pub enum Event {
    // Job comes from control (resume) or recovery retries
    Job(JobSpec),
    // All readers will send a completion event to control upon completion
    Completion(CompletionSpec),
}
