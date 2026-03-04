use anyhow::Result;
use std::collections::HashMap;
use std::fmt;

pub const HANDSHAKE_SIZE: usize = 128;
pub const COMPLETION_MSG_SIZE: usize = 16;
pub const COMPLETION_MSG: [u8; COMPLETION_MSG_SIZE] = [
    0x9f, 0x2a, 0x13, 0x8c, 0x44, 0xd1, 0x7b, 0x0e,
    0x6a, 0x55, 0xc9, 0x30, 0x1d, 0x88, 0x02, 0xfe,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamMode {
    RawBytes,
    Chunked,
}

#[derive(Debug, Clone, Default)]
pub struct LambdaMetadata {
    pub leash_job_id: String,
    pub folder_id: String,
    pub script_id: String,
    pub chunk_start_idx: u32,
    pub is_stateless: bool,
}

impl LambdaMetadata {
    pub fn from_env() -> Self {
        let leash_job_id = std::env::var("LEASH_JOB_ID").unwrap_or_default();
        let folder_id = std::env::var("PASH_FOLDER_ID").unwrap_or_default();
        let script_id = std::env::var("PASH_SCRIPT_ID").unwrap_or_default();
        let chunk_start_idx = std::env::var("PASH_CHUNK_START_IDX")
            .ok()
            .and_then(|v| v.parse::<u32>().ok())
            .unwrap_or(0);
        let is_stateless = std::env::var("PASH_IS_STATELESS")
            .ok()
            .map(|v| matches!(v.as_str(), "1" | "true" | "True" | "TRUE"))
            .unwrap_or(false);

        Self {
            leash_job_id,
            folder_id,
            script_id,
            chunk_start_idx,
            is_stateless,
        }
    }

    pub fn from_metadata_blob(blob: &str, chunk_start_idx: u32, is_stateless: bool) -> Self {
        let kv = parse_metadata_kv(blob);
        Self {
            leash_job_id: kv.get("LEASH_JOB_ID").cloned().unwrap_or_default(),
            folder_id: kv.get("FOLDERS_ID").cloned().unwrap_or_default(),
            script_id: kv.get("SCRIPT_ID").cloned().unwrap_or_default(),
            chunk_start_idx,
            is_stateless,
        }
    }

    pub fn to_metadata_blob(&self) -> String {
        format!(
            "LEASH_JOB_ID={};FOLDERS_ID={};SCRIPT_ID={}",
            self.leash_job_id, self.folder_id, self.script_id
        )
    }
}

impl fmt::Display for LambdaMetadata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LambdaMetadata {{ leash_job_id: {}, folder_id: {}, script_id: {}, chunk_start_idx: {}, is_stateless: {} }}",
            self.leash_job_id, self.folder_id, self.script_id, self.chunk_start_idx, self.is_stateless
        )
    }
}

fn parse_metadata_kv(metadata: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for pair in metadata.split(';') {
        let Some((k, v)) = pair.split_once('=') else {
            continue;
        };
        out.insert(k.to_string(), v.to_string());
    }
    out
}

pub fn encode_handshake(metadata: &LambdaMetadata) -> [u8; HANDSHAKE_SIZE] {
    let mut out = [0u8; HANDSHAKE_SIZE];
    out[0] = u8::from(metadata.is_stateless);
    out[1..5].copy_from_slice(&metadata.chunk_start_idx.to_be_bytes());

    let payload_blob = metadata.to_metadata_blob();
    let payload = payload_blob.as_bytes();
    let payload_len = payload.len().min(HANDSHAKE_SIZE - 5);
    out[5..5 + payload_len].copy_from_slice(&payload[..payload_len]);
    out
}

pub fn decode_handshake(buf: &[u8; HANDSHAKE_SIZE]) -> Result<LambdaMetadata> {
    let is_stateless = buf[0] != 0;
    let chunk_start_idx = u32::from_be_bytes([buf[1], buf[2], buf[3], buf[4]]);
    let metadata_blob = String::from_utf8_lossy(&buf[5..])
        .trim_end_matches('\0')
        .to_string();
    Ok(LambdaMetadata::from_metadata_blob(
        &metadata_blob,
        chunk_start_idx,
        is_stateless,
    ))
}

pub fn decode_completion_msg(buf: &[u8; COMPLETION_MSG_SIZE]) -> bool {
    buf == &COMPLETION_MSG
}
