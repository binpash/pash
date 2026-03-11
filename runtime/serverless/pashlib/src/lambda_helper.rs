use anyhow::{anyhow, Result};
use aws_sdk_lambda::primitives::Blob;
use aws_sdk_lambda::types::InvocationType;
use aws_sdk_lambda::Client;
use serde::Serialize;
use tracing::info;

use crate::metadata::LambdaMetadata;

fn short_rdv_key(rdv_key: &str) -> &str {
    rdv_key.get(..6).unwrap_or(rdv_key)
}

#[derive(Debug, Clone, Serialize)]
pub struct LambdaInvokePayload {
    pub job_id: String,
    pub folder_ids: Vec<String>,
    pub ids: Vec<String>,
    pub chunk_start_idx: u32,
    pub is_stateless: bool,
}

async fn _invoke_lambda(
    lambda_client: &Client,
    function_name: &str,
    invocation_type: InvocationType,
    payload: Vec<u8>,
    rdv_key: &str,
) -> Result<()> {
    let _rdv_key = short_rdv_key(rdv_key);
    let resp = lambda_client
        .invoke()
        .function_name(function_name)
        .invocation_type(invocation_type)
        .payload(Blob::new(payload))
        .send()
        .await?;

    let status = resp.status_code;
    if !(200..=299).contains(&status) {
        return Err(anyhow!(
            "lambda invoke failed with status {} for function {}",
            status,
            function_name
        ));
    }
    Ok(())
}

pub async fn invoke_lambda(
    lambda_client: &Client,
    function_name: &str,
    metadata: &LambdaMetadata,
    chunk_start_idx: u32,
    rdv_key: &str,
) -> Result<()> {
    let payload = LambdaInvokePayload {
        job_id: metadata.leash_job_id.clone(),
        folder_ids: vec![metadata.folder_id.clone()],
        ids: vec![metadata.script_id.clone()],
        chunk_start_idx,
        is_stateless: metadata.is_stateless,
    };
    let rdv_key = short_rdv_key(rdv_key);
    info!(
        function_name = %function_name,
        job_id = %payload.job_id,
        folder_ids = ?payload.folder_ids,
        ids = ?payload.ids,
        chunk_start_idx = payload.chunk_start_idx,
        is_stateless = payload.is_stateless,
        "[lambda_helper.rs][{}] lambda args",
        rdv_key
    );
    let payload_bytes = serde_json::to_vec(&payload)?;
    _invoke_lambda(
        lambda_client,
        function_name,
        InvocationType::Event,
        payload_bytes,
        rdv_key,
    )
    .await
}
