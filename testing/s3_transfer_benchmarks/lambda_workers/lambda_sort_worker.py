"""
Lambda Sort Worker - S3-Based
Downloads a chunk from S3, sorts it, and uploads the result back to S3
"""

import boto3
import subprocess
import os
import json

def lambda_handler(event, context):
    """
    Event format:
    {
        "chunk_key": "chunks/12345/chunk0.txt",
        "folder_id": "12345",
        "worker_id": 0,
        "bucket": "my-bucket"
    }
    """
    print(f"[Lambda Worker] Starting...")
    print(f"Event: {json.dumps(event)}")

    # Extract parameters
    chunk_key = event['chunk_key']
    folder_id = event['folder_id']
    worker_id = event['worker_id']
    bucket = event['bucket']

    s3 = boto3.client('s3')

    # File paths
    input_file = f"/tmp/chunk_{worker_id}.txt"
    output_file = f"/tmp/sorted_{worker_id}.txt"
    result_key = f"results/{folder_id}/sorted{worker_id}.txt"

    try:
        # Step 1: Download chunk from S3
        print(f"[Worker {worker_id}] Downloading {chunk_key}...")
        s3.download_file(bucket, chunk_key, input_file)
        input_size = os.path.getsize(input_file)
        print(f"[Worker {worker_id}] Downloaded {input_size} bytes")

        # Step 2: Sort the chunk
        print(f"[Worker {worker_id}] Sorting...")
        result = subprocess.run(
            ['sort', input_file, '-o', output_file],
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout
        )

        if result.returncode != 0:
            error_msg = f"Sort failed: {result.stderr}"
            print(f"[Worker {worker_id}] ERROR: {error_msg}")
            return {
                'statusCode': 500,
                'body': json.dumps({'error': error_msg})
            }

        output_size = os.path.getsize(output_file)
        print(f"[Worker {worker_id}] Sorted {output_size} bytes")

        # Step 3: Upload sorted result to S3
        print(f"[Worker {worker_id}] Uploading to {result_key}...")
        s3.upload_file(output_file, bucket, result_key)
        print(f"[Worker {worker_id}] Upload complete")

        # Cleanup
        if os.path.exists(input_file):
            os.remove(input_file)
        if os.path.exists(output_file):
            os.remove(output_file)

        print(f"[Worker {worker_id}] SUCCESS")
        return {
            'statusCode': 200,
            'result_key': result_key,
            'input_size': input_size,
            'output_size': output_size
        }

    except subprocess.TimeoutExpired:
        error_msg = "Sort operation timed out"
        print(f"[Worker {worker_id}] ERROR: {error_msg}")
        return {
            'statusCode': 500,
            'body': json.dumps({'error': error_msg})
        }

    except Exception as e:
        error_msg = str(e)
        print(f"[Worker {worker_id}] ERROR: {error_msg}")
        return {
            'statusCode': 500,
            'body': json.dumps({'error': error_msg})
        }
