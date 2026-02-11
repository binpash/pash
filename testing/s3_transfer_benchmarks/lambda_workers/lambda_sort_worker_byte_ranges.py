"""
Lambda Sort Worker - S3 Byte Range Based
Downloads a byte range from S3, sorts it, and uploads the result back to S3
"""

import boto3
import subprocess
import os
import json

def lambda_handler(event, context):
    """
    Event format:
    {
        "input_key": "oneliners/inputs/1M.txt",
        "start_byte": 0,
        "end_byte": 499999,
        "is_last_chunk": false,
        "folder_id": "12345",
        "worker_id": 0,
        "bucket": "my-bucket"
    }
    """
    print(f"[Lambda Worker] Starting...")
    print(f"Event: {json.dumps(event)}")

    # Extract parameters
    input_key = event['input_key']
    start_byte = event['start_byte']
    end_byte = event['end_byte']
    is_last_chunk = event.get('is_last_chunk', False)
    folder_id = event['folder_id']
    worker_id = event['worker_id']
    bucket = event['bucket']

    s3 = boto3.client('s3')

    # File paths
    input_file = f"/tmp/chunk_{worker_id}.txt"
    output_file = f"/tmp/sorted_{worker_id}.txt"
    result_key = f"results/{folder_id}/sorted{worker_id}.txt"

    try:
        # Step 1: Download byte range from S3
        byte_range = f"bytes={start_byte}-{end_byte}"
        print(f"[Worker {worker_id}] Downloading {input_key} range {byte_range}...")

        response = s3.get_object(
            Bucket=bucket,
            Key=input_key,
            Range=byte_range
        )

        # Read the data
        data = response['Body'].read()
        original_size = len(data)
        print(f"[Worker {worker_id}] Downloaded {original_size} bytes")

        # Step 2: Handle line boundaries
        # If not the last chunk, truncate at the last complete newline
        if not is_last_chunk:
            last_newline_pos = data.rfind(b'\n')
            if last_newline_pos != -1:
                data = data[:last_newline_pos + 1]  # Keep the newline
                truncated_size = len(data)
                print(f"[Worker {worker_id}] Truncated at last newline: {original_size} -> {truncated_size} bytes")
            else:
                print(f"[Worker {worker_id}] Warning: No newline found in chunk, keeping all data")
        else:
            print(f"[Worker {worker_id}] Last chunk - keeping all data including final line")

        # Write to input file
        with open(input_file, 'wb') as f:
            f.write(data)

        # Step 3: Sort the chunk
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

        # Step 4: Upload sorted result to S3
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
            'input_size': original_size,
            'processed_size': len(data),
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
