#!/usr/bin/env python3
"""
Manual S3-Based Orchestrator for Serverless Sort (Byte Range Version)
Lambdas read byte ranges directly from S3 input file - no local split needed
"""

import boto3
import subprocess
import tempfile
import os
import time
import json
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

class S3OrchestratorByteRanges:
    def __init__(self, bucket, input_key, output_key, num_workers=2):
        self.s3 = boto3.client('s3')
        self.lambda_client = boto3.client('lambda', region_name='us-east-1')
        self.bucket = bucket
        self.input_key = input_key
        self.output_key = output_key
        self.num_workers = num_workers
        self.folder_id = str(int(time.time()))

        print(f"[Orchestrator] Initialized")
        print(f"  Bucket: {bucket}")
        print(f"  Input: {input_key}")
        print(f"  Output: {output_key}")
        print(f"  Workers: {num_workers}")
        print(f"  Folder ID: {self.folder_id}")

    def get_file_size(self):
        """Get the size of the input file from S3"""
        print(f"[Step 1] Getting file size from S3...")
        start = time.time()

        response = self.s3.head_object(Bucket=self.bucket, Key=self.input_key)
        file_size = response['ContentLength']

        elapsed = time.time() - start
        size_mb = file_size / (1024 * 1024)
        print(f"  ✓ File size: {size_mb:.2f} MB ({file_size} bytes) in {elapsed:.2f}s")
        return file_size

    def calculate_byte_ranges(self, file_size):
        """Calculate byte ranges for each worker"""
        print(f"[Step 2] Calculating byte ranges for {self.num_workers} workers...")

        chunk_size = file_size // self.num_workers
        byte_ranges = []

        for i in range(self.num_workers):
            start_byte = i * chunk_size
            # Last worker gets everything remaining (handles rounding)
            end_byte = file_size - 1 if i == self.num_workers - 1 else (i + 1) * chunk_size - 1
            is_last_chunk = (i == self.num_workers - 1)

            byte_ranges.append({
                'worker_id': i,
                'start_byte': start_byte,
                'end_byte': end_byte,
                'is_last_chunk': is_last_chunk
            })

            range_size_kb = (end_byte - start_byte + 1) / 1024
            print(f"  Worker {i}: bytes {start_byte}-{end_byte} ({range_size_kb:.1f} KB)" +
                  (", last chunk" if is_last_chunk else ""))

        return byte_ranges

    def invoke_lambda_worker(self, byte_range):
        """Invoke Lambda to sort a byte range (synchronous)"""
        worker_id = byte_range['worker_id']
        print(f"[Worker {worker_id}] Invoking Lambda...")
        start = time.time()

        payload = {
            'input_key': self.input_key,
            'start_byte': byte_range['start_byte'],
            'end_byte': byte_range['end_byte'],
            'is_last_chunk': byte_range['is_last_chunk'],
            'folder_id': self.folder_id,
            'worker_id': worker_id,
            'bucket': self.bucket
        }

        response = self.lambda_client.invoke(
            FunctionName='lambda-sort-ranges',
            InvocationType='RequestResponse',  # Synchronous
            Payload=json.dumps(payload)
        )

        elapsed = time.time() - start

        if response['StatusCode'] == 200:
            result = json.loads(response['Payload'].read())
            print(f"[Worker {worker_id}] ✓ Completed in {elapsed:.2f}s")
            return result.get('result_key')
        else:
            error = response.get('FunctionError', 'Unknown error')
            print(f"[Worker {worker_id}] ✗ Failed: {error}")
            raise Exception(f"Lambda invocation failed for worker {worker_id}: {error}")

    def download_sorted_chunks(self):
        """Download sorted results from S3"""
        print(f"[Step 4] Downloading sorted results from S3...")
        start = time.time()

        sorted_chunks = []
        for i in range(self.num_workers):
            result_key = f"results/{self.folder_id}/sorted{i}.txt"
            local_path = tempfile.NamedTemporaryFile(mode='w+b', delete=False, suffix=f'_sorted{i}.txt').name

            try:
                self.s3.download_file(self.bucket, result_key, local_path)
                sorted_chunks.append(local_path)
                size_kb = os.path.getsize(local_path) / 1024
                print(f"  ✓ Downloaded sorted chunk {i}: {size_kb:.1f} KB")
            except Exception as e:
                print(f"  ✗ Failed to download {result_key}: {e}")
                raise

        elapsed = time.time() - start
        print(f"  ✓ Downloaded all sorted chunks in {elapsed:.2f}s")
        return sorted_chunks

    def merge_sorted_chunks(self, chunk_paths, output_path):
        """Merge sorted chunks using sort -m"""
        print(f"[Step 5] Merging {len(chunk_paths)} sorted chunks...")
        start = time.time()

        # Use sort -m for efficient merge
        cmd = ['sort', '-m'] + chunk_paths + ['-o', output_path]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  ✗ Merge failed: {result.stderr}")
            raise Exception(f"Merge failed: {result.stderr}")

        elapsed = time.time() - start
        size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"  ✓ Merged to {size_mb:.2f} MB in {elapsed:.2f}s")

    def upload_result(self, local_path):
        """Upload final result to S3"""
        print(f"[Step 6] Uploading final result to S3...")
        start = time.time()

        self.s3.upload_file(local_path, self.bucket, self.output_key)

        elapsed = time.time() - start
        print(f"  ✓ Uploaded to s3://{self.bucket}/{self.output_key} in {elapsed:.2f}s")

    def execute(self):
        """Main orchestration logic"""
        print("=" * 60)
        print("Starting S3-Based Sort Orchestration (BYTE-RANGE VERSION)")
        print("=" * 60)
        total_start = time.time()
        timings = {}

        # Step 1: Get file size
        step_start = time.time()
        file_size = self.get_file_size()
        timings['1. Get file size'] = time.time() - step_start

        # Step 2: Calculate byte ranges
        step_start = time.time()
        byte_ranges = self.calculate_byte_ranges(file_size)
        timings['2. Calculate byte ranges'] = time.time() - step_start

        # Step 3: Invoke Lambda workers in parallel
        print(f"[Step 3] Invoking {self.num_workers} Lambda workers in parallel...")
        step_start = time.time()

        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            futures = []
            for byte_range in byte_ranges:
                future = executor.submit(self.invoke_lambda_worker, byte_range)
                futures.append(future)

            # Wait for all workers to complete
            results = []
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    print(f"Worker failed: {e}")
                    raise

        worker_elapsed = time.time() - step_start
        timings['3. Lambda workers (parallel)'] = worker_elapsed
        print(f"  ✓ All workers completed in {worker_elapsed:.2f}s")

        # Step 4: Download sorted results
        step_start = time.time()
        sorted_chunks = self.download_sorted_chunks()
        timings['4. Download sorted results'] = time.time() - step_start

        # Step 5: Merge sorted chunks
        step_start = time.time()
        output_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='_output.txt').name
        self.merge_sorted_chunks(sorted_chunks, output_file)
        timings['5. Merge sorted chunks'] = time.time() - step_start

        # Step 6: Upload final result
        step_start = time.time()
        self.upload_result(output_file)
        timings['6. Upload final result'] = time.time() - step_start

        # Cleanup temporary files
        print(f"[Cleanup] Removing temporary files...")
        for path in sorted_chunks + [output_file]:
            if os.path.exists(path):
                os.remove(path)
        print(f"  ✓ Cleaned up temporary files")

        total_elapsed = time.time() - total_start

        # Print timing summary
        print("\n" + "=" * 60)
        print("TIMING SUMMARY - BYTE-RANGE VERSION")
        print("=" * 60)
        for step, duration in timings.items():
            print(f"  {step:.<45} {duration:>6.2f}s")
        print("-" * 60)
        print(f"  {'TOTAL TIME':.<45} {total_elapsed:>6.2f}s")
        print("=" * 60)
        print(f"Output: s3://{self.bucket}/{self.output_key}")
        print("=" * 60)

def main():
    parser = argparse.ArgumentParser(
        description='Manual S3-based sort orchestrator using byte ranges'
    )
    parser.add_argument('--bucket', required=True, help='S3 bucket name')
    parser.add_argument('--input', required=True, help='S3 key for input file')
    parser.add_argument('--output', required=True, help='S3 key for output file')
    parser.add_argument('--workers', type=int, default=2, help='Number of Lambda workers')

    args = parser.parse_args()

    orchestrator = S3OrchestratorByteRanges(
        bucket=args.bucket,
        input_key=args.input,
        output_key=args.output,
        num_workers=args.workers
    )

    orchestrator.execute()

if __name__ == '__main__':
    main()
