#!/usr/bin/env python3
"""
Manual S3-Based Orchestrator for Serverless Sort
Bypasses holepunch/pashlib complexity by using S3 for all data transfer
"""

import boto3
import subprocess
import tempfile
import os
import time
import json
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

class S3Orchestrator:
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

    def download_input(self, local_path):
        """Download input file from S3"""
        print(f"[Step 1] Downloading {self.input_key} from S3...")
        start = time.time()
        self.s3.download_file(self.bucket, self.input_key, local_path)
        elapsed = time.time() - start
        size_mb = os.path.getsize(local_path) / (1024 * 1024)
        print(f"  ✓ Downloaded {size_mb:.2f} MB in {elapsed:.2f}s")

    #TODO. use split script instead
    def split_file(self, input_path, num_chunks):
        """Split file round-robin into N chunks, return chunk paths"""
        print(f"[Step 2] Splitting file into {num_chunks} chunks...")
        start = time.time()

        # Create chunk files
        chunk_files = []
        chunk_handles = []
        for i in range(num_chunks):
            f = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=f'_chunk{i}.txt')
            chunk_files.append(f.name)
            chunk_handles.append(f)

        # Round-robin distribution
        line_count = 0
        with open(input_path, 'r') as infile:
            for line_num, line in enumerate(infile):
                chunk_idx = line_num % num_chunks
                chunk_handles[chunk_idx].write(line)
                line_count += 1

        # Close all handles
        for f in chunk_handles:
            f.close()

        elapsed = time.time() - start
        print(f"  ✓ Split {line_count} lines in {elapsed:.2f}s")
        for i, path in enumerate(chunk_files):
            size_kb = os.path.getsize(path) / 1024
            print(f"    Chunk {i}: {size_kb:.1f} KB")

        return chunk_files

    def upload_chunks(self, chunk_paths):
        """Upload chunks to S3, return S3 keys"""
        print(f"[Step 3] Uploading {len(chunk_paths)} chunks to S3...")
        start = time.time()

        chunk_keys = []
        for i, chunk_path in enumerate(chunk_paths):
            chunk_key = f"chunks/{self.folder_id}/chunk{i}.txt"
            self.s3.upload_file(chunk_path, self.bucket, chunk_key)
            chunk_keys.append(chunk_key)
            print(f"  ✓ Uploaded chunk {i} to s3://{self.bucket}/{chunk_key}")

        elapsed = time.time() - start
        print(f"  ✓ Uploaded all chunks in {elapsed:.2f}s")
        return chunk_keys

    def invoke_lambda_worker(self, chunk_key, worker_id):
        """Invoke Lambda to sort a chunk (synchronous)"""
        print(f"[Worker {worker_id}] Invoking Lambda...")
        start = time.time()

        payload = {
            'chunk_key': chunk_key,
            'folder_id': self.folder_id,
            'worker_id': worker_id,
            'bucket': self.bucket
        }

        response = self.lambda_client.invoke(
            FunctionName='lambda-sort',
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
        print(f"[Step 5] Downloading sorted results from S3...")
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


    #TODO: use merge script instead
    def merge_sorted_chunks(self, chunk_paths, output_path):
        """Merge sorted chunks using sort -m"""
        print(f"[Step 6] Merging {len(chunk_paths)} sorted chunks...")
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
        print(f"[Step 7] Uploading final result to S3...")
        start = time.time()

        self.s3.upload_file(local_path, self.bucket, self.output_key)

        elapsed = time.time() - start
        print(f"  ✓ Uploaded to s3://{self.bucket}/{self.output_key} in {elapsed:.2f}s")

    def execute(self):
        """Main orchestration logic"""
        print("=" * 60)
        print("Starting S3-Based Sort Orchestration")
        print("=" * 60)
        total_start = time.time()

        # Step 1: Download input
        input_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='_input.txt').name
        self.download_input(input_file)

        # Step 2: Split into chunks
        chunk_paths = self.split_file(input_file, self.num_workers)

        # Step 3: Upload chunks to S3
        chunk_keys = self.upload_chunks(chunk_paths)

        # Step 4: Invoke Lambda workers in parallel
        print(f"[Step 4] Invoking {self.num_workers} Lambda workers in parallel...")
        worker_start = time.time()

        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            futures = []
            for i, chunk_key in enumerate(chunk_keys):
                future = executor.submit(self.invoke_lambda_worker, chunk_key, i)
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

        worker_elapsed = time.time() - worker_start
        print(f"  ✓ All workers completed in {worker_elapsed:.2f}s")

        # Step 5: Download sorted results
        sorted_chunks = self.download_sorted_chunks()

        # Step 6: Merge sorted chunks
        output_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='_output.txt').name
        self.merge_sorted_chunks(sorted_chunks, output_file)

        # Step 7: Upload final result
        self.upload_result(output_file)

        # Cleanup temporary files
        print(f"[Cleanup] Removing temporary files...")
        for path in [input_file] + chunk_paths + sorted_chunks + [output_file]:
            if os.path.exists(path):
                os.remove(path)
        print(f"  ✓ Cleaned up temporary files")

        total_elapsed = time.time() - total_start
        print("=" * 60)
        print(f"✓ Orchestration completed successfully in {total_elapsed:.2f}s")
        print(f"  Output: s3://{self.bucket}/{self.output_key}")
        print("=" * 60)

def main():
    parser = argparse.ArgumentParser(description='Manual S3-based sort orchestrator')
    parser.add_argument('--bucket', required=True, help='S3 bucket name')
    parser.add_argument('--input', required=True, help='S3 key for input file')
    parser.add_argument('--output', required=True, help='S3 key for output file')
    parser.add_argument('--workers', type=int, default=2, help='Number of Lambda workers')

    args = parser.parse_args()

    orchestrator = S3Orchestrator(
        bucket=args.bucket,
        input_key=args.input,
        output_key=args.output,
        num_workers=args.workers
    )

    orchestrator.execute()

if __name__ == '__main__':
    main()



"""
   python3 manual_s3_orchestrator.py \
       --bucket inout741448956691 \
       --input oneliners/inputs/1M.txt \
       --output oneliners/outputs/my-sorted-output.txt \
       --workers 2
"""    

"""

TIMING SUMMARY - FULL VERSION
============================================================
  1. Download input............................   0.29s
  2. Split into chunks.........................   0.02s
  3. Upload chunks.............................   0.20s
  4. Lambda workers (parallel).................   1.00s
  5. Download sorted results...................   0.15s
  6. Merge sorted chunks.......................   0.00s
  7. Upload final result.......................   0.12s
------------------------------------------------------------
  TOTAL TIME...................................   1.78s
  """


"""
TIMING SUMMARY - FULL VERSION
============================================================
  1. Download input............................  12.85s
  2. Split into chunks.........................  15.73s
  3. Upload chunks.............................   4.48s
  4. Lambda workers (parallel).................  38.18s
  5. Download sorted results...................  11.12s
  6. Merge sorted chunks.......................   8.44s
  7. Upload final result.......................   4.48s
------------------------------------------------------------
  TOTAL TIME...................................  95.84s
=========================================================
    """
    #TODO. compare with benchmark time to see if makes sense 
    # sequential so not clear if would be diff with streaming 
    # TODO. convert to shell script to make it easier to stream 
    # add timing for the downloading of s3-> lambda 
    # because of streaming have to measure end to end time 
    # when rm 1-3 then get seq time for exp 2



    #in a sequential setting. in a plot. send on discord. 

    # 1. how much time: pull 1gb/ 500 MB from s3 to lambda
    # 2. how much time: pull 1 gb/ 500 mB from s3 to ec2 
    # 3. how much time: byte range of 500 MB from 1GB from s3 to lambda 

    # q: does byte range mess up perf?
    # q: does the network input speed for lambda getting from ec2 differ from lambda getting from s3?


    #Hypotheses:
    # 1. s3 to ec2 ~= s3 to lambda 
    # 2. s3 to lambda with file size ~= byte range pulling of the same size 

    #discuss 

    #then run streaming hardcoded with baseline (no dup s3) and incorrect splitting 
