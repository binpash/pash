#!/usr/bin/env python3
"""
EC2 script to benchmark S3 download performance
Measures time to pull 1G or 500M file from S3
"""
import boto3
import botocore.config
import time
import uuid

BUCKET = "inout741448956691"

KEYS = [
    "unix50/inputs/1_20G.txt",
    "oneliners/inputs/1G.txt",
    "oneliners/inputs/500M.txt",
    "oneliners/inputs/100M.txt",
]

KEY = KEYS[0]  # Change this to select different file sizes

def run_benchmark(event=None, context=None):
    runtimes = []

    for i in range(10):
        # Create a fresh S3 client each iteration to avoid connection reuse
        config = botocore.config.Config(
            max_pool_connections=1,
            retries={'max_attempts': 0}
        )
        s3 = boto3.client("s3", config=config)

        # Add a random query param to avoid any caching
        unique_key = f"{KEY}?nocache={uuid.uuid4().hex}"

        t0 = time.time()
        res = s3.get_object(Bucket=BUCKET, Key=KEY)
        for _ in res['Body'].iter_chunks(chunk_size=1024*1024):
            pass

        t1 = time.time()

        #if i == 0: # first run is cold start garbage
        #    continue

        dt = t1 - t0
        runtimes.append(dt)
        print(f"Run {i+1}: {dt:.2f}s")#, size pulled: {len(data)/1e6:.1f} MB")

        #del data  # free memory

    # Compute statistics
    avg_ = sum(runtimes) / len(runtimes)
    mn = min(runtimes)
    mx = max(runtimes)
    sorted_r = sorted(runtimes)
    p90_index = max(0, int(len(sorted_r) * 0.9) - 1)
    p90 = sorted_r[p90_index]

    print(f"Average: {avg_:.2f}s")
    print(f"Min: {mn:.2f}s")
    print(f"Max: {mx:.2f}s")
    print(f"P90: {p90:.2f}s")

    return {
        "avg": avg_,
        "min": mn,
        "max": mx,
        "p90": p90,
        "runs": runtimes
    }

if __name__ == "__main__":
    # Run benchmark on EC2
    result = run_benchmark()
    print("\nFinal results:", result)
