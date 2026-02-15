import boto3
import sys
import json
import os
import time

AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
BUCKET=os.environ.get("AWS_BUCKET")
object_key, infile = sys.argv[1:][:2]
DEBUG=True

if DEBUG:
    print(f"[s3-put-object.py] Start putting {object_key}")
    print(f"[s3-put-object.py] Input file: {infile}")
    print(f"[s3-put-object.py] S3 bucket: {BUCKET}")
    start_time = time.time()

try:
    # Wrapper class to count bytes as they're read
    class ByteCountingReader:
        def __init__(self, fileobj):
            self.fileobj = fileobj
            self.bytes_read = 0

        def read(self, size=-1):
            data = self.fileobj.read(size)
            self.bytes_read += len(data)
            return data

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return self.fileobj.__exit__(*args)
        
    if object_key == "/dev/null":
        if DEBUG:
            print(f"[s3-put-object.py] Get /dev/null as key, stream data to /dev/null.")
        with open(infile, "rb") as f:
            while f.read(1024 * 1024):
                pass  # Just read through the file without doing anything
            sys.exit(0)

    s3 = boto3.client('s3')

    with open(infile, "rb") as f:
        reader = ByteCountingReader(f)
        s3.upload_fileobj(reader, BUCKET, object_key)
        file_size = reader.bytes_read

    if DEBUG:
        end_time = time.time()
        total_time = end_time - start_time
        total_time = round(total_time, 3)
        print(f"[s3-put-object.py] Successfully uploaded {object_key}")
        print(f"[s3-put-object.py] Upload size: {file_size} bytes in {total_time} seconds")
        if file_size > 0 and total_time > 0:
            throughput = round(file_size / total_time / 1024 / 1024, 2)
            print(f"[s3-put-object.py] Throughput: {throughput} MB/s")
        elif file_size == 0:
            print(f"[s3-put-object.py] WARNING: Uploaded 0 bytes!")

except FileNotFoundError as e:
    print(f"[s3-put-object.py] ERROR: Input file not found: {infile}", file=sys.stderr)
    print(f"[s3-put-object.py] Exception: {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"[s3-put-object.py] ERROR during upload: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
