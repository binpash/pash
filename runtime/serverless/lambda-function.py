import subprocess
import json
import os
import threading
import time
import boto3
BUCKET=os.environ.get("AWS_BUCKET")

def lambda_handler(event, context):
    job_id = event.get('job_id', 'UNKNOWN')
    os.environ['LEASH_JOB_ID'] = job_id

    time_to_crash = event.get('time_to_crash')
    if time_to_crash:
        crash_start = time.time()
        def _crash_timer():
            time.sleep(time_to_crash)
            elapsed = time.time() - crash_start
            print(f"[JOB:{job_id}] [CRASH] Crashing after {time_to_crash}s, actual elapsed: {elapsed:.3f}s", flush=True)
            os._exit(1)
        threading.Thread(target=_crash_timer, daemon=True).start()

    chunk_start_idx = event.get('chunk_start_idx')
    if chunk_start_idx is not None:
        os.environ['PASH_CHUNK_START_IDX'] = str(chunk_start_idx)
    else:
        # Ensure the env var is always set cuz lambda can reuse instances across different jobs.
        os.environ['PASH_CHUNK_START_IDX'] = '0'

    is_stateless = event.get('is_stateless')
    if is_stateless is not None:
        os.environ['PASH_IS_STATELESS'] = str(is_stateless).lower()
    else:
        # Ensure the env var is always set cuz lambda can reuse instances across different jobs.
        os.environ['PASH_IS_STATELESS'] = 'false'
    
    timeout = event.get('timeout')
    if timeout:
        os.environ['PASH_RESUME_TIMEOUT_SEC'] = str(timeout)
    else:
        os.environ['PASH_RESUME_TIMEOUT_SEC'] = '0'

    print(f"[JOB:{job_id}] Args: {json.dumps(event, sort_keys=True)}", flush=True)
    print(
        f"[JOB:{job_id}] Env: LEASH_JOB_ID={os.environ.get('LEASH_JOB_ID', 'unset')} "
        f"PASH_CHUNK_START_IDX={os.environ.get('PASH_CHUNK_START_IDX', 'unset')} "
        f"PASH_IS_STATELESS={os.environ.get('PASH_IS_STATELESS', 'unset')} "
        f"PASH_RESUME_TIMEOUT_SEC={os.environ.get('PASH_RESUME_TIMEOUT_SEC', 'unset')}",
        flush=True,
    )

    for i, folder_id in enumerate(event['folder_ids']):
        id_ = event['ids'][i]
        os.environ['PASH_FOLDER_ID'] = folder_id
        os.environ['PASH_SCRIPT_ID'] = id_
        # load the data from s3
        s3 = boto3.client("s3")
        key = f"sls-scripts/{folder_id}/{id_}.sh"
        print(f"[JOB:{job_id}] Try to pull script from {key}")

        response = s3.get_object(Bucket=BUCKET, Key=key)
        print(f"[JOB:{job_id}] [lambda-function.py] Executing script ID {id_}", flush=True)
        with open(f"/tmp/script-{folder_id}-{id_}.sh", "wb") as f:
            while True:
                x = response["Body"].read(10000)
                if not x:
                    break
                f.write(x)
                f.flush()
        # with open(f"/tmp/script-{id_}.sh", "r") as f:
        #     print(f"Script: {f.read()}", flush=True)
        process = subprocess.run(
            ["/bin/bash", f"/tmp/script-{folder_id}-{id_}.sh", folder_id]
        )
        print(f"[JOB:{job_id}] [lambda-function.py] script {folder_id}/{id_} execution return code: {process.returncode}")
    print(f"[JOB:{job_id}] [lambda-function.py] Finished all execution")
