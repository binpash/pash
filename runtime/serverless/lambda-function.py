import subprocess
import json
import os
import boto3
BUCKET=os.environ.get("AWS_BUCKET")

def lambda_handler(event, context):
    for i, folder_id in enumerate(event['folder_ids']):
        id_ = event['ids'][i]
        # load the data from s3
        s3 = boto3.client("s3")
        key = f"sls-scripts/{folder_id}/{id_}.sh"
        print(f"Try to pull script from {key}")

        response = s3.get_object(Bucket=BUCKET, Key=key)
        print("[lambda-function.py] Executing script ID", id_, flush=True)
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
        print(f"[lambda-function.py] script {folder_id}/{id_} execution return code: {process.returncode}")
    print(f"[lambda-function.py] Finished all execution")