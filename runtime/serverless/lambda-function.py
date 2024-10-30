# import subprocess
# import json

# def lambda_handler(event, context):
#     data = event["data"]
#     id_ = event["id"]
#     scripts_dict = json.loads(data)
#     with open(f"/tmp/data-{id_}", "w") as f:
#         f.write(data)
#     print("[lambda-function.py] Executing script ID", id_, flush=True)
#     with open(f"/tmp/script-{id_}.sh", "w") as f:
#         f.write(scripts_dict[id_])
#     # print(f"Script: {scripts_dict[id_]}", flush=True)
#     process = subprocess.run(
#         ["/bin/bash", f"/tmp/script-{id_}.sh", f"/tmp/data-{id_}"]
#     )
#     print(f"[lambda-function.py] script {id_} execution return code: {process.returncode}")

import subprocess
import json
import os
import boto3
BUCKET=os.environ.get("AWS_BUCKET")

def lambda_handler(event, context):
    folder_id = event["folder_id"]
    id_ = event["id"]
    # load the data from s3
    s3 = boto3.client("s3")
    key = f"sls-scripts/{folder_id}/{id_}.sh"

    response = s3.get_object(Bucket=BUCKET, Key=key)
    print("[lambda-function.py] Executing script ID", id_, flush=True)
    with open(f"/tmp/script-{id_}.sh", "wb") as f:
        while True:
            x = response["Body"].read(10000)
            if not x:
                break
            f.write(x)
            f.flush()
    # with open(f"/tmp/script-{id_}.sh", "r") as f:
    #     print(f"Script: {f.read()}", flush=True)
    process = subprocess.run(
        ["/bin/bash", f"/tmp/script-{id_}.sh", folder_id]
    )
    print(f"[lambda-function.py] script {id_} execution return code: {process.returncode}")