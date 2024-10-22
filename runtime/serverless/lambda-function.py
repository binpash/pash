import subprocess
import json
import boto3
import os

BUCKET=os.environ.get("AWS_BUCKET")

def read_data_s3(s3_filepath: str):
    print("[lambda-function.py, read_data_s3] Reading from file:", s3_filepath, "in", BUCKET)
    session = boto3.Session()
    s3 = session.client("s3", region_name='us-east-1')
    try:
        response = s3.get_object(Bucket=BUCKET, Key=s3_filepath)
        batch = 10000
        result = ""
        while True:
            x = response["Body"].read(batch)
            if not x:
                break
            result += x.decode("utf-8")
        return result
    except Exception as e:
        print(e)
        return None

def lambda_handler(event, context):
    data = read_data_s3(event["data"])
    if data is None: # an error occurred reading the script_id_to_script mapping
        return

    id_ = event["id"]
    scripts_dict = json.loads(data)
    if id_ not in scripts_dict:
        print(f"[lambda-function.py] ERROR: {id_} NOT IN SCRIPTS_DICT")
        print(scripts_dict)
    else:
        print("[lambda-function.py] Executing script ID", id_, flush=True)
        with open(f"/tmp/script-{id_}.sh", "w") as f:
            f.write(scripts_dict[id_])
        # print(f"Script: {scripts_dict[id_]}", flush=True)
        process = subprocess.run(["/bin/bash", f"/tmp/script-{id_}.sh", event["data"]])
        print(f"[lambda-function.py] script {id_} execution return code: {process.returncode}")
