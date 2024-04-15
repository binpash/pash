import subprocess
import json
import time

def lambda_handler(event, context):
    # data = event["data"]
    # id_ = event["id"]
    # scripts_dict = json.loads(data)
    # with open(f"/tmp/data-{id_}", "w") as f:
        # f.write(data)
    # print("[lambda-function.py] Executing script ID", id_, flush=True)
    # with open(f"/tmp/script-{id_}.sh", "w") as f:
        # f.write(scripts_dict[id_])
    # print(f"Script: {scripts_dict[id_]}", flush=True)
    start_time = time.time()
    process = subprocess.run(
        ["/bin/bash", f"demo-spell.sh"]
    )
    end_time = time.time()
    total_time = end_time - start_time
    total_time = round(total_time, 3)
    print(f"[lambda-function.py] script execution return code: {process.returncode}")
    print(f"[lambda-function.py] script execution time: {total_time}")
