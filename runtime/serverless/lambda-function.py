import subprocess
import json

def lambda_handler(event, context):
    data = event["data"]
    id_ = event["id"]
    scripts_dict = json.loads(data)
    with open(f"/tmp/data-{id_}", "w") as f:
        f.write(data)
    print("Executing script ID", id_, flush=True)
    with open(f"/tmp/script-{id_}.sh", "w") as f:
        f.write(scripts_dict[id_])
    # print(f"Script: {scripts_dict[id_]}", flush=True)
    process = subprocess.run(
        ["/bin/bash", f"/tmp/script-{id_}.sh", f"/tmp/data-{id_}"]
    )
    print(f"script {id_} execution return code: {process.returncode}")