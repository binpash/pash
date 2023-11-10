import subprocess
import json
import os

def lambda_handler(event, context):
    data = event["data"]
    id_ = event["id"]
    scripts_dict = json.loads(data)
    with open(f"/tmp/{id_}", "w") as f:
        f.write(data)
    print("Executing script ID", id_)
    with open("/tmp/script.sh", "w") as f:
        f.write(scripts_dict[id_])
    script = f"/tmp/script.sh"
    # with open(script, "r") as f:
    #     for line in f.readlines():
    #         print(line)
    # print(os.getcwd())
    process = subprocess.Popen(
        ["/bin/bash", script, f"/tmp/{id_}"], cwd=os.getcwd()
    )

    process.wait()
    output, _ = process.communicate()
    print(output)
    return output