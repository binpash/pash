import subprocess


def lambda_handler(event, context):
    # num = event["num"]
    # data = event["data"]
    # id = event["id"]
    script = f"scripts/curl.sh"

    process = subprocess.Popen(
        ["/bin/bash", script],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    process.wait()

    output, stderr = process.communicate()

    print(f"stdout={output}, stderr={stderr}")

    return output
