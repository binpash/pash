import subprocess


def lambda_handler(event, context):
    id = int(event["id"])

    if id == 0:
        script = "sender.sh"
    elif id == 1:
        script = "receiver.sh"

    process = subprocess.Popen(
        ["/bin/bash", script],
        # stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE,
    )

    process.wait()

    # output, stderr = process.communicate()

    # return f"stdout={output}, stderr={stderr}"
    return "Done"
