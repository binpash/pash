import json
import os
import re
import socket
import sys
import time
import boto3
import logging

logger = logging.getLogger(__name__)
AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
QUEUE=os.environ.get("AWS_QUEUE")
EC2_IP=os.environ.get("AWS_EC2_IP")

def exec():
    pass


def invoke_lambda_ec2(script_id_to_script, script_id):
    EC2_PORT = 9999
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((EC2_IP, EC2_PORT))
    json_data = json.dumps({"id": script_id, "data": json.dumps(script_id_to_script)})
    # send length first?
    # s.sendall(str(len(encoded_json_data))+"\r\n".encode("utf-8"))
    try:
        s.sendall(json_data.encode("utf-8"))
        s.close()
    except Exception as e:
        logger.error(f"EC2 {script_id} invocation error: {e}")
    logger.info(f"EC2 {script_id} invocation")

def wait_msg_done():
    while True:
        sqs = boto3.client('sqs')
        queue_url = f'https://sqs.us-east-1.amazonaws.com/{AWS_ACCOUNT_ID}/{QUEUE}'
        # Receive message from SQS queue
        response = sqs.receive_message(
            QueueUrl=queue_url,
            AttributeNames=[
                'SentTimestamp'
            ],
            MaxNumberOfMessages=1,
            MessageAttributeNames=[
                'All'
            ],
            VisibilityTimeout=30,
            WaitTimeSeconds=20
        )
        try:
            message = response['Messages'][0]
            receipt_handle = message['ReceiptHandle']

            # Delete received message from queue
            sqs.delete_message(
                QueueUrl=queue_url,
                ReceiptHandle=receipt_handle
            )
            logger.info('Received and deleted message: %s' % message)
            break
        except:
            time.sleep(1)

def invoke_lambda(script_id_to_script, script_id):
    lambda_client = boto3.client("lambda",region_name='us-east-1')
    response = lambda_client.invoke(
        FunctionName="lambda",
        InvocationType="Event",
        LogType="None",
        Payload=json.dumps({"data": json.dumps(script_id_to_script) , "id": script_id}),
    )
    return response

# Function to replace all $variable occurrences in the script with their environment variable values
def replace_variables_in_script(script):

    # Find all $x and ${x} variables in the script
    variables = re.findall(r'\$([A-Za-z_][A-Za-z0-9_]*)|\${([A-Za-z_][A-Za-z0-9_]*)}', script)

    # Replace each variable with its corresponding environment variable value
    for var1, var2 in variables:
        var = var1 if var1 else var2  # Choose the one that isn't empty
        env_value = os.getenv(var, f"<{var}>")  # Default to <var> if the env variable is not found
        script = script.replace(f"${{{var}}}", env_value)  # Replace ${x} syntax
        script = script.replace(f"${var}", env_value)  # Replace $x syntax
    
    return script


def main(bash_path, instance):
    with open(bash_path, 'r') as f:
        bash_script = f.read()
    export_prefix = "export PATH=$PATH:runtime \n export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib\n"
    script_id_to_script = {
        "bash": export_prefix+replace_variables_in_script(bash_script)
    }
    main_subgraph_script_id = "bash"
    logger.info(script_id_to_script)
    # start serverless execution by invoking the first lambda
    if instance == 'lambda':
        response = invoke_lambda(script_id_to_script, main_subgraph_script_id)
        logger.info(response)
        wait_msg_done()
    elif instance == 'hybrid':
        invoke_lambda_ec2(script_id_to_script, main_subgraph_script_id)
        wait_msg_done()

if __name__ == '__main__':
    bash_path = sys.argv[1]
    instance = sys.argv[2]
    main(bash_path, instance)


#  python aws/s3-get-object.py \"oneliners/inputs/1G.txt\" \"/dev/stdout\" | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | python aws/s3-put-object.py oneliners/outputs/wf.bash/stdout.txt \"/dev/stdin\"

# python aws/s3-get-object.py \"oneliners/inputs/1G.txt\" \"/dev/stdout\" | col -bx | tr -cs A-Za-z '\n' | tr A-Z a-z | tr -d '[:punct:]' | sort | uniq | comm -23 - <(python aws/s3-get-object.py \"oneliners/inputs/dict.txt\" \"/dev/stdout\") | python aws/s3-put-object.py oneliners/outputs/spell.bash/stdout.txt \"/dev/stdin\"

# export PATH=$PATH:runtime\n export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib\n rm -rf \"/tmp/fifo1\" \n mkfifo \"/tmp/fifo1\" \n python aws/s3-get-object.py \"oneliners/inputs/dict.txt\" \"/tmp/fifo1\" & \n python aws/s3-get-object.py \"oneliners/inputs/1G.txt\" \"/dev/stdout\" | col -bx | tr -cs A-Za-z '\n' | tr A-Z a-z | tr -d '[:punct:]' | sort | uniq | comm -23 - \"/tmp/fifo1\" | python aws/s3-put-object.py oneliners/outputs/spell.bash/stdout.txt \"/dev/stdin\"