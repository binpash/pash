import sys
import time
import boto3
import os
sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
import config
from util import log

def wait_msg_done_dynamo(key):
    while True:
        dynamo = boto3.client('dynamodb')
        table_name = 'sls-result'
        try:
            response = dynamo.get_item(
                TableName=table_name,
                Key={
                    'key': {
                        'S': key
                    }
                }
            )
            item = response['Item']
            print('[Serverless Client] Received message: %s' % item, file=sys.stderr)
            break
        except:
            # log("Waiting for message")
            time.sleep(1)
    return

if __name__ == "__main__":
    key = sys.argv[1:][0]
    wait_msg_done_dynamo(key)