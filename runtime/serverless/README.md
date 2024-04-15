# Serverless runtime components

To deploy the serverless runtime as a function on AWS Lambda:
```sh
# First, export the following environment variables.
# Ideally, you should put these to ~/.bashrc.
export AWS_ACCOUNT_ID="Your AWS account id here"
export AWS_QUEUE="Your queue id here"
export AWS_BUCKET="Your bucket id here"

# Then, deploy to AWS:
sls deploy
```
