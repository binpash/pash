# Serverless runtime components

To deploy the serverless runtime as a function on AWS Lambda:
```sh
# First, export the following environment variables.
# Ideally, you should put these to ~/.bashrc.
export AWS_ACCOUNT_ID="Your AWS account id here"
export AWS_QUEUE="Your queue id here"
export AWS_BUCKET="Your bucket id here"

# Second, prepare all necessary runtime binary in lambda
cd $PASH_TOP/runtime/serverless
cp ../{r_merge,r_split,r_unwrap,r_wrap,split} runtime/

# Then, deploy to AWS:
sls deploy
```
