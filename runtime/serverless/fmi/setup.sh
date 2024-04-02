#!/bin/bash

rm -f lambda.zip

cd lambda
zip -r ../lambda.zip .
cd ..


for functionName in lambda; do
	aws lambda delete-function --function-name $functionName
	aws lambda create-function \
		--function-name $functionName \
		--zip-file fileb://lambda.zip \
		--handler lambda-function.lambda_handler \
		--runtime python3.7 \
		--role arn:aws:iam::347768412644:role/lambda-ex \
		--memory 2048 \
		--layers arn:aws:lambda:us-east-1:347768412644:layer:fmi-python37:1 \
		--timeout 20 | cat
done
