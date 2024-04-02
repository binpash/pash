aws iam create-policy --policy-name EC2Policy --policy-document file://policies/ec2-policy.json
aws iam attach-role-policy --role-name lambda-ex --policy-arn
