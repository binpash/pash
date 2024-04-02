#!/bin/bash

# aws ec2 create-key-pair --key-name key --query 'KeyMaterial' --output text > key.pem
# chmod 400 key.pem

vpcId=$(aws ec2 create-vpc --cidr-block "10.0.0.0/16" --output text --query 'Vpc.VpcId')
subnetId=$(aws ec2 create-subnet --vpc-id $vpcId --cidr-block "10.0.0.0/24" --availability-zone us-east-1a --output text --query 'Subnet.SubnetId')
ec2GroupId=$(aws ec2 create-security-group --group-name ec2 --description "For ec2 instances" --vpc-id $vpcId --output text --query 'GroupId')

internetGatewayId=$(aws ec2 create-internet-gateway --output text --query 'InternetGateway.InternetGatewayId')
routeTableId=$(aws ec2 describe-route-tables --filters "Name=vpc-id,Values=$vpcId" --output text --query 'RouteTables[*].Associations[*].RouteTableId')

aws ec2 attach-internet-gateway --vpc-id $vpcId --internet-gateway-id $internetGatewayId
aws ec2 create-route --route-table-id $routeTableId --destination-cidr-block "0.0.0.0/0" --gateway-id $internetGatewayId

instanceId=$(aws ec2 run-instances \
	--image-id ami-01bc990364452ab3e \
	--count 1 \
	--instance-type t2.micro \
	--key-name key \
	--tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=test}]' \
	--security-group-ids $ec2GroupId \
	--subnet-id $subnetId \
	--associate-public-ip-address \
	--output text \
	--query 'Instances[*].InstanceId')

aws ec2 wait instance-running --instance-ids $instanceId

aws ec2 authorize-security-group-ingress --group-id $ec2GroupId --protocol tcp --port 22 --cidr "0.0.0.0/0"
aws ec2 authorize-security-group-ingress --group-id $ec2GroupId --protocol tcp --port 80 --cidr "10.0.0.0/24"
aws ec2 authorize-security-group-ingress --group-id $ec2GroupId --protocol tcp --port 10000 --cidr "10.0.0.0/24"
aws ec2 authorize-security-group-ingress --group-id $ec2GroupId --protocol tcp --port "0-65535" --cidr "10.0.0.0/24"

instanceIp=$(aws ec2 describe-instances --instance-ids $instanceId --output text --query 'Reservations[*].Instances[*].PublicIpAddress')

ssh -i key.pem ec2-user@$instanceIp

lambdaGroupId=$ec2GroupId
aws lambda update-function-configuration --function-name lambda --vpc-config SubnetIds=$subnetId,SecurityGroupIds=$lambdaGroupId
