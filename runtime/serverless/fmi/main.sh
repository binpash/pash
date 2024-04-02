#!/bin/bash

rm -f sender receiver

aws lambda invoke \
	--function-name lambda \
	--invocation-type RequestResponse \
	--cli-binary-format raw-in-base64-out \
	--payload '{"id": "'0'"}' \
	sender &

sleep 1

aws lambda invoke \
	--function-name lambda \
	--invocation-type RequestResponse \
	--cli-binary-format raw-in-base64-out \
	--payload '{"id": "'1'"}' \
	receiver &
