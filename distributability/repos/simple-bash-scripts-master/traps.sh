#!/bin/bash

function finish() {
	# your cleanup here.
	echo "CTL+C pressed"
	echo "clean ..."
	sleep 1
}
trap finish EXIT

echo 'runing ...'
until false; do
	sleep 1
done
