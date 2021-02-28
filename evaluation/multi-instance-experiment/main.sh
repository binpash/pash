#!/usr/bin/env bash

## This script is just an experiment trying to:
##  1. spawn EC2 instances
##  2. run some tests on them
##  3. collect the results
##  4. stop the instances
##  5. plot the results 

## running a script on an instance is done using with-ec2.sh
##
## the ssh wrapper to run on an ubuntu ec2 instance is in ssh-script.sh
##
## the actual script that is run on the ec2 instance is experiment-script.sh
##
## the instance profiles to use when creating an instance are saves in `instance-profiles`

## Q: We need to make sure to have configurable steps that are already done (e.g. reinstalling, running tests)

## Q: What if we need to create a new key pair

## Q: Is it safe to save the instance ids in public?

## TODO: Refactor this to run for an arbitrary experiment so that it can be reused.

key_path="~/.ssh/aws_pash.pem"
standard_disk_instance_id="i-0347068fae17c256e"
standard_disk_pretty_name="standard_disk_instance"
fast_disk_instance_id="i-07814ed42a1ddd013"
fast_disk_pretty_name="fast_disk_instance"

eval_dir="$PASH_TOP/evaluation/multi-instance-experiment/"
local_res_dir="$eval_dir/results"

## Assumes that script is installed and instance exists.
## This can be done by sourcing a profile:
## ```sh
## source instance-profiles/fast-disk.sh
## ```
## and then running `scripts/aws/ec2/make-ec2.sh`

execute_on_instance_and_collect_results()
{
    local instance_id=$1
    local instance_res_dir="$local_res_dir/$2"

    ## Execute and Collect results
    $PASH_TOP/scripts/with-ec2.sh "$instance_id" "$eval_dir/execute-and-collect.sh" "$instance_id" "$key_path" "$instance_res_dir"
}

execute_on_instance_and_collect_results "$standard_disk_instance_id" "$standard_disk_pretty_name"
execute_on_instance_and_collect_results "$fast_disk_instance_id" "$fast_disk_pretty_name"

## TODO: Add plotting
