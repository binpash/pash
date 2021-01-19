This system handles networked Pash operations, with a focus on CI.

The system consists of a controller, and workers.


## Worker

A worker is an EC2 instance that:

1. has a script located at `~/worker-script.sh` for the AMI's default user.
2. has a tag named `Pash`. The value of the tag reflects the role of that instance in the project.
3. will not run two instances of its script concurrently.

The `workers` directory contains scripts that are each suitable for
use in the first invariant. None of those scripts should depend on
another.


## Controller

The controller is a web server that controls workers.
To start it, run `node controller/main.js`.


## Changelog

### Jan 19 2021

- Decouple service from routes so that implementations can be swapped.

### Jan 18 2021

- Move active code for route implementations from `scripts`.
- Write script to issue SSM command that distribute works to EC2 instances.
- Write script to run tests using Docker and upload the reports to shared storage.
- (External) Provision SSM, EC2, and S3 resources in AWS to prototype implementation.


### Jan 15 2021

- Created `ci` directory
- Wrote `install-aws-cli-v2.sh`, a helper script to install AWS CLI v2
- Wrote draft HMAC comparison to protect the server.
- Refactored code to control more functionality solely in terms of arguments.
