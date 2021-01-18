This is a draft for an AWS-based CI system.


## Use

To start the controller service, run `node main.js`.

To add a worker for the controller to use, see [here](#add-a-worker).


## High-Level Design

The system consists of a controller, and workers.

A worker is an EC2 instance that runs an assigned script on demand.

The controller is a web server that tells at least one worker to run
their scripts, and manages the status information reported by all
of the workers. It can be hosted wherever.

To be clear, the controller does not care _what_ the workers do, only
that they identify themselves and write up reports. This frees up
programmers to create a worker for performance tests, a worker
for correctness tests, etc.

Workers upload information to S3, and the controller helps users
navigate and interpret the results.

This is functionally not much different than the prior version,
but the implementation leverages more cloud computing features.


## Add a Worker

To add a worker to the system:

1. Create the EC2 instance.
2. Connect to the instance and write its instructions in an executable `~/worker-script.sh`
3. In AWS, create a `Pash` tag on the instance and set the value to `worker`.


## Changelog

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
