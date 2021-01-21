This system handles networked Pash operations using a controller, and
workers. Use for delegating project tasks (like CI) to remote hosts.


## Worker

A worker is an SSH-enabled host that:

1. has a script located at `~/worker-script.sh` for the AMI's default user.
2. will not run two instances of its script concurrently.

The `workers` directory contains suitable worker scripts.


## Controller

The controller is a web server that controls workers.  To start it,
run `node controller/main.js`.


## Changelog

### Jan 21 2021

- Rename `ci` to `remote`
- Use only vanilla SSH routes in the name of faster turnaround
- Add directory for downloading reports


### Jan 20 2021

- Got `ci.sh` working on EC2 instance
- Add `ssh-routes.js` to issue commands directly to a host.


### Jan 19 2021

- Decouple service from routes so that implementations can be swapped.
- Implement new versions of `/ci` and `/now` to leverage AWS SSM
- Use AWS JS SDK to reduce dependence on Bash subprocesses


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
