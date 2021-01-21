This system handles networked Pash operations using a controller, and
workers. Use for delegating project tasks (like CI) to remote hosts.


## Worker

A worker is a host listening for SSH connections. It must have an
executable `~/work.sh`. The home directory depends on the user account
given to the controller.

The `workers` directory contains useful scripts to run in the context
of a worker. Either link or copy the one you need to `~/work.sh`.


## Controller

The controller is a web server that issues commands to workers.  To
start it, run `node controller/main.js` with desired values for these
environment variables:

* `PASH_REMOTE_HOST`: A hostname or IP address on which an SSH daemon runs. Defaults to `localhost`.
* `PASH_REMOTE_USER`: The username used for authentication. Defaults to `$USER`.
* `PASH_REMOTE_PRIVATE_KEY`: A path to the private key used for authentication.
  Defaults to `/home/${PASH_REMOTE_USER}/.ssh/id_rsa`.

Note that the defaults make it such that the controller and remote
worker are the same machine. Leverage this when prototyping.

Currently, the controller only allows one command to run at a time for
_any_ host. Furthermore, it only protects this invariant using local
state, so restarting the server will make it forget that a worker is
already running. For redundancy, the worker should also try to prevent
any unwanted concurrent executions. This situation should improve with
the addition of shared state.


### API

* `/run` runs `~/work.sh`.

* `/now` dumps the current STDOUT and STDERR for the running command.
   If no command is running, it dumps the STDOUT and STDERR for the
   last-issued command. Note that this shows values accumulated by the
   _controller_ in memory! If you restart the server, it will clear
   the response of this endpoint.

* `/stdout`: like `/now`, but only shows STDOUT.

* `/stderr`: like `/now`, but only shows STDERR.

* `/ci`: like `/run`, with an added invariant. It assumes the worker
   script is `workers/run-ci.sh`, or a script based on it. It
   downloads a `reports` directory from the remote host. *This is for
   backwards compatibility to make sure reports are visible to the
   NGINX server using `scripts/pash-nginx.conf`*.


## Quick Start on Your Machine

Make sure `sshd` is running on your machine (`systemctl start sshd`),
and that you add your public key to `~/.ssh/authorized_keys`. That way
you can log in as yourself using `ssh localhost`.

1. Start the server with `node controller/main.js` in one terminal.
2. Put a script at `~/work.sh`.
3. Run `curl [::]:2047/run`
4. Check `curl [::]:2047/now` to monitor STDOUT and STDERR.


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
