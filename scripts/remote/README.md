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
start it, run `node controller/main.js` with a desired [configuration](#runtime-configuration).

Note that the defaults make it such that the controller and remote
worker are the same machine. Leverage this when prototyping.


### Limitations

A controller process only allows one command to run at a time for its
worker. Furthermore, it only protects this invariant using local
state, so restarting the server will make it forget A) if the worker
is already running, and B) anything it knows about a command running
on the worker.

You should not (yet?) rely on the controller to protect a worker's
invariants. For now, that's a worker's problem. This situation will
improve with the addition of shared state. Scaling to use more hosts
is a matter of spawning more controller processes with different
runtime configurations.


### Runtime Configuration

The `PASH_REMOTE_RCFILE` environment variable specifies the location
of the runtime configuration file ("rcfile"). Make sure this is an
absolute path, otherwise the process will try to load configuration
w.r.t. the working directory. By default, `PASH_REMOTE_RCFILE` is set
to the absolute path of `controller/rc.json` on the host system.
It does not have to exist.

The rcfile is a JSON document. Some useable examples are in
`controller/configs`. A controller uses the following keys:

* `host`: A hostname or IP address on which an SSH daemon runs. Defaults to `localhost`.

* `user`: The username used for authentication. Defaults to `process.env.USER`.

* `private_key`: A path to the private key used for authentication. Defaults to `${process.env.HOME}/.ssh/id_rsa`.

* `ci_reports_path`: The destination path used by `/ci` to download reports from its worker.
  Defaults to the (completed) path of `controller/reports` on the host system.

* `command_module`: A path to a (trusted!) CommonJS module used to
  configure the controller's relationship with the worker. You need to
  set this if you want to define when the controller considers the
  command "done", along with any custom behavior to run after that. If
  this is a relative path, then the path is resolves w.r.t. the
  complete path to the `controller` directory. You will likely set
  this to a file in `controller/commands`,
  e.g. `./commands/correctness.js`. See [Command
  Modules][cm] for built-in modules.

If the key is not defined in the JSON document--or the rcfile does not
exist--,then the process will check the current environment variables
instead. The envvar name is the upcased version of the rcfile key,
with a `PASH_REMOTE_` prefix, and underscores replacing all
non-alphanumeric characters.

e.g. `host` -> `PASH_REMOTE_HOST`, or `my-conf1g-key` -> `PASH_REMOTE_MY_CONF1G_KEY`


### API

* `/run` runs the worker's `~/work.sh`, with added invariants
  depending on the server's [command module][cm].

* `/now` dumps the current STDOUT and STDERR for the running command.
   If no command is running, it dumps the STDOUT and STDERR for the
   last-issued command. The actual response body depends on
   [Limitations](#limitations).

* `/stdout`: like `/now`, but only shows STDOUT.

* `/stderr`: like `/now`, but only shows STDERR.


## Quick Start

Make sure `sshd` is running on your machine (`systemctl start sshd`),
and that you add your public key to `~/.ssh/authorized_keys`. That way
you can log in as yourself using `ssh localhost`.

Follow these instructions _without_ making an SSH connection.

1. Start the server with `node controller/main.js` in one terminal.
2. Put a script at `~/work.sh`.
3. Run `curl [::]:2047/run`
4. Check `curl [::]:2047/now` to monitor/review STDOUT and STDERR.
5. Repeat steps 3-4 to iterate until you get desired output.

Once you are satisfied with how the script works, update the
`PASH_REMOTE_*` environment variables to the (presumably more
powerful) machine you want to run the script from that point on. Copy
`work.sh` over, then repeat step 5 to iterate on the behavior of the
worker.


## Command Modules
[cm]: #command-modules

`commands/default.js` simply runs the worker script without waiting
for it to complete.

`commands/correctness.js` assumes the worker script is
`workers/run-ci.sh`, or a script based on it. If the controller
detects relevant output, it downloads a `reports` directory from the
remote host. *This is for backwards compatibility to make sure reports
are visible to the NGINX server using `scripts/pash-nginx.conf`*.


## Changelog

### Jan 26 2021

- Small refactor to allow configurable `/run` endpoint.
- Remove `/ci`, now that it was made redundant due to the above change.
- Add `command_module` and `port` settings


### Jan 25 2021

- Add runtime configuration module
- Use `rsync` to fix issues with downloading reports directory
- Add setting to control where reports directory gets downloaded


### Jan 21 2021

- Rename files to clarify their purpose.
- Use only vanilla SSH routes in the name of faster turnaround
- Add documentation
- Add submodule to protect command invariants.


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
