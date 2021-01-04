This directory contains CI and installation scripts for PaSh:

* `ci.sh`: Majority of CI job, including packaging and testing. N.b.: Not re-entrant.
* `install.sh`: Install dependencies for running PaSh locally. Assumes clean install, otherwise can mess with system-wide config.
* `now.sh`: Prints information about the currently executing job.
* `pash-nginx.conf`: Web server config, including routes. Assumes a certain dir structure in `/www` and that `webhook` is running.
* `pkg.sh`: Package PaSh for use by `up.sh` and other downloads.
* `up.sh`: Download and install PaSh.
* `webhook.js`: Endpoints for `ctrl.pash.ndr.md`.
* `clone_compress_repo.sh`: Clone and package repo. Assumes access to GitHub.
* `ssh-install.sh`: Push to a remote directory
