# Installation

On Ubuntu, Fedora, Debian, or Arch, run the following to get PaSh up and running.
```sh
wget https://raw.githubusercontent.com/binpash/pash/main/scripts/up.sh
sh up.sh
export PASH_TOP="$PWD/pash/"
## Run PaSh with echo hi
"$PASH_TOP/pa.sh" -c "echo hi"
```

If on other environments or prefer manual setup, there are essentially three steps required to set PaSh up:
1. Clone repo: `git clone git@github.com:binpash/pash.git`
2. Run `distro-deps.sh` (with `sudo`) and `setup-pash.sh`
3. Export `PASH_TOP` and, optionally, add it to your `PATH`

Most scripts support a `--help` flag that documents their options. These three steps are described in detail below, depending on the environment.


Quick Jump: [Clone & Setup](#) | [Manual Setup](#manual-setup) | [Docker Setup](#docker-setup) | [Windows using WSL](#windows-using-wsl) :


### Clone & Setup

The following steps clone the repo, set up dependencies (e.g., compilers), and then build PaSh:

```sh
git clone git@github.com:binpash/pash.git
./pash/scripts/distro-deps.sh
./pash/scripts/setup-pash.sh
```

These scripts have been tested on Ubuntu, Fedora, Debian, and Arch.


### Docker Setup

If Docker [is available](https://docs.docker.com/get-docker/), PaSh can be used as part of a Docker container.
Depending on the setup, PaSh running on Docker may or may not be able to exploit all available hardware resources.
There are two main options for setting up PaSh on Docker:

_Pull from DockerHub (Major Releases):_
To `pull` the docker image [from Docker Hub](https://hub.docker.com/r/binpash/pash), run:
```sh
docker pull binpash/pash
```
We refresh this image (as well as other images) on every major release.

[//]: # "TODO(@nvasilakis, @dkarnikis): Need to automate this per release."

_Build Image (Latest Commit):_
To build the latest Docker container, run `docker build` in [scripts/docker](https://github.com/binpash/pash/tree/main/scripts/docker):
```sh
git clone git@github.com:binpash/pash.git
cd pash/scripts/docker/
docker build -f ./ubuntu/Dockerfile -t "pash:latest" .
```
In the last line, replace `ubuntu` with `fedora` or `debian` to get the corresponding base.


To run the container:
```sh
docker run --name pash-play -it pash:latest
```

PaSh can be found in the container's `/opt/pash` directory, so run `cd pash; git pull` to fetch the latest updates.
More information in the [pash-on-docker guide](../contributing/contrib.md#pash-on-docker-a-pocket-guide).

### Windows using WSL

To run PaSh on windows without Docker, install [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10).
A short tutorial is included in the [contributing](../contributing/contrib.md) guide.

[//]: # "TODO(@nvasilakis, @dkarnikis): Need to add instructions for OS X."
