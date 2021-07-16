
# PaSh Docker Support

Docker files for [PaSh](https://github.com/binpash/pash) published on Dockerhub:

* Ubuntu 18.04: [binpash/pash:ubuntu-18.04](https://hub.docker.com/r/binpash/pash)
* Fedora 35 : [binpash/pash:fedora-35](https://hub.docker.com/r/binpash/pash)
* Debian 10: [binpash/pash:debian-10](https://hub.docker.com/r/binpash/pash)

## Run locally

The following will pull and launch the latest image (tag: `pash:latest`) using the name `pash-play`

```sh
docker pull binpash/pash
docker run --name pash-play -it pash
```

To restart after you exit, run `docker start -i pash-play`

## From PaSh's Official Dockerfiles

```sh
git clone git@github.com:binpash/pash.git
cd pash/scripts/docker
docker build -t "pash:latest" .
```

## As a Dockerfile

```
FROM ubuntu:18.04

# ARG CACHEBUST=1

RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y git sudo locales locales-all curl wget wamerican-insane

RUN git clone https://github.com/andromeda/pash.git

RUN sudo bash /pash/scripts/install.sh -p

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

ENV PASH_TOP=/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

CMD ["/bin/bash"]
```

and build and run this locally:

```sh
$ docker build -t app . && docker run -it --init -p 1993:1993 app
```

## (optional) Add `pash-docker` as an alias to your shell

You can add `pash` command to your shell init file (e.g.`.bashrc`) to launch scripts in the Docker environment.

```sh
pash () {
  docker start -i pash-play
    --interactive \
    --tty \
    --rm \
    --volume $PWD:/pwd \
    --workdir /pwd \
    pash-play \
    "$@"
}
```

and in your terminal:

```sh
source ~/.bashrc
pash --version
pash -w=2 hello.sh
```

_This will not work if the scipts depend on environment variables on the environment hosting the Docker container._

### Preparing Docker Images

* build the image: `docker build -t "binpash/pash:$distro-$version"`
* push the image `docker push binpash/pash:$distro-$version`

[//]: # (TODO(@nvasilakis, @dkarnikis): Just add a script.)
