
# PaSh Docker Support

Docker files for [PaSh](https://github.com/binpash/pash) published on Dockerhub:

* Ubuntu 18.04: [binpash/pash:ubuntu-18.04](https://hub.docker.com/r/binpash/pash)
* Fedora 35: [binpash/pash:fedora-35](https://hub.docker.com/r/binpash/pash)
* Debian 10: [binpash/pash:debian-10](https://hub.docker.com/r/binpash/pash)

## Run locally

The following will pull and launch the latest image (tag: `pash:latest`) using the name `pash-play`

```sh
docker pull binpash/pash
docker run --name pash-play -it pash
```

To restart after you exit, run `docker start -i pash-play`

## From PaSh's Official docker

```sh
git clone --depth 1 git@github.com:binpash/pash.git
cd pash/scripts/docker
docker build -t "pash:latest" .
```

## As a Dockerfile

* Ubuntu 18.04: [rendered](https://github.com/binpash/pash/blob/main/scripts/docker/ubuntu/Dockerfile); [raw](https://raw.githubusercontent.com/binpash/pash/main/scripts/docker/ubuntu/Dockerfile)
* Fedora 35: [rendered](https://github.com/binpash/pash/blob/main/scripts/docker/ubuntu/Dockerfile); [raw](https://raw.githubusercontent.com/binpash/pash/main/scripts/docker/ubuntu/Dockerfile)
* Debian 10: [rendered](https://github.com/binpash/pash/blob/main/scripts/docker/ubuntu/Dockerfile); [raw](https://raw.githubusercontent.com/binpash/pash/main/scripts/docker/ubuntu/Dockerfile)

To build any of these containers and run it locally:

```sh
# URL points to ubuntu, pick any of the "raw" URLs above
wget https://raw.githubusercontent.com/binpash/pash/main/scripts/docker/ubuntu/Dockerfile
docker build -t "pash:stable" . && docker run -it "pash:stable"
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

_This will not work if the scripts depend on environment variables on the environment hosting the Docker container._

### Preparing Docker Images

* build the image: `docker build -t "binpash/pash:$distro-$version"`
* push the image `docker push binpash/pash:$distro-$version`

[//]: # "TODO(@nvasilakis, @dkarnikis): Just add a script."
