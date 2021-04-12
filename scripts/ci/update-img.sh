#!/bin/bash

###
# Repackage updated pash docker image to latest commit
###

cd $(dirname $0)

# Assumes a pash image exists already
# curl img.pash.ndr.md | docker load; docker run --name pash-playground -it pash/18.04

docker start pash-playground
docker exec  pash-playground bash -c 'cd /pash; git pull'
docker stop  pash-playground


docker commit $(docker ps -a | grep pash-playground | cut -f1 -d' ')  pash/18.04:latest
docker save pash/18.04:latest | gzip > pash-docker.tar.gz

if [[ "$(hostname)" == "beta" ]]; then
  # This assumes you're on beta
  mv pash-docker.tar.gz /var/www/pash-web/
fi

docker build -t pash-play ../

if [[ ./token.txt ]]; then
  cat ~/token.txt | docker login https://docker.pkg.github.com -u nvasilakis --password-stdin
fi

docker push docker.pkg.github.com/andromeda/pash/play:latest