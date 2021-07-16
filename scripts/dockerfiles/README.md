# build the image
docker build -t "binpash/pash:$distro-$version" .
# push the image
docker push binpash/pash:$distro-$version
