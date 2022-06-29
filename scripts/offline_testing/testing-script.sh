rm -f Output.txt
touch Output.txt
# Creating tar file and ready for uploading to docker image building.
tar -cvf pash.tar ../../../pash
versions=(debian:latest ubuntu:latest ubuntu:20.04 ubuntu:18.04) 
cp Dockerfile-Ubuntu Dockerfile
#For every version of Ubuntu, do following things:
for version in ${versions[@]}
do
  # Substituting for every version of Ubuntu inside Dockerfile.
  sed -i "" "s/FROM.*/FROM $version\n/g" Dockerfile
	echo "${version}"
  # When the building fails, output the error message and break.
	if ! sudo docker build -t testing . ;then
    echo "There is an error with current $version."
    break
  fi
  #Image Cleaning
	sudo docker image rm testing
done

#For every version of Fedora, do the same thing as the above.
cp Dockerfile-Fedora Dockerfile
versions=(fedora:latest)
for version in ${versions[@]}
do
  sed -i "" "s/FROM.*/FROM $version\n/g" Dockerfile
	echo "${version}"
	sudo docker build -t testing .
	sudo docker image rm testing
done
rm Dockerfile
