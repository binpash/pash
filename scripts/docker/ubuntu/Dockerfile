FROM ubuntu:18.04
RUN apt-get update -y
RUN apt-get upgrade -y
# fetch all the require deps
RUN apt-get install -y git sudo locales locales-all curl wget wamerican-insane software-properties-common python3-pip
RUN python3 -m pip install -U --force-reinstall pip
# download PaSh
RUN git clone https://github.com/binpash/pash.git /opt/pash
# install PaSh
RUN sudo bash /opt/pash/scripts/install.sh -p
# install to /usr/bin directory
RUN sudo cp /opt/pash/pa.sh /usr/bin/
# add to the bashrc
RUN echo "export PASH_TOP=/opt/pash" >> ~/.bashrc
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
# export PaSh directory
ENV PASH_TOP=/opt/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
CMD ["/bin/bash"]
