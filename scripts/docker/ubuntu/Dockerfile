FROM ubuntu:18.04
RUN apt-get update -y && apt-get install -y git
RUN apt-get update -y && apt-get install -y sudo
ENV PASH_TOP=/opt/pash
# download PaSh
RUN git clone --branch caruca-eval https://github.com/binpash/pash.git /opt/pash
RUN bash /opt/pash/scripts/distro-deps.sh -o 
RUN pip3 install Cython
RUN yes | bash /opt/pash/scripts/setup-pash.sh -o
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV LC_COLLATE C
ENV LC_ALL en_US.UTF-8
CMD ["/bin/bash"]
