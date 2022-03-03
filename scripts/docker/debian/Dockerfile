FROM debian:10
RUN apt-get update -y && apt-get install -y git
ENV PASH_HOST="docker"
ENV PASH_TOP=/opt/pash
# download PaSh
RUN git clone https://github.com/binpash/pash.git /opt/pash --depth 1
RUN bash /opt/pash/scripts/distro-deps.sh -o 
RUN yes | bash /opt/pash/scripts/setup-pash.sh -o
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV LC_COLLATE C
ENV LC_ALL en_US.UTF-8
CMD ["/bin/bash"]
