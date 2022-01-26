FROM debian:10
RUN apt-get update -y
RUN apt-get install -y git 
ENV PASH_HOST="docker"
# download PaSh
RUN git clone https://github.com/binpash/pash.git /opt/pash --depth 1
# install PaSh
RUN bash /opt/pash/scripts/install.sh -p
CMD ["/bin/bash"]
