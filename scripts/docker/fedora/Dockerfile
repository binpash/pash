FROM fedora:35
RUN dnf install git -y
ENV PASH_HOST="docker"
# download PaSh
RUN git clone https://github.com/binpash/pash.git /opt/pash --depth 1
# install PaSh and setup shell paths
RUN yes | bash /opt/pash/scripts/install.sh -p
CMD ["/bin/bash"]
