FROM fedora:latest
RUN dnf install git -y
RUN git clone https://github.com/qishenli/pash.git && \
   pash/scripts/distro-deps.sh && \
#   pash/scripts/setup-pash.sh && \
   pash/pa.sh pash/evaluation/intro/hello-world.sh
