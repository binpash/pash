FROM debian:10

# ARG CACHEBUST=1

RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y git sudo locales locales-all curl wget wamerican-insane

RUN git clone https://github.com/andromeda/pash.git

RUN sudo bash /pash/scripts/install.sh -p

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

ENV PASH_TOP=/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

CMD ["/bin/bash"]
