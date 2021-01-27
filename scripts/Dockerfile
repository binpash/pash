FROM ubuntu:18.04

RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install git sudo -y

RUN git clone https://github.com/andromeda/pash.git
WORKDIR "/pash"

RUN sed -i 's#git@github.com:angelhof/libdash.git#https://github.com/angelhof/libdash/#g' .gitmodules
RUN sudo bash scripts/install.sh -p
RUN export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

ENV PASH_TOP=/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
CMD ["/bin/bash"]
