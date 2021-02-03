FROM ubuntu:18.04


RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y git sudo locales locales-all curl

# According to an email, Node v8.x is close to lowest stakeholder version.
# Adjust PPA accordingly.
RUN curl -sL https://deb.nodesource.com/setup_8.x -o nodesource_setup.sh
RUN bash nodesource_setup.sh
RUN apt-get update -y
RUN apt-get install -y nodejs npm

RUN git clone https://github.com/andromeda/pash.git
WORKDIR "/pash"

RUN sed -i 's#git@github.com:angelhof/libdash.git#https://github.com/angelhof/libdash/#g' .gitmodules
RUN sudo bash scripts/install.sh -p

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

ENV PASH_TOP=/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
CMD ["/bin/bash"]
