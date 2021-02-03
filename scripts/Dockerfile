FROM ubuntu:18.04


RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y git sudo locales locales-all curl

# NVM
ENV NVM_DIR ~/.nvm
ENV NODE_VERSION 8.9.4 # Reportedly the lowest NodeJS version used by stakeholders.
RUN curl https://raw.githubusercontent.com/creationix/nvm/v0.37.2/install.sh | bash \
    && . $NVM_DIR/nvm.sh \
    && nvm install $NODE_VERSION \
    && nvm alias default $NODE_VERSION \
    && nvm use default

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
