FROM debian:10
SHELL ["/bin/bash", "-c"]
## Necessary to avoid interactive responses in installations
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y && apt-get install -y git

## Install pyenv to install python3.8
RUN apt install -y make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl
RUN curl https://pyenv.run | bash

## Add the necessary pyenv post-installation commands in .bashrc
RUN echo 'export PATH="$HOME/.pyenv/bin:$PATH"' >> /root/.bashrc
RUN echo 'eval "$(pyenv init -)"' >> /root/.bashrc
RUN echo 'eval "$(pyenv virtualenv-init -)"' >> /root/.bashrc

## Note that after this point we need to source /root/.bashrc before executing something with python3.8 in the Dockerfile

## Install and set python 3.8.13 as the global python version
RUN source /root/.bashrc && pyenv install -v 3.8.13 && pyenv global 3.8.13

ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV LC_COLLATE C
ENV LC_ALL en_US.UTF-8

ENV PASH_TOP=/opt/pash
# download PaSh
RUN git clone https://github.com/binpash/pash.git /opt/pash
RUN source /root/.bashrc && bash /opt/pash/scripts/distro-deps.sh -o 
RUN source /root/.bashrc && yes | bash /opt/pash/scripts/setup-pash.sh -o
CMD ["/bin/bash"]
