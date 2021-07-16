FROM fedora:35

# ARG CACHEBUST=1

RUN dnf install -y git sudo curl diffutils python-devel pip python3-setuptools zlib-devel libjpeg-devel gcc g++
RUN python3 -m pip install wheel --upgrade Pillow

RUN git clone https://github.com/andromeda/pash.git 
RUN sudo bash /pash/scripts/install.sh -p 

ENV PASH_TOP=/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

CMD ["/bin/bash"]
