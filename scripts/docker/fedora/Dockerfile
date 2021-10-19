FROM fedora:35
RUN dnf install -y git sudo curl diffutils python-devel pip python3-setuptools zlib-devel libjpeg-devel gcc g++
RUN python3 -m pip install wheel --upgrade Pillow
# download PaSh
RUN git clone https://github.com/binpash/pash.git /opt/pash
# install PaSh
RUN bash /opt/pash/scripts/install.sh -p
# install to /usr/bin directory
RUN cp /opt/pash/pa.sh /usr/bin/
# add to the bashrc
RUN echo "export PASH_TOP=/opt/pash" >> ~/.bashrc
# export PaSh directory
ENV PASH_TOP=/opt/pash
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
