# install dependencies
pkgs='ffmpeg unrtf imagemagick libarchive-tools zstd liblzma-dev libbz2-dev'
sudo apt-get install $pkgs -y

wget https://github.com/samtools/samtools/archive/refs/tags/1.7.zip
unzip 1.7.zip
cd samtools-1.7
wget https://github.com/samtools/htslib/archive/refs/tags/1.7.zip
unzip 1.7.zip
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make

wget http://pac-n4.csail.mit.edu:81/pash_data/makedeb.deb
sudo dpkg -i makedeb.deb
