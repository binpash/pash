# install dependencies
pkgs='ffmpeg unrtf imagemagick'
if ! dpkg -s $pkgs >/dev/null 2>&1; then
  sudo apt-get install $pkgs -y
fi

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

sudo apt-get install libarchive-tools zstd
wget https://proget.hunterwittenborn.com/debian-feeds/makedeb/main/makedeb-makepkg/all/makedeb-makepkg_all_8.8.0-1.deb
sudo dpkg -i makedeb-makepkg_all_8.8.0-1.deb
