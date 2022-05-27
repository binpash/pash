IN=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/
mkdir -p ${IN}/deps/
# install dependencies
pkgs='ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump'

if ! dpkg -s $pkgs >/dev/null 2>&1 ; then
    sudo apt-get install $pkgs -y
    echo 'Packages Installed'
fi

if [ ! -d ${IN}/deps/samtools-1.7 ]; then
    cd ${IN}/deps/
    wget https://github.com/samtools/samtools/archive/refs/tags/1.7.zip
    unzip 1.7.zip
    rm 1.7.zip
    cd samtools-1.7
    wget https://github.com/samtools/htslib/archive/refs/tags/1.7.zip
    unzip 1.7.zip
    autoheader            # Build config.h.in (this may generate a warning about
    # AC_CONFIG_SUBDIRS - please ignore it).
    autoconf -Wno-syntax  # Generate the configure script
    ./configure           # Needed for choosing optional functionality
    make
    rm -rf 1.7.zip
    echo 'Samtools installed'
fi

if ! dpkg -s "makedeb-makepkg" >/dev/null 2>&1 ; then
    cd ${IN}/deps/
    wget http://pac-n4.csail.mit.edu:81/pash_data/makedeb.deb
    sudo dpkg -i makedeb.deb
    echo 'Makedeb installed'
fi
