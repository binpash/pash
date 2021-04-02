if [[ $1 == "-c" ]]; then
    rm -rf input/
    exit 1
fi
mkdir -p input
cd input
if [[ ! -f dblp.xml ]]; then
    wget https://dblp.uni-trier.de/xml/dblp.xml.gz
    gunzip dblp.xml.gz
    cat dblp.xml | head -n 35 > mini.xml
fi

if [[ ! -f fid ]]; then
    wget http://www.bmrb.wisc.edu/ftp/pub/bmrb/timedomain/bmr6443/timedomain_data/c13-hsqc/june11-se-6426-CA.fid/fid
    wget https://www.stats.govt.nz/assets/Uploads/International-trade/International-trade-December-2020-quarter/Download-data/international-trade-december-2020-quarter-csv.zip
    unzip *.zip
    mv final/* .
    rm -rf final
fi
#!/bin/sh

# MIT License

# Copyright (c) 2017 Stockholm University

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

install_dir=/opt/nmrpipe
mkdir -p /opt/nmrpipe
mkdir -p nmr
download_dir=$PWD/nmr
if [[ ! -f $download_dir/install.com ]]; then
    URL_list="https://www.ibbr.umd.edu/nmrpipe/install.com \
        https://www.ibbr.umd.edu/nmrpipe/binval.com \
        https://www.ibbr.umd.edu/nmrpipe/NMRPipeX.tZ \
        https://www.ibbr.umd.edu/nmrpipe/s.tZ \
        https://www.ibbr.umd.edu/nmrpipe/dyn.tZ \
        https://www.ibbr.umd.edu/nmrpipe/talos.tZ \
        http://spin.niddk.nih.gov/bax/software/smile/plugin.smile.tZ"

    for i in $URL_list; do
        echo downloading $i
        wget "--directory-prefix=$download_dir" --no-directories $i
    done
    sudo sh -c "cd $install_dir ; ./install.com +src $download_dir option +type linux212_64"
    sudo mkdir "$install_dir"
    sudo find $download_dir -maxdepth 1 -mindepth 1 -type f -name '*.com' -exec cp {} "$install_dir" \;
    sudo find "$install_dir" -maxdepth 1 -mindepth 1 -type f -name '*.com' -exec chmod 755 {} \;

    ../install-deps.sh
fi

