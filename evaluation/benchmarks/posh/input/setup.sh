if [[ "$1" == "-c" ]]; then
  rm -rf input
  rm -rf output
  exit 1
fi
mkdir output
./install-deps.sh
if [[ ! -d cr_data ]]; then
  echo "Downloading 37GB"
  wget https://sing.stanford.edu/deeptir/posh_data/posh_data.tar.gz
  echo "extracting 140GB"
  tar -xf posh_data.tar.gz
fi

