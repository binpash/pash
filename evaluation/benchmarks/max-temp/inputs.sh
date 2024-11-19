cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

if [ ! -f ./temperatures.2015.txt ]; then
    bash ../max-temp-preprocess.sh
fi

if [ ! -f ./temperatures_1G.txt  ]; then
    touch temperatures_1G.txt

    tail -n 2000000 temperatures.2015.txt >>temperatures_1G.txt
    head -n 2000000 temperatures.2015.txt >>temperatures_1G.txt

    "$PASH_TOP/scripts/append_nl_if_not.sh" temperatures_1G.txt
fi