cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

if [ ! -f ./temperatures.2015.txt ]; then
    bash ../max-temp-preprocess.sh
fi

# if [ ! -f ./temperatures_200M.txt  ]; then
#     touch temperatures_200M.txt

#     tail -n 1200000 temperatures.2015.txt >>temperatures_200M.txt
#     head -n 1200000 temperatures.2015.txt >>temperatures_200M.txt

#     "$PASH_TOP/scripts/append_nl_if_not.sh" temperatures_200M.txt
# fi


# if [ ! -f ./temperatures_500M.txt  ]; then
#     touch temperatures_500M.txt

#     tail -n 2000000 temperatures.2015.txt >>temperatures_500M.txt
#     head -n 2000000 temperatures.2015.txt >>temperatures_500M.txt

#     "$PASH_TOP/scripts/append_nl_if_not.sh" temperatures_500M.txt
# fi

if [ ! -f ./temperatures_1G.txt  ]; then
    touch temperatures_1G.txt

    tail -n 2000000 temperatures.2015.txt >>temperatures_1G.txt
    head -n 2000000 temperatures.2015.txt >>temperatures_1G.txt

    "$PASH_TOP/scripts/append_nl_if_not.sh" temperatures_1G.txt
fi
