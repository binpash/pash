if [[ -z "$PASH_TOP" ]]; then
    export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
fi

export IN1=$PASH_TOP/evaluation/intro/input/100M.txt
export IN2=$PASH_TOP/evaluation/intro/input/words

if [ ! -f $IN1 ]; then
    $PASH_TOP/evaluation/intro/input/setup.sh
fi

if [[ $(uname -s) == 'Linux' ]]; then
    ./test-linux.sh
elif [[ $(uname -s) == 'FreeBSD' ]]; then
    ./test-bsd.sh
fi
