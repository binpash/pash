#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: ${0##*/} program target"
    exit 2
fi

p=$1
tgt=$2

orig=$(${p} ${tgt})
if [ "$?" -ne 0 ];
then echo ${tgt} FAILED, couldn\'t run; exit 2
fi

rt=$(${p} ${tgt} | ${p})
if [ "$?" -ne 0 ];
then echo ${tgt} FAILED round trip, couldn\'t run; exit 3
fi

if [ "${orig}" = "${rt}" ];
then echo ${tgt} OK; exit 0
else
    echo ${tgt} FAILED
    echo ${orig}
    echo ==========
    echo ${rt}
    exit 1
fi
