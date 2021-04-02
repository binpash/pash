#!/bin/bash
# tag: nuclear
# from: https://www2.dmst.aueb.gr/dds/sw/dgsh/#NMRPipe
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input/fid}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}
rm -f a b  A B
mkfifo a b A B

# IP/AP channel conversion
# See http://tech.groups.yahoo.com/group//opt/bin/opt/nmrbin.linux9/nmrPipe/message/389
cat a | /opt/nmrpipe/nmrbin.linux9/nmrPipe |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn SOL |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn ZF -auto |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn FT |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn PS -p0 177 -p1 0.0 -di |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn EXT -left -sw -verb |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn TP |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn COADD -cList 1 0 -time |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5 |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn ZF -auto |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn FT |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn PS -p0 0 -p1 0 -di |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn TP |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn POLY -auto -verb >A &

cat b | /opt/nmrpipe/nmrbin.linux9/nmrPipe |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn SOL |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn ZF -auto |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn FT |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn PS -p0 177 -p1 0.0 -di |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn EXT -left -sw -verb |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn TP |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn COADD -cList 0 1 -time |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5 |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn ZF -auto |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn FT |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn PS -p0 -90 -p1 0 -di |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn TP |
    /opt/nmrpipe/nmrbin.linux9/nmrPipe -fn POLY -auto -verb >B & 


/opt/nmrpipe/nmrbin.linux9/var2pipe -in ${IN}            \
    -xN            1280            -yN     256    \
    -xT            640             -yT     128    \
    -xMODE         Complex -yMODE  Complex      \
    -xSW           8000    -ySW    6000      \
    -xOBS          599.4489584     -yOBS   60.7485301      \
    -xCAR          4.73    -yCAR   118.000      \
    -xLAB          1H      -yLAB   15N      \
    -ndim          2       -aq2D   States      \
    -verb  | tee a b > /dev/null 
#FIXME

# We use temporary files rather than streams, because
# addNMR mmaps its IN files. The diagram displayed in the
# example shows the notional data flow.
if [ -z "${DGSH_DRAW_EXIT}" ]
then
    /opt/nmrpipe/nmrbin.linux9/addNMR -in1 A -in2 B -out A+B.dgsh.ft2 -c1 1.0 -c2 1.25 -add
    /opt/nmrpipe/nmrbin.linux9/addNMR -in1 A -in2 B -out A-B.dgsh.ft2 -c1 1.0 -c2 1.25 -sub
fi
rm -f a b A B
