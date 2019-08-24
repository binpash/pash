#!/usr/bin/env dgsh
#
# SYNOPSIS Nuclear magnetic resonance processing
# DESCRIPTION
# Nuclear magnetic resonance in-phase/anti-phase channel conversion and
# processing in heteronuclear single quantum coherence spectroscopy.
# Demonstrate processing of NMR data using the NMRPipe family of programs.
#
# See also F. Delaglio, S. Grzesiek, G. W. Vuister, G. Zhu, J. Pfeifer
# and A. Bax: NMRPipe: a multidimensional spectral processing system based
# on UNIX pipes. J. Biomol. NMR. 6, 277-293 (1995).
# http://spin.niddk.nih.gov/NMRPipe/
#

# The conversion is configured for the following file:
# http://www.bmrb.wisc.edu/ftp/pub/bmrb/timedomain/bmr6443/timedomain_data/c13-hsqc/june11-se-6426-CA.fid/fid
var2pipe -in $1            \
 -xN            1280            -yN     256    \
 -xT            640             -yT     128    \
 -xMODE         Complex -yMODE  Complex      \
 -xSW           8000    -ySW    6000      \
 -xOBS          599.4489584     -yOBS   60.7485301      \
 -xCAR          4.73    -yCAR   118.000      \
 -xLAB          1H      -yLAB   15N      \
 -ndim          2       -aq2D   States      \
-verb  |
tee |
{{
  # IP/AP channel conversion
  # See http://tech.groups.yahoo.com/group/nmrpipe/message/389
  nmrPipe |
  nmrPipe -fn SOL |
  nmrPipe -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 |
  nmrPipe -fn ZF -auto |
  nmrPipe -fn FT |
  nmrPipe -fn PS -p0 177 -p1 0.0 -di |
  nmrPipe -fn EXT -left -sw -verb |
  nmrPipe -fn TP |
  nmrPipe -fn COADD -cList 1 0 -time |
  nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5 |
  nmrPipe -fn ZF -auto |
  nmrPipe -fn FT |
  nmrPipe -fn PS -p0 0 -p1 0 -di |
  nmrPipe -fn TP |
  nmrPipe -fn POLY -auto -verb >A

  nmrPipe |
  nmrPipe -fn SOL |
  nmrPipe -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 |
  nmrPipe -fn ZF -auto |
  nmrPipe -fn FT |
  nmrPipe -fn PS -p0 177 -p1 0.0 -di |
  nmrPipe -fn EXT -left -sw -verb |
  nmrPipe -fn TP |
  nmrPipe -fn COADD -cList 0 1 -time |
  nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5 |
  nmrPipe -fn ZF -auto |
  nmrPipe -fn FT |
  nmrPipe -fn PS -p0 -90 -p1 0 -di |
  nmrPipe -fn TP |
  nmrPipe -fn POLY -auto -verb >B

}}

# We use temporary files rather than streams, because
# addNMR mmaps its input files. The diagram displayed in the
# example shows the notional data flow.
if [ -z "${DGSH_DRAW_EXIT}" ]
then
	addNMR -in1 A -in2 B -out A+B.dgsh.ft2 -c1 1.0 -c2 1.25 -add
	addNMR -in1 A -in2 B -out A-B.dgsh.ft2 -c1 1.0 -c2 1.25 -sub
fi
