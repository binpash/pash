#!/bin/sh
# Automatically generated file
# Source file example/NMRPipe.sh
#!/usr/bin/env sgsh
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
set -e 
IN=${FID:-$PASH_TOP/evaluation/benchmarks/dgsh/input/fid}

export 
/opt/nmrpipe/nmrbin.linux212_64/var2pipe -in $IN						\
 -xN            1280            -yN     256		\
 -xT            640             -yT     128		\
 -xMODE         Complex -yMODE  Complex			\
 -xSW           8000    -ySW    6000			\
 -xOBS          599.4489584     -yOBS   60.7485301      \
 -xCAR          4.73    -yCAR   118.000			\
 -xLAB          1H      -yLAB   15N			\
 -ndim          2       -aq2D   States			\
-verb  |
(

	export SGDIR=/tmp/sg-$$.0

	rm -rf $SGDIR

	# Cleanup on exit or interrupt
	cleanup()
	{
		SIGNAL=$1
		[ $SIGNAL = EXIT ] || echo sgsh interrupted. Cleaning up... 1>&2

		# Stop key-value stores
		
		# Kill processes we have launched in the background
		kill $SGPID 2>/dev/null

		# Remove temporary directory
		rm -rf "$SGDIR"

		# Propagate real signals and exit with non-0
		if [ $SIGNAL != EXIT ]
		then
			trap - $SIGNAL EXIT
			kill -s $SIGNAL $$
		fi

		# Exit with the original exit value
		exit

	}

	for sig in HUP INT QUIT TERM EXIT
	do
		trap "cleanup $sig" $sig
	done

	mkdir $SGDIR
	cat <&3 3<&-  >$SGDIR/npi-0.0.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.1.0
 {  /opt/nmrpipe/nmrbin.linux212_64/nmrPipe |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn SOL |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn ZF -auto |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn FT |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn PS -p0 177 -p1 0.0 -di |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn EXT -left -sw -verb |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn TP |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn COADD -cList 1 0 -time |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5 |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn ZF -auto |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn FT |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn PS -p0 0 -p1 0 -di |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn TP |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn POLY -auto -verb > A
: ; } <$SGDIR/npi-0.0.0 >$SGDIR/npfo-none-0.0.0
 {  /opt/nmrpipe/nmrbin.linux212_64/nmrPipe |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn SOL |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn ZF -auto |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn FT |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn PS -p0 177 -p1 0.0 -di |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn EXT -left -sw -verb |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn TP |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn COADD -cList 0 1 -time |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5 |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn ZF -auto |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn FT |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn PS -p0 -90 -p1 0 -di |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn TP |
	   /opt/nmrpipe/nmrbin.linux212_64/nmrPipe -fn POLY -auto -verb > B
: ; } <$SGDIR/npi-0.1.0 >$SGDIR/npfo-none-0.1.0

# Gather the results
./sgsh-tee  -i $SGDIR/npfo-none-0.0.0  -i $SGDIR/npfo-none-0.1.0 > /dev/null

)  3<&0 

# We use temporary files rather than streams, because
# addNMR mmaps its input files. The diagram displayed in the
# example shows the notional data flow.
/opt/nmrpipe/nmrbin.linux212_64/addNMR -in1 A -in2 B -out A+B.sgsh.ft2 -c1 1.0 -c2 1.25 -add 
/opt/nmrpipe/nmrbin.linux212_64/addNMR -in1 A -in2 B -out A-B.sgsh.ft2 -c1 1.0 -c2 1.25 -sub 
rm -f A B A+B.sgsh.ft2 A-B.sgsh.ft2
