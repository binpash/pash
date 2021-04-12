#!/bin/sh
# Automatically generated file
# Source file example/ft2d.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Waves: 2D Fourier transforms
# DESCRIPTION
# Create two graphs:
# 1) a broadened pulse and the real part of its 2D Fourier transform, and
# 2) a simulated air wave and the amplitude of its 2D Fourier transform.
# Demonstrates using the tools of the Madagascar shared research environment
# for computational data analysis in geophysics and related fields.
# Also demonstrates the use of two scatter blocks in the same script,
# and the used of named streams.
#
# Adapted from: http://www.reproducibility.org/RSF/book/bei/ft1/ft2d.html
# Description: http://www.reproducibility.org/RSF/book/bei/ft1/paper_html/node14.html
# Madagascar project: http://www.reproducibility.org
#

mkdir -p Fig

# The SConstruct SideBySideIso "Result" method
side_by_side_iso()
{
	vppen size=r vpstyle=n gridnum=2,1 $*
}

# A broadened pulse and the real part of its 2D Fourier transform
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
	 {  sfspike n1=64 n2=64 d1=1 d2=1 nsp=2 k1=16,17 k2=5,5 mag=16,16 \
        label1='time' label2='space' unit1= unit2= |
        sfsmooth rect2=2 |
	sfsmooth rect2=2
}</dev/null  >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
 {  sfgrey pclip=100 wanttitle=n
} <$SGDIR/npi-1.0.0 >$SGDIR/npfo-pulse.vpl.0
 {  sffft1 | sffft3 axis=2 pad=1 | sfreal
} <$SGDIR/npi-1.1.0 >$SGDIR/npi-2.0.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.1.0
 {  sgsh-tee -I
} <$SGDIR/npi-2.0.0 >$SGDIR/npfo-ft2d.0
 {  sfwindow f1=1 |
			    sfreverse which=3 |
			    sfcat axis=1 $SGDIR/npfo-ft2d.0 |
			    sfgrey pclip=100 wanttitle=n \
			     label1="1/time" label2="1/space"
} <$SGDIR/npi-2.1.0 >$SGDIR/npfo-ft2d.vpl.0
 {  side_by_side_iso $SGDIR/npfo-pulse.vpl.0 $SGDIR/npfo-ft2d.vpl.0 \
	   yscale=1.25 >Fig/ft2dofpulse.vpl
: ; }</dev/null  >$SGDIR/npfo-none-0.1.0

# Gather the results
sgsh-tee  -i $SGDIR/npfo-none-0.1.0 >/dev/null

)  3<&0 

# A simulated air wave and the amplitude of its 2D Fourier transform
(

	export SGDIR=/tmp/sg-$$.1

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
	 {  sfspike n1=64 d1=1 o1=32 nsp=4 k1=1,2,3,4 mag=1,3,3,1 \
		label1='time' unit1= |
	   sfspray n=32 d=1 o=0 |
	   sfput label2=space |
	   sflmostretch delay=0 v0=-1
}</dev/null  >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
 {  sgsh-tee -I
} <$SGDIR/npi-1.0.0 >$SGDIR/npfo-air.0
 {  sfwindow f2=1 |
		   sfreverse which=2 |
		   sfcat axis=2 $SGDIR/npfo-air.0
} <$SGDIR/npi-1.1.0 >$SGDIR/npi-2.0.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.1.0
 {  sfgrey pclip=100 wanttitle=n
} <$SGDIR/npi-2.0.0 >$SGDIR/npfo-airtx.vpl.0
 {  sffft1 |
			   sffft3 sign=1
} <$SGDIR/npi-2.1.0 >$SGDIR/npi-3.0.0
ln $SGDIR/npi-3.0.0 $SGDIR/npi-3.1.0
 {  sfreal
} <$SGDIR/npi-3.0.0 >$SGDIR/npfo-airftr.0
 {  sfimag
} <$SGDIR/npi-3.1.0 >$SGDIR/npfo-airfti.0
 {  sfmath re=$SGDIR/npfo-airftr.0 im=$SGDIR/npfo-airfti.0 output="sqrt(re*re+im*im)"
}</dev/null  >$SGDIR/npi-4.0.0
ln $SGDIR/npi-4.0.0 $SGDIR/npi-4.1.0
 {  sgsh-tee -I
} <$SGDIR/npi-4.0.0 >$SGDIR/npfo-airft1.0
 {  sfwindow f1=1 |
		   sfreverse which=3 |
		   sfcat axis=1 $SGDIR/npfo-airft1.0 |
		   sfgrey pclip=100 wanttitle=n label1="1/time" \
		   	label2="1/space"
} <$SGDIR/npi-4.1.0 >$SGDIR/npfo-airfk.vpl.0
 {  side_by_side_iso $SGDIR/npfo-airtx.vpl.0 $SGDIR/npfo-airfk.vpl.0 \
	   yscale=1.25 >Fig/airwave.vpl
: ; }</dev/null  >$SGDIR/npfo-none-0.2.0

# Gather the results
sgsh-tee  -i $SGDIR/npfo-none-0.2.0 >/dev/null

)  3<&0 
