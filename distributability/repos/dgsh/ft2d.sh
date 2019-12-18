#!/usr/bin/env dgsh
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
	vppen size=r vpstyle=n gridnum=2,1 /dev/stdin $*
}

export -f side_by_side_iso

# A broadened pulse and the real part of its 2D Fourier transform
sfspike n1=64 n2=64 d1=1 d2=1 nsp=2 k1=16,17 k2=5,5 mag=16,16 \
	label1='time' label2='space' unit1= unit2= |
sfsmooth rect2=2 |
sfsmooth rect2=2 |
tee |
{{
	sfgrey pclip=100 wanttitle=n

	sffft1 |
	sffft3 axis=2 pad=1 |
	sfreal |
	tee |
	{{
		sfwindow f1=1 | sfreverse which=3
		cat
	}} |
	sfcat axis=1 "<|" |
	sfgrey pclip=100 wanttitle=n label1="1/time" label2="1/space"
}} |
call_with_stdin side_by_side_iso '<|' yscale=1.25 >Fig/ft2dofpulse.vpl

# A simulated air wave and the amplitude of its 2D Fourier transform
sfspike n1=64 d1=1 o1=32 nsp=4 k1=1,2,3,4 mag=1,3,3,1 \
	label1='time' unit1= |
sfspray n=32 d=1 o=0 |
sfput label2=space |
sflmostretch delay=0 v0=-1 |
tee |
{{
	sfwindow f2=1 | sfreverse which=2
	cat
}} |
sfcat axis=2 "<|" |
tee |
{{
	sfgrey pclip=100 wanttitle=n

	sffft1 |
	sffft3 sign=1 |
	tee |
	{{
		sfreal
		sfimag
	}} |
	dgsh-wrap -e sfmath nostdin=y re="<|" im="<|" \
	  output="sqrt(re*re+im*im)" |
	tee |
	{{
		sfwindow f1=1 | sfreverse which=3
		cat
	}} |
	sfcat axis=1 "<|" |
	sfgrey pclip=100 wanttitle=n label1="1/time" label2="1/space"
}} |
call_with_stdin side_by_side_iso '<|' yscale=1.25 >Fig/airwave.vpl

wait
