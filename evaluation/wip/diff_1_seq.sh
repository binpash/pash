#!/bin/bash
# https://crashingdaily.wordpress.com/2008/03/06/diff-two-stdout-streams/
# This would work with coreutils.

diff -B <(cat $IN | awk 'BEGIN {srand(); OFMT="%.17f"} {print rand(), $0}' "$@" | sort -k1,1n | cut -d ' ' -f2- | sort | tr [:lower:] [:upper:]) <(cat $IN | awk 'BEGIN {srand(); OFMT="%.17f"} {print rand(), $0}' "$@" | sort -k1,1n | cut -d ' ' -f2- | sort | tr [:upper:] [:lower:] )

