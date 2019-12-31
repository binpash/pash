#!/bin/bash

# Schematic for distributing computations

A="158.130.4.119"  #TOBLIND
B="158.130.4.113"  #TOBLIND
C="158.130.4.114"  #TOBLIND
D="158.130.4.120"  #TOBLIND
DSTAR="158.130.4.212" #TOBLIND

# Collect results
# Implement: `socat` can listen for multiple connections
nc -l -p 5555 > r1 &
nc -l -p 5556 > r2 &
nc -l -p 5557 > r3 &
nc -l -p 5558 > r4 &

# Receiver should run
nc -l -p 5000 | process | nc $DSTAR 5555

# Distribute load
cat ./a/b | tr 'x' 'x' | nc $B 5000
cat ./a/b | tr 'x' 'x' | nc $C 5000
cat ./a/b | tr 'x' 'x' | nc $D 5000
