#!/bin/bash
cat $1 | sed 's/T..:..:..//' | cut -d ',' -f 3,1
