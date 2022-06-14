#!/bin/bash
cat $1 | sed 's/T..:..:..//' | cut -d ',' -f 1,3
