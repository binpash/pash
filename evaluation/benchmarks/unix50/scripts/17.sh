#!/bin/bash

# 7.3: all the decades in which a unix version was released
cat $IN | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/
