#!/bin/bash
cut -c 3-3 | uniq | sed s/\$/'0s'/
