#!/bin/bash
cat $1 | awk "{print \$2, \$0}"
