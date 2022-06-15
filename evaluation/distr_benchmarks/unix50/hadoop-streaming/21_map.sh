#!/bin/bash
cat $1 | tr -c "[a-z][A-Z]" '\n'
