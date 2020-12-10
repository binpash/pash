#!/usr/bin/python
import sys, os, functools, utils

def agg(a, b):
  return b + a

utils.help()
res = functools.reduce(agg, utils.read_all(), [])
# Python3 syntax? print("".join(res), end=' ')
# Note comma to avoid newline
utils.out("".join(res))
