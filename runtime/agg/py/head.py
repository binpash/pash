#!/usr/bin/python
import sys, os, functools, utils

def agg(a, b):
  # print(a, b)
  if not a:
    return b
  return a

utils.help()
res = functools.reduce(agg, utils.read_all(), [])
utils.out("".join(res))
