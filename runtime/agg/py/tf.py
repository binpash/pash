#!/usr/bin/python
import sys, os, functools, utils
from collections import OrderedDict

PAD_LEN = 4

def parseLine(s):
  global PAD_LEN
  # FIXME: This could identify padding number---out of band?
  res = s.split()
  # print(res)
  return (res[0], int(res[1])) #, res[2]) 

def emitLine(t):
  global PAD_LEN
  return " ".join(t)

def update_index(index, lst):
  # { index[i]:(index[i] if index[i] else j) for (i, j) in lst}
  # print(index)
  # print(lst)
  print(index, lst)
  index.update(lst)
  return index # {index[i]: j for i, j in lst}

def agg(a, b):
  # print(a, b)
  # TODO: take and emit out of fold!
  return update_index(a, map(parseLine, b));

utils.help()
res = functools.reduce(agg, utils.read_all(), OrderedDict())
utils.out("\n".join([ a + "  " + str(b) for a,b in res.items()]))

