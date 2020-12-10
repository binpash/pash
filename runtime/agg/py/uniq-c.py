#!/usr/bin/python
import sys, os, functools, utils

PAD_LEN = 7 # FIXME OSX vs Linux

def parseLine(s):
  global PAD_LEN
  # FIXME: This could identify padding number
  res = s.split(None, 1)
  # print(res)
  return (int(res[0]), res[1])

def emitLine(t):
  global PAD_LEN
  return " ".join([str(t[0]).rjust(PAD_LEN, ' '), t[1]])

def combinePair(a, b):
  # print(a, b)
  az = parseLine(a)
  bz = parseLine(b)
  if (az[1] == bz[1]):
    return [emitLine([az[0] + bz[0], az[1]])]
  else:
    return [a, b] # already strings, no need to emit

def agg(a, b):
  # print(a, b)
  if not a:
    return b
  return a[:-1]  + combinePair(a[-1], b[0]) + b[1:]

utils.help()
res = functools.reduce(agg, utils.read_all(), [])
utils.out("".join(res))
