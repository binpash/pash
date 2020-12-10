#!/usr/bin/python
import sys, os, functools, utils

PAD_LEN = 7 # needs to add a space character for when they exceed

def parseLine(s):
  global PAD_LEN
  # FIXME: This could identify padding number
  return map(int, s.split())

def emitLine(t):
  global PAD_LEN
  return [" ".join(map(lambda e: str(e).rjust(PAD_LEN, ' '), t))]

def combiner(a, b):
  # print(a, b)
  if not a:
    return b
  az = parseLine(a[0])
  bz = parseLine(b[0])
  return emitLine([ (i+j) for (i,j) in zip(az, bz) ])

utils.help()
res = functools.reduce(combiner, utils.read_all(), [])
utils.out("".join(res))
