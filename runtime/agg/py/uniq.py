#!/usr/bin/python
import sys, os, functools, utils

# not used
def parseLine(s):
  return s

# not used
def emitLine(t):
  return t

def combinePair(a, b):
  # print(a, b)
  if (a == b):
    return [a]
  else:
    return [a, b] # already strings, no need to emit

def combineBlackbox(a, b):
  return [utils.execute(utils.cmd(), a + b)]

def combiner(a, b):
  if not a:
    return b
  pair = combinePair(a[-1], b[0])
  # pair = combineBlackbox(a[-1], b[0])
  # print pair
  return a[:-1]  + pair + b[1:]

utils.help()
res = functools.reduce(combiner, utils.read_all(), [])
utils.out("".join(res))
