#!/usr/bin/python
import sys, os, functools, utils

PAD_LEN = 6

def parseLine(s):
  global PAD_LEN
  # FIXME: This could identify padding number---out of band
  res = s.split(None, 1) # if s else (None, s)
  # if can't parse, just wrap "s"
  return (int(res[0]), res[1]) if len(res) > 0 else (None, s)

def emitLine(tup):
  global PAD_LEN
  return (str(tup[0]).rjust(PAD_LEN, ' ') + '\t' + tup[1]) if tup[0] else tup[1]

def augment(total, b):
    def augment_aux(s):
      (n, l) = parseLine(s)
      # print("augmenting:", n, l)
      t = (total + n, l) if n else (None, l)
      return emitLine(t)
    return map(augment_aux, b)

# There's a chance the last line is empty, ie, with no counter; this searches
# recursively to find a line that has a counter.
def lastCount(a):
    # print("lastCount of", a)
    if not a:
      # print("twice", a)
      return 0
    # print a[-1]
    (n, l) = parseLine(a[-1])
    return n if n else lastCount(a[:-1])

def agg(a, b):
  # print(a, b)
  if not a:
    return b
  # print("combining", a, b)
  return a + augment(lastCount(a), b)

utils.help()
res = functools.reduce(agg, utils.read_all(), [])
utils.out("".join(res))

