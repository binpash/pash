#!/usr/bin/python
import sys, os, functools, utils

def agg(a, b):
  return a + b

utils.help()
utils.out("".join(functools.reduce(agg, utils.read_all(), [])))
