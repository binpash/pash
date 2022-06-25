#!/bin/sh
func_emit_tests_Makefile_am ()
{
  ofd=3
  {
    echo hi
  } >&$ofd
}
fd=1 
echo hi >&$fd
