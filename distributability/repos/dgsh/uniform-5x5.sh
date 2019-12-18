#!/usr/bin/env dgsh

row()
{
  dgsh-parallel -n 5 'echo C{}' | paste
}

matrix()
{
  dgsh-parallel -n 5 row
}

export -f row

matrix | cat
