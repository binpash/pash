#!/usr/bin/env bash

## This command drains a stream. It is used if we want a prefix of a
## stream that was written by tee. Since tee writes in both streams
## "almost" in lockstep, if we get a prefix on one side, the other
## side cannot progress.
dd of=/dev/null > /dev/null 2>&1
