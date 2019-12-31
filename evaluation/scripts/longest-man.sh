#!/bin/bash

# Find the 10 largest man pages

find /usr/share/man -type f | xargs du -scb | sort -rn | head -n 10
