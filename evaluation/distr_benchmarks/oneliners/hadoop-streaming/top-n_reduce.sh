#!/bin/bash
uniq -c | sort -rn | sed 100q
