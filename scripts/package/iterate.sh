#! /usr/bin/env bash
set -u
./deploy.sh "$1" 0.0.1 "$2"
read -p "Repeat? [Y/n]" -n 1 -r
echo
[[ $REPLY =~ ^[Yy]$ ]] && "$0"
