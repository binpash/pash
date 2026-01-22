#!/usr/bin/env bash

sed 's/</{/g' "$1" | \
    sed 's/>/}/g' | \
    sed 's/(/[/g' | \
    sed 's/)/]/g' | \
    python -m json.tool
