#!/bin/bash
sort -n | head -15 | sed -r 's/^[0-9]+//g'
