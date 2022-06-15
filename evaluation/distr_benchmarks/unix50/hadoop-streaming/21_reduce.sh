#!/bin/bash
sed 's/[[:space:]]*$//' | awk "length >= 16"
