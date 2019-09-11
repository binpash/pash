#!/bin/bash

echo "Hello $USER"
echo "$(uptime)" >>"$(date)".txt
echo "Your File is being saved to $(pwd)"
