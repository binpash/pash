#!/bin/bash
read -p "Enter String Uppercase : " i
o=$(echo "$i" | tr '[:upper:]' '[:lower:]')
echo $o