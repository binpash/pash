#!/bin/bash
# A bash script for determining how many pages are in a folder of OpenOffice documents
# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7
echo "$(exiftool *.odt | grep Page-count | cut -d ":" -f2 | tr '\n' '+')""0" | bc
