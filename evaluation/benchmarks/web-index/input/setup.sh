#!/bin/bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

# another solution for capturing HTTP status code
# https://superuser.com/a/590170

if [[ "$1" == "-c" ]]; then
  rm -f wikipedia-en-html.tar.7z wikipedia-en-html.tar.7z
  rm -rf en/
  rm -rf npm_modules/
  exit
fi

## FIXME Isn't this input huge? If so let's warn the users before running it, or maybe require a --full flag.
if [ ! -d ./en/ ]; then
  wget https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z
  7za x wikipedia-en-html.tar.7z
  tar -xvf wikipedia-en-html.tar > /dev/null
fi

if [ ! -f index.txt ]; then
  wget http://ndr.md/data/wikipedia/index.txt
  head -n 5 index.txt > 5.txt
  head -n 100 index.txt > 100.txt
fi

if [ ! -d ./node_modules/ ]; then
  npm install
fi



