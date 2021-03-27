#!/bin/bash

wget https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z
7za x wikipedia-en-html.tar.7z
tar -xvf wikipedia-en-html.tar > /dev/null
wget http://ndr.md/data/wikipedia/index.txt

head -n 5 index.txt > 5.txt
head -n 100 index.txt > 100.txt

npm install



