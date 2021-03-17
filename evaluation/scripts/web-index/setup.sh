#!/bin/bash

mkdir input
wget https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z
7za x wikipedia-en-html.tar.7z
tar -xvf wikipedia-en-html.tar
wget http://ndr.md/data/wikipedia/index.txt
