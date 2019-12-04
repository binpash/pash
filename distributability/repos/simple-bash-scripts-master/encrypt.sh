#!/bin/bash

echo "Welcome, I am ready to encrypt a file/folder for you"
echo "currently I have a limitation, Place me to the same folder, where a file to be encrypted is present"
echo "Enter the Exact File Name with extension"
read file
# decryption command
# gpg -d filename.gpg > filename
gpg -c $file
echo "I have encrypted the file sucessfully..."
echo "Now I will be removing the original file"
rm -rf $file
