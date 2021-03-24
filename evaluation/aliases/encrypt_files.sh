# compress and encrypt all files in a directory 
tar_encrypt_files()
(
    cd $1
    find . -name "*.rtf" | xargs -I {} sh -c "tar -czf - {} | openssl enc -e -pbkdf2 -out {}.enc; mv {}.enc $2" sh {}
)
tar_encrypt_files input/rtf $PWD/output
