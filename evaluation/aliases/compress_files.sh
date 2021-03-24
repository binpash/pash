# compress all the files in a directory using dd and tar
compress_files()
(
    cd $1
    find . -name "*.rtf" | xargs -P16 -I {} sh -c "dd if={} bs=1 status=none > '{}f'; tar -zcf {}.tar.gz {}f; rm {}f; mv {}.tar.gz $2" sh {}
)

compress_files input/rtf $PWD/output
