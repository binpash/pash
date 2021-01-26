timestamp=`date +"%Y%m%d%H%M%S"`

tar cf cdash.tar \
    *.c *.h *.sh Makefile

cp -p cdash.tar "cdash-${timestamp}.tar"
