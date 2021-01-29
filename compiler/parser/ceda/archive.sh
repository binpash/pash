timestamp=`date +"%Y%m%d%H%M%S"`

tar cf cdash.tar \
    *.c *.h *.sh *.py Makefile

cp -p cdash.tar "cdash-${timestamp}.tar"
