myfunction() {
    env | sort > tmp1.txt
}
shellvar1=123456
shellvar2="This is several words"
shellvar3="                        xxx                  "
export shellvar2
trap myfunction EXIT
env | sort > tmp2.txt
