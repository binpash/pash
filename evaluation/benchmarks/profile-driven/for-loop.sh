mkdir -p temp-out

for i in $(seq 100)
do
	cat $PASH_TOP/README.md | grep pash | grep pash > temp-out/$i.out 
done

echo done
