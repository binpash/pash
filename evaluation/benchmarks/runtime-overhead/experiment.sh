rm -f total-times.txt
touch total-times.txt

for i in `seq 10`;
do
    ./run.sh | tee -a total-times.txt
done
