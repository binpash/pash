rm -rf logs && mkdir logs
for i in {1..84};do 
    for j in real_logs/*;do
        n=$(basename $j)
        cat $j > logs/log${i}_${n}.log; 
    done
done