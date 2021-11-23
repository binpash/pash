rm -rf pcaps && mkdir pcaps
# generates 20G
for i in {1..15};do
    for j in input/*;do
        n=$(basename $j)
        cat $j > pcaps/pcap${i}_${n}; 
    
    done
done
