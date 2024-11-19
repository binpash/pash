N=40000000 # 40M

mkdir -p inputs
./gen_data.py "$N" > inputs/logs-popcount-org.txt
