import sys
import re

filename = sys.argv[1]
# print(filename)

with open(filename) as f:
    lines = f.readlines()

no_header = [line for line in lines[1:]
             if line != '\n']

data = "".join(no_header)

# print(data)

re_matcher = re.compile('[0-9]+\.[0-9]+\n')

rest = data
while True:
    temp = rest.split(":")
    if len(temp) == 1:
        break
    bench = temp[0]
    # print("Bench:", bench)
    rest = ":".join(temp[1:])
    perf = re_matcher.findall(rest)[0]
    # perf = re.findall('[0-9]+\.[0-9]+\n', rest)[0]
    # print("Perf:", perf)
    rest = rest.split(perf)[1]
    # print(rest)
    print(f'{bench}:\t{perf.rstrip()}')
    
    # break
