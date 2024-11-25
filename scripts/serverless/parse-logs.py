# [2024-11-20 01:07:34.051000] START RequestId: 71ddca67-3d9d-49d3-af4f-7f98bd7400a5 Version: $LATEST
# [2024-11-20 01:07:53.229000] END RequestId: 71ddca67-3d9d-49d3-af4f-7f98bd7400a5
# [2024-11-20 01:07:53.229000] REPORT RequestId: 71ddca67-3d9d-49d3-af4f-7f98bd7400a5	Duration: 19178.04 ms	Billed Duration: 19179 ms	Memory Size: 1769 MB	Max Memory Used: 719 MB	

from datetime import datetime
import matplotlib.pyplot as plt
import os
record = []
dir_path = "/Users/yxie91/Documents/pash/evaluation/benchmarks/oneliners/fix-split-10G-w4-eager-2ary-4GM-10GS-logs/logs"
for file in os.listdir(dir_path):
    with open(f"{dir_path}/{file}", "r") as f:
        for line in f.readlines():
            if "START RequestId" in line:
                time = line.split("]")[0][1:]
                if "." not in time:
                    time += ".000000"
                seconds_since_epoch = datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f").timestamp()
                record.append((seconds_since_epoch, 1))
            if "END RequestId" in line:
                time = line.split("]")[0][1:]
                if "." not in time:
                    time += ".000000"
                seconds_since_epoch = datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f").timestamp()
                record.append((seconds_since_epoch, -1))

record.sort(key=lambda x: x[0])

counter = [1]
for i in range(1, len(record)):
    counter.append(counter[-1] + record[i][1])

time_series = [int(item[0])-int(record[0][0]) for item in record]
for i in range(1, len(time_series)):
    print(time_series[i], counter[i])

plt.plot(time_series, counter)
plt.xlabel("Time")
plt.ylabel("Concurrent Workers")
plt.title(dir_path)
plt.show()
