import matplotlib.pyplot as plt

# Daemon:  [Serverless Manager] Counter value: 47 at time: 1731966273.1897645

def collect_data(filepath):
    data = []
    with open(filepath, "r") as f:
        for line in f.readlines():
            if "Daemon:  [Serverless Manager] Counter value: " in line:
                counter_str, time_str = line.split("Daemon:  [Serverless Manager] Counter value: ")[1].split(" at time: ")
                data.append((int(counter_str), float(time_str)))
    return data

filepath  = "/Users/yxie91/Documents/pash/evaluation/benchmarks/media-conv/times/to_mp3_heavy.sh__envAWS__mem2048M__sysSplash__w1.time"
data = collect_data(filepath)

x = [d[1]-data[0][1] for d in data]
y = [d[0] for d in data]

plt.plot(x, y)
plt.xlabel("Time (s)")
plt.ylabel("Concurrent Workers")
plt.title(filepath)
plt.show()