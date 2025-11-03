#!/usr/bin/env python3
import re
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

# === 1) Parse input directory ===
parser = argparse.ArgumentParser(description="Extract and plot Lambda billed durations from log files.")
parser.add_argument("directory", help="Path to directory containing log files")
parser.add_argument("--pattern", default="*.log", help="Glob pattern for log files (default: *.log)")
args = parser.parse_args()

log_dir = args.directory
pattern = args.pattern

if not os.path.isdir(log_dir):
    raise NotADirectoryError(f"{log_dir} is not a valid directory")

# === 2) Gather log files ===
log_files = glob.glob(os.path.join(log_dir, pattern))

if not log_files:
    raise FileNotFoundError(f"No files matching {pattern} found in {log_dir}")

durations = []

# === 3) Regex for billed duration ===
pattern = re.compile(r"Billed Duration:\s*(\d+)\s*ms")

for f in log_files:
    with open(f, "r", encoding="utf-8") as fh:
        for line in fh:
            m = pattern.search(line)
            if m:
                durations.append(int(m.group(1)))

if not durations:
    raise ValueError("No Billed Duration entries found!")

# === 4) Convert to DataFrame ===
df = pd.DataFrame(durations, columns=["billed_ms"])

# === 5) Plot distribution ===
plt.figure(figsize=(8, 5))
plt.hist(df["billed_ms"], bins=40, edgecolor="black", alpha=0.7)
plt.xlabel("Billed Duration (ms)")
plt.ylabel("Count")
plt.title("Distribution of Lambda Billed Durations")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("billed_duration_distribution.png")
print("Plot saved as billed_duration_distribution.png")