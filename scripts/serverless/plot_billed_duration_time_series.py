#!/usr/bin/env python3
import re
import glob
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

# ---- 1) Args ----
parser = argparse.ArgumentParser(
    description="Extract Lambda Billed Duration with timestamps and plot a scatter time series (PDF only)."
)
parser.add_argument("directory", help="Directory containing log files")
parser.add_argument("--pattern", default="*.log", help="Glob for log files (default: *.log)")
parser.add_argument("--outfile", default="billed_duration_timeseries.pdf",
                    help="Output PDF filename (default: billed_duration_timeseries.pdf)")
parser.add_argument("--csv", default="billed_duration_timeseries.csv",
                    help="Output CSV filename (default: billed_duration_timeseries.csv)")
parser.add_argument("--annotate-top", type=int, default=0,
                    help="Annotate top-N longest requests (default: 0 = off)")
args = parser.parse_args()

log_dir = args.directory
if not os.path.isdir(log_dir):
    raise NotADirectoryError(f"{log_dir} is not a valid directory")

files = glob.glob(os.path.join(log_dir, args.pattern))
if not files:
    raise FileNotFoundError(f"No files matching {args.pattern} in {log_dir}")

# ---- 2) Regexes ----
ts_re = re.compile(r"^\[(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{6})\]")
bd_re = re.compile(r"Billed Duration:\s*(\d+)\s*ms")

rows = []
for fp in files:
    with open(fp, "r", encoding="utf-8") as fh:
        for line in fh:
            m_bd = bd_re.search(line)
            if not m_bd:
                continue
            m_ts = ts_re.search(line)
            if not m_ts:
                # Skip if the billed line has no timestamp on the same line
                continue
            rows.append((m_ts.group("ts"), int(m_bd.group(1)), os.path.basename(fp)))

if not rows:
    raise ValueError("No timestamped 'Billed Duration' entries found.")

# ---- 3) DataFrame ----
df = pd.DataFrame(rows, columns=["timestamp", "billed_ms", "source_file"])
df["timestamp"] = pd.to_datetime(df["timestamp"], format="%Y-%m-%d %H:%M:%S.%f", utc=False)
df = df.sort_values("timestamp").reset_index(drop=True)

# ---- 4) Scatter plot (no connecting line) ----
plt.figure(figsize=(10, 5))
plt.scatter(df["timestamp"], df["billed_ms"], s=12)  # scatter only
plt.xlabel("Time")
plt.ylabel("Billed Duration (ms)")
plt.title("Lambda Billed Duration Over Time")
plt.grid(alpha=0.3)

ax = plt.gca()
ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d\n%H:%M:%S"))
plt.tight_layout()

# Optional: annotate top-N longest
if args.annotate_top and args.annotate_top > 0:
    top = df.nlargest(args.annotate_top, "billed_ms")
    for _, row in top.iterrows():
        ax.annotate(f"{row['billed_ms']} ms",
                    (row["timestamp"], row["billed_ms"]),
                    xytext=(0, 8),
                    textcoords="offset points",
                    fontsize=8)

# ---- 5) Save outputs ----
out_pdf = os.path.join(log_dir, args.outfile)
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close()

out_csv = os.path.join(log_dir, args.csv)
df.to_csv(out_csv, index=False)

print(f"Wrote scatter PDF: {out_pdf}")
print(f"Wrote CSV: {out_csv}")
print(f"Entries: {len(df)}, Time range: {df['timestamp'].min()} → {df['timestamp'].max()}")
