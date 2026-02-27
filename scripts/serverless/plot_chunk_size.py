import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

covid_csv = "evaluation/benchmarks/covid-mts/benchmark_results/2026-02-26_00-56-36/results.csv"
unix50_csv = "evaluation/benchmarks/unix50/benchmark_results/2026-02-26_01-05-25/results.csv"

covid = pd.read_csv(covid_csv).sort_values("chunk_size_mb")
unix50 = pd.read_csv(unix50_csv).sort_values("chunk_size_mb")

fig, ax = plt.subplots(figsize=(18, 6))

ax.semilogx(unix50["chunk_size_mb"], unix50["wall_time_sec"], marker="o", label="unix50/6.sh (width=64)")
ax.semilogx(covid["chunk_size_mb"], covid["wall_time_sec"], marker="s", label="covid-mts/1.sh (width=16)")

# Annotate minimums
for df, label in [(unix50, "unix50"), (covid, "covid-mts")]:
    idx = df["wall_time_sec"].idxmin()
    x = df.loc[idx, "chunk_size_mb"]
    y = df.loc[idx, "wall_time_sec"]
    ax.annotate(
        f"min {y:.1f}s\n@ {x:.1f} MB",
        xy=(x, y),
        xytext=(x * 1.5, y + 3),
        arrowprops=dict(arrowstyle="->", color="gray"),
        fontsize=9,
    )

all_x = sorted(set(unix50["chunk_size_mb"]).union(covid["chunk_size_mb"]))

ax.set_xscale("log")
ax.set_xticks(all_x)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.1f}"))
plt.xticks(rotation=90, fontsize=8)
ax.set_xlabel("Chunk Size (MB, log scale)")
ax.set_ylabel("Wall Time (s)")
ax.set_title("Chunk Size vs Wall Time (s3approxdynamic mode)")
ax.legend()
ax.grid(True, which="both", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("chunk_size_vs_wall_time.png", dpi=150)
print("Saved chunk_size_vs_wall_time.png")

try:
    plt.show()
except Exception:
    pass
