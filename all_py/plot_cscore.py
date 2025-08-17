import sys, pandas as pd, matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("用法: python plot_cscore.py <bedgraph> [\"自定义标题\"]")
    sys.exit(1)

bedgraph = sys.argv[1]
title    = sys.argv[2] if len(sys.argv) >= 3 else "A/B Compartment (100 kb resolution)"

df = pd.read_csv(bedgraph, sep=r"\s+", header=None, comment="t", engine="python")
df.columns = ["chr", "start", "end", "cscore"]

plt.figure(figsize=(14, 4))
plt.plot(df["start"], df["cscore"], color="blue", linewidth=1)
plt.axhline(0, ls="--", color="gray")
plt.title(title)                       # 使用可变标题
plt.xlabel("Genomic Position")
plt.ylabel("C-score")
plt.tight_layout()
out_png = bedgraph.replace(".bedgraph", ".png")
plt.savefig(out_png, dpi=300)
print(f"已保存: {out_png}")
