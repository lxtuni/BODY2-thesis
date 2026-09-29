#!/usr/bin/env python3
"""
Wang 等 (2025) Biomimetics 10(3):148 的参数扫描结果整理 + 与 BODY2 的对照 + 步态时序重叠图
    python3 wang2025_compare.py [--out .]
数据全部录自论文 Table 2 / Table 3（开放获取）；F̄x,avg 与 C_T 为论文的无量纲量。
输出 fig_wang2025_compare.png
"""
import os, argparse, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle

ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__))); a = ap.parse_args()
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False

# 配色：分类色经 validate_palette 六项检查全通过（#2a78d6 / #eb6834 / #7a5cc0）
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
BLUE, ORANGE, PURPLE, RED = "#2a78d6", "#eb6834", "#7a5cc0", "#e34948"
RAMP = ["#8fb9ea", "#2a78d6", "#16416f"]                        # 顺序（幅度）：b0 = 20 / 30 / 40°

A0, B0, DR = [60, 70, 80], [20, 30, 40], [50, 33, 25]
Fx = {50: {20: [0.75, 0.59, 0.43], 30: [0.62, 0.57, 0.47], 40: [0.53, 0.45, 0.40]},
      33: {20: [1.27, 0.89, 0.72], 30: [0.82, 0.70, 0.54], 40: [0.73, 0.56, 0.42]},
      25: {20: [2.35, 1.80, 0.86], 30: [1.31, 1.17, 0.84], 40: [1.22, 1.11, 0.76]}}
CT = {50: {20: [0.24, 0.19, 0.14], 30: [0.35, 0.32, 0.26], 40: [0.41, 0.34, 0.30]},
      33: {20: [0.19, 0.13, 0.11], 30: [0.21, 0.18, 0.14], 40: [0.27, 0.21, 0.16]},
      25: {20: [0.18, 0.14, 0.06], 30: [0.17, 0.15, 0.11], 40: [0.23, 0.21, 0.15]}}
T3 = [((60, 20), 1.12, 1.20, 1.00), ((70, 20), 1.18, 1.42, 1.16), ((80, 20), 1.09, 1.52, 1.87),
      ((60, 30), 1.13, 1.33, 1.26), ((70, 30), 1.06, 1.37, 1.34), ((80, 30), 1.02, 1.50, 1.56),
      ((60, 40), 1.20, 1.17, 1.04), ((70, 40), 1.26, 1.44, 1.20), ((80, 40), 1.19, 1.71, 1.47)]

fig = plt.figure(figsize=(15.4, 10.2), facecolor=SURF)
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.22, left=0.058, right=0.985, top=0.875, bottom=0.058)
def style(ax):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9.5)
def title2(ax, t1, t2):
    ax.set_title(t1 + "\n" + t2, loc="left", fontsize=11, color=INK, pad=9)

# ---------------- A：推力随初始摆角单调下降 ----------------
ax = fig.add_subplot(gs[0, 0]); style(ax)
for k, b in enumerate(B0):
    ax.plot(A0, Fx[50][b], "-o", color=RAMP[k], lw=2.2, ms=8, mec=SURF, mew=2, zorder=3)
    ax.annotate(f"摆幅 b0 = {b}°", (A0[0], Fx[50][b][0]), xytext=(-9, 0), textcoords="offset points",
                fontsize=9.5, color=RAMP[k], va="center", ha="right", fontweight="bold")
ax.annotate("", xy=(79.4, 0.445), xytext=(60.6, 0.735), arrowprops=dict(arrowstyle="->", color=RED, lw=1.4, ls=(0, (4, 3))))
ax.text(70.5, 0.635, "α0 每 +20° → 推力 −43 %", fontsize=10, color=RED, ha="left", va="bottom", rotation=-19)
ax.set_xticks(A0); ax.set_xlim(50, 83); ax.set_ylim(0.33, 0.82)
ax.set_xlabel("初始摆角 α0 [°]（峰值推力时 α = α0 + 10°）", color=INK2, fontsize=10)
ax.set_ylabel("单腿周期平均推力 Fx,avg [−]", color=INK2, fontsize=10)
title2(ax, "A  行程中心角越低，推力越大（DR = 50 %）", "     → 与 BODY2 的 V3（中心角 106° → 95°）方向一致")

# ---------------- B：缩短划水相 → 推力涨、效率跌 ----------------
ax = fig.add_subplot(gs[0, 1]); style(ax)
fm = [np.mean([Fx[d][b][i] for b in B0 for i in range(3)]) for d in DR]
cm = [np.mean([CT[d][b][i] for b in B0 for i in range(3)]) for d in DR]
x = np.arange(3); fi, ci = np.array(fm) / fm[0], np.array(cm) / cm[0]
ax.axhline(1.0, color=MUTED, lw=1, ls=(0, (5, 4)), zorder=1)
ax.plot(x, fi, "-o", color=BLUE, lw=2.4, ms=9, mec=SURF, mew=2, zorder=3)
ax.plot(x, ci, "-s", color=ORANGE, lw=2.4, ms=8.5, mec=SURF, mew=2, zorder=3)
ax.annotate("周期平均推力 Fx,avg", (1, fi[1]), xytext=(0, 15), textcoords="offset points",
            fontsize=10, color=BLUE, ha="center", fontweight="bold")
ax.annotate("推力系数 C_T（≈ 效率）", (1, ci[1]), xytext=(0, -20), textcoords="offset points",
            fontsize=10, color=ORANGE, ha="center", fontweight="bold")
for i in (0, 2):
    ax.annotate(f"×{fi[i]:.2f}", (x[i], fi[i]), xytext=(0, 14), textcoords="offset points", ha="center", fontsize=10, color=INK2)
    ax.annotate(f"×{ci[i]:.2f}", (x[i], ci[i]), xytext=(0, -19), textcoords="offset points", ha="center", fontsize=10, color=INK2)
ax.set_xticks(x); ax.set_xticklabels([f"DR = {d} %" for d in DR]); ax.set_xlim(-0.35, 2.35); ax.set_ylim(0.30, 2.80)
ax.set_ylabel("相对 DR = 50 %（9 个 α0×b0 组合均值）", color=INK2, fontsize=10)
ax.set_xlabel("划水相占比 DR（周期 T 不变 → DR 越小，划水越快）", color=INK2, fontsize=10)
title2(ax, "B  划水相压缩 → 推力 ×2.4、效率 ×0.55", "     → 与 BODY2 的 V1 vs V2（效率 2.84 vs 3.45）同一个取舍")

# ---------------- C：整机 / 4×单腿 ----------------
ax = fig.add_subplot(gs[1, 0]); style(ax)
y = np.arange(len(T3))[::-1]
for j, (col, nm, mk) in enumerate([(PURPLE, "两相 小跑（对角对）", "o"),
                                   (BLUE, "四相 对角序列 DR 33 %", "s"),
                                   (ORANGE, "四相 对角序列 DR 25 %", "^")]):
    ax.scatter([r[j + 1] for r in T3], y, s=80, color=col, marker=mk, zorder=3, edgecolors=SURF, linewidths=1.6, label=nm)
ax.axvline(1.0, color=INK2, lw=1.2, zorder=2)
ax.text(1.012, len(T3) - 0.45, "1.0 = 恰好等于\n4 条独立腿", fontsize=9, color=INK2, va="top")
ax.set_yticks(y); ax.set_yticklabels([f"({k[0]}, {k[1]})" for k, *_ in T3], fontsize=9.5)
ax.set_ylabel("(α0, b0) [°]", color=INK2, fontsize=10)
ax.set_ylim(-0.75, len(T3) - 0.25); ax.set_xlim(0.90, 2.00)
ax.set_xlabel("整机推力 / （4 × 单腿推力）[−]", color=INK2, fontsize=10)
ax.legend(frameon=False, fontsize=9.5, loc="lower right", labelcolor=INK2, handletextpad=0.4)
title2(ax, "C  腿间干涉是正的：整机总比 4 条独立腿多",
       "     → BODY2 用单腿 ×4 估整机偏保守；四相增益更大，但力矩也更大")

# ---------------- D：BODY2 V3（DUTY = 0.45）的两相 / 四相时序 ----------------
ax = fig.add_subplot(gs[1, 1]); style(ax); ax.grid(False)
for s in ("left", "bottom"): ax.spines[s].set_visible(False)
DUTY, H = 0.45, 0.62
rows2 = [("LF + RH", 0.00, 9.10), ("RF + LH", 0.50, 8.25)]          # (名, 相位偏移, y)
rows4 = [("LF", 0.00, 5.05), ("LH", 0.25, 4.20), ("RF", 0.50, 3.35), ("RH", 0.75, 2.50)]
def block(rows, npl, ycurve, title, ytitle, col_note, note):
    for nm, off, yb in rows:
        ax.add_patch(Rectangle((0, yb), 1, H, facecolor="#f1f0ed", edgecolor="none", zorder=1))
        segs = [(off, min(off + DUTY, 1.0))] + ([(0.0, off + DUTY - 1.0)] if off + DUTY > 1.0 else [])
        for s0, s1 in segs:
            ax.add_patch(Rectangle((s0, yb), s1 - s0, H, facecolor=BLUE, alpha=0.88, edgecolor=SURF, lw=1.5, zorder=3))
        ax.text(-0.02, yb + H / 2, nm, ha="right", va="center", fontsize=9.5, color=INK2)
    tt = np.linspace(0, 1, 2000, endpoint=False); n = np.zeros_like(tt)
    for nm, off, _ in rows: n += npl * (((tt - off) % 1.0) < DUTY)
    ax.plot(tt, ycurve + 0.34 * n, drawstyle="steps-post", color=INK, lw=1.7, zorder=4)
    for lv in (0, 1, 2):
        ax.plot([0, 1], [ycurve + 0.34 * lv] * 2, color=GRID, lw=0.8, zorder=1)
        ax.text(-0.02, ycurve + 0.34 * lv, str(lv), ha="right", va="center", fontsize=8.5, color=MUTED)
    ax.text(1.015, ycurve + 0.34, f"同时划水的腿数\n（平均 {n.mean():.1f}）", fontsize=8.5, color=INK2, va="center")
    ax.text(0, ytitle, title, fontsize=10.5, color=INK, fontweight="bold", va="bottom")
    ax.text(0.5, ytitle, note, fontsize=9.5, color=col_note, va="bottom", ha="left")
block(rows2, 2, 6.55, "两相：对角对，偏移 T/2", 9.85, BLUE, "几乎无重叠，仅 2 × 0.05 T 空档")
block(rows4, 1, 0.80, "四相：对角序列，偏移 T/4", 5.80, RED, "相邻腿重叠 0.20 T → Hindering 模式")
for s, e in ((0.0, 0.20), (0.25, 0.45), (0.50, 0.70), (0.75, 0.95)):     # 四相的重叠带
    ax.add_patch(Rectangle((s, 2.42), e - s, 5.05 - 2.42 + H + 0.08, facecolor=RED, alpha=0.11, edgecolor="none", zorder=2))
for s, e in ((0.45, 0.50), (0.95, 1.00)):                                # 两相的空档
    ax.add_patch(Rectangle((s, 8.17), e - s, 9.10 - 8.17 + H + 0.08, facecolor=MUTED, alpha=0.22, edgecolor="none", zorder=2))
ax.set_xlim(-0.13, 1.20); ax.set_ylim(0.35, 10.45)
ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.45, 0.75, 1.0])
ax.set_xticklabels(["0", "0.25", "0.45\nDUTY", "0.75", "1"])
ax.set_xlabel("周期相位 τ（蓝 = 划水相）", color=INK2, fontsize=10)
title2(ax, "D  把 V3（DUTY = 0.45）放进整机会怎样",
       "     → 0.45 恰好铺满周期：两相无重叠，四相必然互相干涉")

fig.suptitle("Wang 等 (2025) Biomimetics 10(3):148 —— 参数扫描结果与 BODY2 的对照", fontsize=14, x=0.058, ha="left", y=0.972, color=INK)
fig.text(0.058, 0.930, "数据录自论文 Table 2（单腿，27 个算例）与 Table 3（整机 ÷ 4×单腿）。论文条件：浸没边界法、机体固定、刚性腿、无自由面、Re ≈ 3×10⁴（按肢宽 46.5 mm 与桨尖均速 0.64 m/s）",
         fontsize=9.5, color=MUTED, ha="left")
fn = os.path.join(a.out, "fig_wang2025_compare.png"); plt.savefig(fn, dpi=120, facecolor=SURF); print("→", fn)
