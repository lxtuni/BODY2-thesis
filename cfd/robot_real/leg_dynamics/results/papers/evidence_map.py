#!/usr/bin/env python3
"""
四篇论文 × BODY2：三张证据图（力方向利用率 / 腿间干涉增益 / 回收相折叠面积比）
    python3 evidence_map.py [--out .]
数据来源：Qu 2025 Fig 4g、7g；Wang 2025 Table 2/3 与 3.2 节；Fish 1984 p.194；BODY2 verify_analysis.json。输出 fig_evidence_map.png
"""
import os, json, argparse, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__))); a = ap.parse_args()
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
BLUE, ORANGE, PURPLE, RED = "#2a78d6", "#eb6834", "#7a5cc0", "#e34948"        # 分类：BODY2 / Qu / Wang（经 validate_palette 通过）
HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(HERE, "..", "..", "results_verify", "verify_analysis.json")))
b2 = {k: S[k]["Fz_mean"] / S[k]["Fx_mean"] for k in ("current", "V1", "V2", "V3")}

fig, axs = plt.subplots(1, 3, figsize=(16, 5.6), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.25, 1.05, 0.95], wspace=0.36, left=0.125, right=0.985, top=0.80, bottom=0.14))
def style(ax):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8, zorder=0, axis="x")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID); ax.tick_params(colors=INK2, labelsize=9.5, length=0)

# ---------- A：F̄z / F̄x（力方向利用率）----------
ax = axs[0]; style(ax)
rows = [("BODY2 V3  120→70°", b2["V3"], BLUE), ("BODY2 V2", b2["V2"], BLUE), ("BODY2 V1", b2["V1"], BLUE), ("BODY2 现步态", b2["current"], BLUE),
        ("Qu 整机 PP33 实测", 0.53 / 0.92, ORANGE), ("Qu 整机 PP25 实测", 0.63 / 0.86, ORANGE), ("Qu 整机 PP50 实测", 0.66 / 0.69, ORANGE),
        ("Wang 小跑 α0 60° CFD", 3.23 / 3.34, PURPLE), ("Wang 小跑 α0 70° CFD", 3.42 / 2.76, PURPLE),
        ("Qu 单腿 PP50 实测", 0.17 / 0.12, ORANGE), ("Qu 单腿 PP33 实测", 0.23 / 0.15, ORANGE), ("Qu 单腿 PP25 实测", 0.32 / 0.17, ORANGE),
        ("Wang 小跑 α0 80° CFD", 3.97 / 1.83, PURPLE)]
y = np.arange(len(rows))[::-1]
for yi, (nm, v, c) in zip(y, rows):
    ax.barh(yi, v, height=0.62, color=c, alpha=0.9 if c == BLUE else 0.75, zorder=3)
    ax.text(v + 0.03, yi, f"{v:.2f}", va="center", fontsize=9, color=INK2)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
ax.axvline(1.0, color=INK2, lw=1, ls=(0, (4, 3)), zorder=2); ax.text(1.02, len(rows) - 0.55, "垂向力 = 推力", fontsize=8.5, color=INK2, va="top")
ax.set_xlim(0, 2.45); ax.set_xlabel("周期平均 F̄z / F̄x [−]（越小，力越“正”）", color=INK2, fontsize=10)
ax.set_title("A  力方向利用率：平行四边形 + 扇形桨 ≫ 两连杆狗腿\n     BODY2 V3 = 0.28，狗腿机器人 0.6–2.2", loc="left", fontsize=11, color=INK, pad=9)

# ---------- B：整机 / (4 × 单腿) ----------
ax = axs[1]; style(ax); ax.grid(False); ax.grid(color=GRID, lw=0.8, zorder=0, axis="y")
ax.spines["left"].set_visible(True); ax.spines["left"].set_color(GRID); ax.spines["bottom"].set_visible(False)
T3 = {"两相小跑": [1.12, 1.18, 1.09, 1.13, 1.06, 1.02, 1.20, 1.26, 1.19], "四相 DR33": [1.20, 1.42, 1.52, 1.33, 1.37, 1.50, 1.17, 1.44, 1.71], "四相 DR25": [1.00, 1.16, 1.87, 1.26, 1.34, 1.56, 1.04, 1.20, 1.47]}
qu = {"两相 PP50": 0.69 / (4 * 0.12), "四相 PP33": 0.92 / (4 * 0.15), "四相 PP25": 0.86 / (4 * 0.17)}
xs = [0, 1, 2]; rng = np.random.default_rng(3)
for xi, (k, v) in zip(xs, T3.items()):
    ax.scatter(xi + rng.uniform(-0.13, 0.13, len(v)), v, s=42, color=PURPLE, alpha=0.75, zorder=3, edgecolors=SURF, linewidths=1)
    ax.plot([xi - 0.22, xi + 0.22], [np.mean(v)] * 2, color=PURPLE, lw=2.2, zorder=4)
    ax.text(xi, 0.815, f"Wang CFD 均值 {np.mean(v):.2f}", ha="center", va="bottom", fontsize=8.5, color=INK2)
for xi, (k, v) in zip(xs, qu.items()):
    ax.scatter([xi + 0.3], [v], s=110, marker="D", color=ORANGE, zorder=5, edgecolors=SURF, linewidths=1.4)
    ax.text(xi + 0.3, v + 0.06, f"Qu 实测\n{v:.2f}", ha="center", va="bottom", fontsize=8.5, color=INK2)
ax.axhline(1.0, color=INK2, lw=1, ls=(0, (4, 3)), zorder=2); ax.text(-0.42, 1.015, "1.0 = 四条独立腿", fontsize=8.5, color=INK2, va="bottom")
ax.set_xticks(xs); ax.set_xticklabels(["两相（对角对）", "四相 DR/PP 33 %", "四相 DR/PP 25 %"]); ax.set_xlim(-0.5, 2.6); ax.set_ylim(0.78, 2.05); ax.tick_params(axis="x", pad=14)
ax.set_ylabel("整机推力 / （4 × 单腿推力）[−]", color=INK2, fontsize=10)
ax.set_title("B  腿间干涉是正的（两个独立来源）\n     → BODY2 按单腿 ×4 估整机是下限", loc="left", fontsize=11, color=INK, pad=9)

# ---------- C：回收相折叠面积比 ----------
ax = axs[2]; style(ax); ax.grid(False); ax.grid(color=GRID, lw=0.8, zorder=0, axis="y")
ax.spines["left"].set_visible(True); ax.spines["left"].set_color(GRID); ax.spines["bottom"].set_visible(False)
STEM, STEM_W, HEAD, FAN_R, FAN = 0.05, 0.003, 0.042, 0.042, 120
A_o = STEM * STEM_W + 0.5 * FAN_R ** 2 * np.radians(FAN)
def ratio(cw): return (STEM * STEM_W + cw * HEAD - (1 - np.pi / 4) * cw * 0.025 / 2) / A_o
items = [("麝鼠\nFish 1984", 1 - 0.555, ORANGE), ("BODY2 CFD\n收拢 30 mm", ratio(0.030), BLUE), ("收拢 20 mm", ratio(0.020), BLUE), ("收拢 15 mm", ratio(0.015), BLUE), ("收拢 10 mm\n（真扇估计 8–12）", ratio(0.010), BLUE)]
for i, (nm, v, c) in enumerate(items):
    ax.bar(i, v, width=0.6, color=c, alpha=0.9 if i == 1 else (0.5 if c == BLUE else 0.8), zorder=3)
    ax.text(i, v + 0.015, f"{v:.2f}", ha="center", fontsize=9.5, color=INK2)
ax.axhline(1 - 0.555, color=ORANGE, lw=1.1, ls=(0, (4, 3)), zorder=2)
ax.set_xticks(range(len(items))); ax.set_xticklabels([it[0] for it in items], fontsize=8.8); ax.set_ylim(0, 0.8)
ax.set_ylabel("回收相迎流面积 / 划水相面积 [−]", color=INK2, fontsize=10)
ax.set_title("C  回收折叠：CFD 用的 30 mm 收拢宽偏保守\n     麝鼠 0.45；量一下真扇，跑 1 个 closed 算例", loc="left", fontsize=11, color=INK, pad=9)

fig.suptitle("四篇论文 × BODY2：三条最有用的定量对照", fontsize=13.5, x=0.125, ha="left", y=0.965, color=INK)
fig.text(0.125, 0.905, "A：Qu 2025 Fig 4g/7g（1.2 s 周期实测）、Wang 2025 §3.2（β0 = 20°，小跑）、BODY2 verify_analysis.json；B：Wang Table 3、Qu Fig 4g/7g；C：Fish 1984 p.194，BODY2 桨几何（扇 120°、R 42 mm、桨头 42 mm）",
         fontsize=9, color=MUTED, ha="left")
fn = os.path.join(a.out, "fig_evidence_map.png"); plt.savefig(fn, dpi=120, facecolor=SURF); print("→", fn)
