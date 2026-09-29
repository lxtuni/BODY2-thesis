#!/usr/bin/env python3
"""CFD 算例矩阵图示：七组算例，每组把各变体并排画一个周期的频闪 (stroboscopic) 叠影。
→ fig_cfd_matrix.png / .pdf
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
from cfd_matrix_kin import SCENES, foil_poly, C, U, SERVO_RATE

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
GOLD, RED = "#d4a017", "#c2185b"

NSTROBE = 15
DX   = 0.125        # 每个变体横向占的宽度（= 一个周期）
GAPX = 0.048
YHALF = 0.056       # 竖向半高


def draw_strobe(ax, case, x0):
    taus = np.linspace(0, 1, NSTROBE, endpoint=False)
    xs = x0 + np.linspace(0, DX, NSTROBE, endpoint=False)
    tt = np.linspace(0, 1, 200)
    ax.plot(x0 + np.linspace(0, DX, 200), [case.at(t)["h"] for t in tt],
            color=case.color, lw=0.9, alpha=0.30, zorder=1)
    for tau, x in zip(taus, xs):
        s = case.at(tau)
        P = np.array([x, s["h"]])
        poly, LE, TE = foil_poly(P, s["th"], case.c, case.pivot)
        a = 0.32 + 0.68*(0.5 + 0.5*np.cos(2*np.pi*tau))**0.8
        ax.add_patch(plt.Polygon(poly, closed=True, color=case.color, alpha=a, lw=0, zorder=4))
        if case.le == "round":
            ax.add_patch(Circle(LE, 0.0019, color=case.color, alpha=a, lw=0, zorder=5))
        ax.plot(*P, "o", color=GOLD, ms=2.8, mec="none", zorder=6)


fig = plt.figure(figsize=(18.0, 12.2), facecolor=SURF)
gs = fig.add_gridspec(4, 2, left=0.035, right=0.988, top=0.862, bottom=0.028,
                      hspace=0.62, wspace=0.10)
axes = [fig.add_subplot(gs[i, j]) for i in range(4) for j in range(2)]
WMAX = 4*DX + 3*GAPX     # 最宽的一组（4 个变体）

for sc, ax in zip(SCENES, axes):
    cases = sc["cases"]; n = len(cases)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    span = n*DX + (n - 1)*GAPX
    xoff = (WMAX - span)/2
    ax.add_patch(FancyBboxPatch((xoff - 0.012, -YHALF*0.92), span + 0.024, YHALF*1.84,
                                boxstyle="round,pad=0.004", facecolor=WATER, edgecolor="none", zorder=0))
    for k, cs in enumerate(cases):
        x0 = xoff + k*(DX + GAPX)
        draw_strobe(ax, cs, x0)
        ax.text(x0 + DX/2, YHALF*1.02, cs.label, ha="center", va="bottom",
                fontsize=10.6, color=cs.color, weight="bold")
        marg = SERVO_RATE/cs.rate_pk
        ax.text(x0 + DX/2, -YHALF*1.06,
                f"f {cs.f:.2f} Hz · χ {cs.chi:.2f}\n"
                f"θ0 {np.degrees(cs.th0):.0f}° · 俯仰峰速 {cs.rate_pk:.0f} °/s\n"
                f"舵机余量 {marg:.1f}×",
                ha="center", va="top", fontsize=8.2,
                color=RED if marg < 1.5 else INK2, weight="bold" if marg < 1.5 else "normal")
    ax.add_patch(FancyArrowPatch((0.014, YHALF*0.80), (-0.002, YHALF*0.80), arrowstyle="-|>",
                                 mutation_scale=9, color="#9fc3e6", lw=1.1, zorder=3))
    ax.text(0.018, YHALF*0.80, f"U = {U} m/s", fontsize=8, color="#4f86c6", va="center", zorder=3)
    ax.set_xlim(-0.012, WMAX + 0.012)
    ax.set_ylim(-YHALF*1.95, YHALF*1.95)
    ax.set_aspect("equal")
    ax.set_title(sc["title"] + "\n" + sc["why"], loc="left", fontsize=12.0, color=INK, pad=8,
                 linespacing=1.9)
    ax.title.set_color(INK)

# 第 8 格：读图说明
ax = axes[7]
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
ax.add_patch(FancyBboxPatch((0.01, 0.02), 0.98, 0.96, boxstyle="round,pad=0.012",
                            transform=ax.transAxes, facecolor="#f5f4f1", edgecolor=GRID, lw=1.0))
ax.text(0.045, 0.955, "怎么读这张图", transform=ax.transAxes, fontsize=12.8, color=INK,
        va="top", weight="bold")
ax.text(0.045, 0.855,
    "· 每一小块是一个算例，从左到右是一个拍动周期的 15 个频闪快照；颜色越深 = 越靠近周期起点。\n"
    "· 金色小点 = 俯仰轴（pitch axis）；淡色曲线 = 脚蹼中心的沉浮轨迹；来流从右向左。\n"
    "· χ = 羽化系数（feathering parameter）= θ0 / 诱导角。χ 越小 → 攻角越大。\n"
    "· 全部算例：c = 42 mm、h0/c = 0.83（验证算例 0.75）、U = 0.25 m/s、折扇 AR = 2.87。\n\n"
    "两个在开算之前就暴露的舵机约束：\n"
    "① St = 0.60（纯推力最优点）需 642 °/s，7465W 只剩 1.3× 余量 —— 顶到上限。\n"
    "② 在 St = 0.50 上做任何 α 波形整形都顶满舵机：强制正弦 α 854 °/s（1.0×）、\n"
    "　 方波型 α 973 °/s（0.8×）。要么降 St、要么换更快的舵机，开算前必须定。",
    transform=ax.transAxes, fontsize=8.6, color=INK2, va="top", ha="left", linespacing=1.7)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

fig.text(0.035, 0.972, "BODY2 · 海豚尾鳍式拍动 —— 待验证的 CFD 算例矩阵（按文献定）",
         fontsize=18, color=INK, ha="left")
fig.text(0.035, 0.940,
    "基准工况：弦 c = 42 mm（120° 展开态折扇，AR = 2.87）、h0/c = 0.83、ψ = 90°、俯仰轴 1/4 弦、U = 0.25 m/s。"
    "① 与 ② 是必做项，③–⑦ 按时间排。\n"
    "参数依据：Anderson 1998 (JFM 360:41)、Read 2003 (JFS 17:163)、Hover 2004 (JFS 19:37)、"
    "Fish 1993 (JEB 185:179)、Rohr & Fish 2004 (JEB 207:1633)、Fish et al. 2007 (Anat. Rec. 290:614)。\n"
    "推力、效率、涡结构全部交给 CFD —— 本图只锁运动学，不含任何受力估计。",
    fontsize=9.8, color=MUTED, ha="left", va="top", linespacing=1.7)

for ext in ("png", "pdf"):
    p = os.path.join(HERE, f"fig_cfd_matrix.{ext}")
    fig.savefig(p, dpi=165 if ext == "png" else None, facecolor=SURF)
    print("→", p)
