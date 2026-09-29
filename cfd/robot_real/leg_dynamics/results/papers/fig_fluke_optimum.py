#!/usr/bin/env python3
"""最优脚蹼形态 + 最佳尾鳍轨迹 总览图（读 sweep_A/B/C.csv）→ fig_fluke_optimum.png/.pdf"""
import os, json, math, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon, FancyArrowPatch
from scipy.interpolate import griddata

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
BLUE, RED, GOLD, GREEN, PUR, TEAL = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#0f7b8a"
Q = 0.5 * 1000 * 0.25**2 * 18.5e-4; CD = 0.03; U = 0.25

def load(n):
    d = pd.read_csv(os.path.join(HERE, f"sweep_{n}.csv"))
    if "error" in d: d = d[d.error.isna()].copy()
    d["T_v"] = (d.CT - CD) * Q * 1e3; d["eta_v"] = (d.CT - CD) / d.CP
    d = d[d.eta <= 0.9].copy()          # 无黏 η > 0.9 的点（多为 ψ=75°、强羽化）功率账不可信，剔除
    return d
A, B, C = load("A"), load("B"), load("C")
geo = json.load(open(os.path.join(HERE, "sweep_C_geo.json")))

# 最优点（C：最大推力；knee：η ≥ 0.6 里推力最大）
Cp = C[C.psi == 90]
best = C.loc[C.T_v.idxmax()]                                             # 约束盒角点（上界）
knee = C[C.St <= 0.80].sort_values("T_v", ascending=False).iloc[0]       # 可信最优：St ≤ 0.8（验证范围外推有限）
eff  = C[(C.T_v >= 40)].sort_values("eta_v", ascending=False).iloc[0]    # 高效点

fig = plt.figure(figsize=(18, 11.5), facecolor=SURF)
gs = fig.add_gridspec(2, 3, left=0.045, right=0.985, top=0.86, bottom=0.06, hspace=0.42, wspace=0.26)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, 0]); axE = fig.add_subplot(gs[1, 1]); axF = fig.add_subplot(gs[1, 2])
def style(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("bottom", "left"): ax.spines[s].set_color(GRID)
    ax.grid(color=GRID, lw=0.8); ax.tick_params(colors=INK2, labelsize=9)
for ax in (axB, axC, axD, axE, axF): style(ax)

# ---- A 平面形 + 截面
pl = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_planform_mm.csv"), delimiter=",", skiprows=1)
sec = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_section_NACA0021.csv"), delimiter=",", skiprows=1)
axA.set_facecolor(WATER); axA.set_aspect("equal")
for s in axA.spines.values(): s.set_visible(False)
axA.set_xticks([]); axA.set_yticks([])
axA.add_patch(Polygon(np.column_stack([pl[:, 1], -pl[:, 0]]), closed=True, color=BLUE, alpha=0.85, zorder=3))
b = pl[:, 1].max() - pl[:, 1].min(); c_mid = pl[:, 0].max() - pl[:, 0].min()
cm = 28.7
axA.plot([0], [-0.30 * cm], "o", color=GOLD, ms=10, mec=INK, mew=1.4, zorder=6)
axA.plot([-45, 45], [-0.30 * cm, -0.30 * cm], color=GOLD, lw=1.2, ls=(0, (3, 2)), zorder=5)
axA.text(0, 5, "来流 ↓", ha="center", fontsize=9, color="#4f86c6")
axA.annotate("", (-b / 2, -42), (b / 2, -42), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.1))
axA.text(0, -46, f"展 b = {b:.0f} mm", ha="center", va="top", fontsize=9, color=INK2)
axA.text(b / 2 + 4, -cm / 2, f"展中弦 {cm:.1f} mm\n翼尖弦 14.3 mm\n后掠 35°\nAR 4 · S 18.5 cm²", fontsize=9, color=INK2, va="center")
# 截面小图
sx = -46; sy = -68; sc = 60
axA.add_patch(Polygon(np.column_stack([sx + sec[:, 0] * sc, sy + sec[:, 1] * sc]), closed=True, color=INK, alpha=0.85, zorder=4))
axA.plot(sx + 0.30 * sc, sy, "o", color=GOLD, ms=7, mec=INK, mew=1.2, zorder=6)
axA.text(sx + sc + 4, sy, "截面 NACA 0021\nt_max 21% @ 0.30c（Kelly 2023 最优 = 海豚 0.218c @ 0.285c）\n前缘 R 1.4 mm · 尾缘闭合 · 俯仰轴 0.30c（金点）", fontsize=8.6, color=INK2, va="center")
axA.set_xlim(-60, 92); axA.set_ylim(-84, 12)
axA.set_title("A  最优脚部形态（俯视 + 截面）", loc="left", fontsize=12, color=INK)

# ---- B 推力地图 T(f, h0)（C 阶段，ψ90，各点取 α_max 与俯仰轴最优）
g = Cp.loc[Cp.groupby(["f", "h0"]).T_v.idxmax()]
fs, hs = sorted(g.f.unique()), sorted(g.h0.unique())
Z = np.array([[g[(g.f == f_) & (g.h0 == h_)].T_v.iloc[0] for h_ in hs] for f_ in fs])
im = axB.imshow(Z, origin="lower", aspect="auto", cmap="Blues", extent=[hs[0] * 1e3 - 5, hs[-1] * 1e3 + 5, fs[0] - 0.25, fs[-1] + 0.25])
for i, f_ in enumerate(fs):
    for j, h_ in enumerate(hs):
        r = g[(g.f == f_) & (g.h0 == h_)].iloc[0]
        axB.text(h_ * 1e3, f_, f"{r.T_v:.0f} mN\nη {r.eta_v:.2f}\nSt {r.St:.2f}", ha="center", va="center", fontsize=8.2,
                 color="white" if r.T_v > 0.55 * Z.max() else INK)
axB.set_xlabel("沉浮半幅 h0 [mm]", color=INK2); axB.set_ylabel("频率 f [Hz]", color=INK2)
axB.set_xticks([h * 1e3 for h in hs]); axB.set_yticks(fs); axB.grid(False)
axB.set_title("B  单桨推力地图（ψ = 90°，α_max、俯仰轴取最优）", loc="left", fontsize=12, color=INK)

# ---- C Pareto 推力–功率
sc_ = axC.scatter(C.P_mW, C.T_v, c=C.eta_v, cmap="viridis", s=18, vmin=0.3, vmax=0.85, alpha=0.75, edgecolors="none")
cb = plt.colorbar(sc_, ax=axC, pad=0.01); cb.set_label("效率 η（含 C_D 0.03）", color=INK2, fontsize=9); cb.ax.tick_params(labelsize=8, colors=INK2)
for r, col, lab in ((best, RED, "盒角上界"), (knee, GOLD, "可信最优 St≤0.8"), (eff, GREEN, "高效点")):
    axC.plot(r.P_mW, r.T_v, "o", color=col, ms=11, mec=INK, mew=1.3, zorder=6)
    axC.annotate(f"{lab}\n{r.T_v:.0f} mN · {r.P_mW:.0f} mW · η {r.eta_v:.2f}\nf {r.f:g} Hz · h0 {r.h0*1e3:.0f} · α {r.amax:g}° · 轴 {r['pivot']:g}c",
                 (r.P_mW, r.T_v), {"盒角上界": (r.P_mW - 62, r.T_v + 8), "可信最优 St≤0.8": (r.P_mW - 70, r.T_v - 62), "高效点": (r.P_mW + 6, r.T_v + 45)}[lab], fontsize=8.4, color=col,
                 arrowprops=dict(arrowstyle="-", color=col, lw=0.8))
axC.axhspan(21, 36, color="#fde8ee", zorder=0); axC.text(2, 30, "整机阻力（腿整流后）21–36 mN ÷ 4 桨 → 每桨只需 5–9 mN", fontsize=8, color=RED, va="center")
axC.set_xlabel("输入功率 P [mW]（沉浮 + 俯仰做功）", color=INK2); axC.set_ylabel("周期平均推力 T [mN]", color=INK2)
axC.set_xlim(0, C.P_mW.max() * 1.05); axC.set_ylim(-20, C.T_v.max() * 1.12)
axC.set_title("C  推力 vs 功率（最优几何，324 组运动学）", loc="left", fontsize=12, color=INK)

# ---- D 几何效应（B 阶段，maxT 运动学）
Bm = B[B.tag == "maxT"]
for fam, col in (("rect", MUTED), ("delta", GOLD), ("fishtail", GREEN), ("ellipse", PUR), ("lunate", RED)):
    x = Bm[Bm.family == fam].groupby("AR").T_v.max()
    axD.plot(x.index, x.values, "-o", color=col, lw=2, ms=5, label={"rect": "矩形", "delta": "前缘直/尾缘前掠", "fishtail": "尾缘直/前缘后掠", "ellipse": "椭圆", "lunate": "月牙（前后缘后掠 35°）"}[fam])
axD.axvline(4, color=RED, lw=1, ls=(0, (4, 3))); axD.text(2.05, Bm.T_v.max() - 8, "选 AR 4：Lagopoulos 2023（黏性，Re 8500）\nC_T 到 AR 4 饱和；再高弦太短、Re 掉\n（UVLM 无黏所以这里不饱和）", fontsize=8, color=RED, va="top")
axD.set_xlabel("展弦比 AR（面积固定 18.5 cm²）", color=INK2); axD.set_ylabel("推力 T [mN]（f 2 Hz, h0 50, α 30°）", color=INK2)
axD.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower right"); axD.set_xticks([2, 3, 4, 5, 6])
axD.set_title("D  平面形与展弦比（UVLM，无黏 + C_D 0.03）", loc="left", fontsize=12, color=INK)

# ---- E 最优轨迹频闪 + α 时程
from cfd_matrix_kin import Case, foil_poly
kb = Case("best", "上界（盒角）", St=float(best.St), al_max=float(best.amax), psi=float(best.psi), pivot=float(best["pivot"]), c=cm * 1e-3, h0c=float(best.h0) / (cm * 1e-3), u=U, color=RED)
kk = Case("knee", "可信最优", St=float(knee.St), al_max=float(knee.amax), psi=float(knee.psi), pivot=float(knee["pivot"]), c=cm * 1e-3, h0c=float(knee.h0) / (cm * 1e-3), u=U, color=GOLD)
axE.set_facecolor(WATER); axE.set_aspect("equal"); axE.grid(False)
for s in axE.spines.values(): s.set_visible(False)
axE.set_xticks([]); axE.set_yticks([])
for k, (cs, y0) in enumerate(((kb, 0.06), (kk, -0.06))):
    dx = 0.19
    taus = np.linspace(0, 1, 17, endpoint=False)
    axE.plot(np.linspace(0, dx, 200), y0 + np.array([cs.at(t)["h"] for t in np.linspace(0, 1, 200)]), color=cs.color, lw=0.9, alpha=0.3)
    for tau, x in zip(taus, np.linspace(0, dx, 17, endpoint=False)):
        s = cs.at(tau); P = np.array([x, y0 + s["h"]]); poly, LE, TE = foil_poly(P, s["th"], cs.c, cs.pivot, thk=0.006)
        a_ = 0.3 + 0.7 * (0.5 + 0.5 * np.cos(2 * np.pi * tau))**0.8
        axE.add_patch(Polygon(poly, closed=True, color=cs.color, alpha=a_, lw=0)); axE.plot(*P, "o", color=GOLD, ms=2.5, mec="none", zorder=5)
    axE.text(-0.005, y0, f"{cs.label}\nf {cs.f:.2f} Hz · h0 ±{cs.h0*1e3:.0f} mm\nθ0 ±{np.degrees(cs.th0):.0f}° · ψ {cs.psi_d:.0f}°\nα_max {np.degrees(np.abs(cs.al).max()):.0f}° · St {cs.St:.2f}\n俯仰峰速 {cs.rate_pk:.0f} °/s", ha="right", va="center", fontsize=8.6, color=cs.color)
axE.set_xlim(-0.10, 0.20); axE.set_ylim(-0.135, 0.135)
axE.set_title("E  最佳轨迹（一个周期频闪，来流从右向左）", loc="left", fontsize=12, color=INK)

# ---- F ψ / α_max / 俯仰轴 的影响（C 阶段，best 的 f,h0）
sub = C[(np.isclose(C.f, knee.f)) & (np.isclose(C.h0, knee.h0))]
for pv, col in ((0.2, BLUE), (0.3, GREEN), (0.4, PUR)):
    for ps, ls in ((75, ":"), (90, "-"), (105, "--")):
        x = sub[(np.isclose(sub["pivot"], pv)) & (sub.psi == ps)].sort_values("amax")
        if len(x): axF.plot(x.amax.values, x.T_v.values, ls=ls, color=col, lw=1.8, marker="o", ms=4, label=f"轴 {pv:g}c · ψ {ps}°")
axF.set_xlabel("α_max [°]", color=INK2); axF.set_ylabel(f"推力 T [mN]（f {knee.f:g} Hz, h0 {knee.h0*1e3:.0f} mm, St {knee.St:.2f}）", color=INK2)
axF.legend(frameon=False, fontsize=7.8, labelcolor=INK2, ncol=2); axF.set_xticks([20, 25, 30])
axF.set_title("F  攻角上限与俯仰轴（ψ = 90°；此 St 下 ψ 75/105° 压不住攻角上限）", loc="left", fontsize=12, color=INK)

fig.text(0.045, 0.965, "BODY2 · 海豚尾鳍式脚蹼：最优形态 + 最佳轨迹（U = 0.25 m/s 巡航，S = 18.5 cm²，f ≤ 3 Hz，h0 ≤ 50 mm，α_max ≤ 30°）", fontsize=16, color=INK)
fig.text(0.045, 0.925,
    f"方法：三维非定常涡格法 UVLM（PteraSoftware 5.1，对 Schouveiler 2005 验证 C_T 0.354 vs 0.32±0.01、η 0.75 vs 0.73±0.04）+ 黏性偏置 C_D 0.03。"
    f"无黏方法不含失速：α_max 30° 的推力上偏约 20%（Read 2003 对照），St > 0.8 的点当上界看。\n"
    f"几何最优：月牙形 AR 4、后掠 35°、锥度 0.5、NACA 0021 截面。轨迹最优：ψ = 90°、α_max 顶到上限、俯仰轴 0.3c；推力随 St 单调升到 ≈1.0 后被羽化幅值卡住。",
    fontsize=9.6, color=MUTED, va="top", linespacing=1.7)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"fig_fluke_optimum.{ext}"), dpi=160 if ext == "png" else None, facecolor=SURF)
json.dump(dict(best=best.to_dict(), knee=knee.to_dict(), eff=eff.to_dict(), geo=geo), open(os.path.join(HERE, "fluke_optimum.json"), "w"), indent=1, ensure_ascii=False, default=float)
print("best", best[["f", "h0", "amax", "psi", "pivot", "St", "theta0", "T_v", "P_mW", "eta_v", "rate_pk"]].round(3).to_dict())
print("knee", knee[["f", "h0", "amax", "psi", "pivot", "St", "theta0", "T_v", "P_mW", "eta_v", "rate_pk"]].round(3).to_dict())
print("eff ", eff[["f", "h0", "amax", "psi", "pivot", "St", "theta0", "T_v", "P_mW", "eta_v", "rate_pk"]].round(3).to_dict())
