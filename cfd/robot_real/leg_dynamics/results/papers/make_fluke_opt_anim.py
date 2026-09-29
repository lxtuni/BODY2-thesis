#!/usr/bin/env python3
"""最佳尾鳍轨迹演示：上界（盒角）vs 可信最优（St 0.8），并排实时拍动 + 俯视平面形 → fluke_optimum_anim.mp4"""
import os, json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, animation
from matplotlib.patches import Polygon, FancyArrowPatch, FancyBboxPatch
from cfd_matrix_kin import Case, foil_poly
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
RED, GOLD, BLUE = "#c2185b", "#d4a017", "#2a78d6"
U = 0.25; CM = 0.0287
opt = json.load(open(os.path.join(HERE, "fluke_optimum.json")))
pl = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_planform_mm.csv"), delimiter=",", skiprows=1)
sec = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_section_NACA0021.csv"), delimiter=",", skiprows=1)
def mk(k, lab, col):
    r = opt[k]; return Case(k, lab, St=r["St"], al_max=r["amax"], psi=r["psi"], pivot=r["pivot"], c=CM, h0c=r["h0"]/CM, u=U, color=col), r
cases = [mk("best", "上界（约束盒角）", RED), mk("knee", "可信最优（St ≤ 0.8）", GOLD)]

FPS, SEC, SLOW = 30, 10.0, 2.0     # 放慢 2×
NF = int(SEC*FPS)
fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.75], height_ratios=[2.2, 1], left=0.04, right=0.98, top=0.80, bottom=0.08, hspace=0.35, wspace=0.18)
axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]; axP = fig.add_subplot(gs[0, 2]); axA = fig.add_subplot(gs[1, :])
for ax in axs + [axP]:
    ax.set_facecolor(WATER); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
for s in ("top", "right"): axA.spines[s].set_visible(False)
axA.grid(color=GRID, lw=0.8); axA.tick_params(colors=INK2, labelsize=9); axA.axhline(0, color=INK2, lw=0.9)
axA.set_xlim(0, SEC/SLOW); axA.set_ylim(-36, 36); axA.set_xlabel("真实时间 t [s]（视频放慢 2×）", color=INK2); axA.set_ylabel("攻角 α [°]", color=INK2)
# 俯视平面形
axP.add_patch(Polygon(np.column_stack([pl[:, 1], -pl[:, 0]]), closed=True, color=BLUE, alpha=0.85))
axP.plot([0], [-0.3*28.7], "o", color=GOLD, ms=9, mec=INK, mew=1.2); axP.text(0, 6, "来流 ↓", ha="center", fontsize=9, color="#4f86c6")
axP.add_patch(Polygon(np.column_stack([-40 + sec[:, 0]*50, -62 + sec[:, 1]*50]), closed=True, color=INK, alpha=0.85))
axP.text(14, -62, "NACA 0021", fontsize=9, color=INK2, va="center")
axP.set_xlim(-50, 50); axP.set_ylim(-75, 12)
axP.set_title("最优形态：月牙 AR 4 · 后掠 35° · S 18.5 cm²", loc="left", fontsize=10.5, color=INK)
fig.text(0.04, 0.955, "BODY2 · 最佳尾鳍轨迹（U = 0.25 m/s 巡航，UVLM 优化结果）", fontsize=18, color=INK, va="top")
fig.text(0.04, 0.905, "两条轨迹：左 = 约束盒角（f 3 Hz、h0 50 mm、St 1.2，无黏上界，α 30°）；右 = 可信最优（f 2 Hz、h0 50 mm、St 0.8）。"
         "共同点：ψ = 90°（俯仰领先 1/4 周期）、攻角顶到 30° 上限、俯仰轴 0.3–0.4c。\n推力：346 / 260 mN 每桨（含 C_D 0.03；无失速模型，上偏 ~20%）；四桨远超整机阻力 21–36 mN。", fontsize=10, color=MUTED, va="top", linespacing=1.6)
arts = []
for ax, (cs, r) in zip(axs, cases):
    ax.set_xlim(-0.16, 0.06); ax.set_ylim(-0.11, 0.11)
    ax.set_title(f"{cs.label}：f {cs.f:.1f} Hz · h0 ±{cs.h0*1e3:.0f} mm · θ0 ±{np.degrees(cs.th0):.0f}° · St {cs.St:.2f}", loc="left", fontsize=10.5, color=cs.color)
    ax.text(-0.155, -0.10, f"推力 {r['T_v']:.0f} mN · 功率 {r['P_mW']:.0f} mW · η {r['eta_v']:.2f} · 俯仰峰速 {r['rate_pk']:.0f} °/s", fontsize=9, color=INK2)
    ax.add_patch(FancyArrowPatch((-0.10, 0.095), (-0.13, 0.095), arrowstyle="-|>", mutation_scale=11, color="#9fc3e6", lw=1.3)); ax.text(-0.095, 0.095, "来流 U", fontsize=9, color="#4f86c6", va="center")
    foil = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=cs.color, zorder=8); ax.add_patch(foil)
    piv, = ax.plot([], [], "o", color=GOLD, ms=8, mec=INK, mew=1.2, zorder=10)
    tr, = ax.plot([], [], color=cs.color, lw=1.4, alpha=0.4)
    txt = ax.text(0.03, 0.08, "", fontsize=10, color=INK, ha="right")
    ln, = axA.plot([], [], color=cs.color, lw=2, label=cs.label); mk_, = axA.plot([], [], "o", color=cs.color, ms=6, mec="white")
    arts.append(dict(cs=cs, ax=ax, foil=foil, piv=piv, tr=tr, txt=txt, ln=ln, mk=mk_, tx=[], ty=[], at=[], av=[]))
axA.legend(frameon=False, fontsize=9, ncol=2, labelcolor=INK2, loc="upper right")
DT = 1/FPS/SLOW
def draw(fi):
    t = fi*DT
    for a in arts:
        cs = a["cs"]; s = cs.at((t/cs.T) % 1.0); P = np.array([0.0, s["h"]])
        poly, LE, TE = foil_poly(P, s["th"], cs.c, cs.pivot, thk=0.006); a["foil"].set_xy(poly); a["piv"].set_data([0], [s["h"]])
        for i in range(len(a["tx"])): a["tx"][i] -= U*DT
        a["tx"].append(0.0); a["ty"].append(s["h"])
        while a["tx"] and a["tx"][0] < -0.15: a["tx"].pop(0); a["ty"].pop(0)
        a["tr"].set_data(a["tx"], a["ty"]); a["txt"].set_text(f"α = {np.degrees(s['al']):+5.1f}°  θ = {np.degrees(s['th']):+5.1f}°")
        a["at"].append(t); a["av"].append(np.degrees(s["al"])); a["ln"].set_data(a["at"], a["av"]); a["mk"].set_data([t], [a["av"][-1]])
    return []
anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000/FPS, blit=False)
out = os.path.join(HERE, "fluke_optimum_anim.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=5000, extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out)
