#!/usr/bin/env python3
"""CFD 算例矩阵演示视频：七个场景依次播放，每个场景把该组的变体并排实时拍动。
只画运动学（沉浮 + 俯仰 + 尾迹轨迹），不含任何受力估计。
→ cfd_matrix_demo.mp4
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, animation
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch, Polygon
from cfd_matrix_kin import SCENES, foil_poly, U, SERVO_RATE

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
GOLD, RED = "#d4a017", "#c2185b"

FPS, SEC = 30, 5.0
NF_SCENE = int(SEC*FPS)
DT = 1/FPS
SLOTW, GAP, YH = 0.150, 0.052, 0.062      # 每个变体的槽宽、间距、半高
WMAX = 4*SLOTW + 3*GAP
TRAIL_S = 0.62                            # 尾迹长度（相对槽宽）

fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gsp = fig.add_gridspec(2, 1, height_ratios=[1.75, 1], left=0.048, right=0.975,
                       top=0.775, bottom=0.085, hspace=0.24)
axM = fig.add_subplot(gsp[0]); axA = fig.add_subplot(gsp[1])

axM.set_xticks([]); axM.set_yticks([])
for s in axM.spines.values(): s.set_visible(False)
axM.set_xlim(-0.016, WMAX + 0.016); axM.set_ylim(-YH*1.72, YH*1.72); axM.set_aspect("equal")
for s in ("top", "right"): axA.spines[s].set_visible(False)
for s in ("bottom", "left"): axA.spines[s].set_color(GRID)
axA.grid(color=GRID, lw=0.8); axA.tick_params(colors=INK2, labelsize=9)
axA.set_xlabel("时间 t [s]（各算例按各自频率实时拍动）", color=INK2, fontsize=9.5)
axA.set_ylabel("攻角 α [°]", color=INK2, fontsize=9.5)

TIT = fig.text(0.048, 0.955, "", fontsize=19, color=INK, ha="left", va="top")
SUB = fig.text(0.048, 0.905, "", fontsize=11.2, color=MUTED, ha="left", va="top", linespacing=1.7)
FOOT = fig.text(0.048, 0.022, "BODY2 · 海豚尾鳍式拍动 —— CFD 算例矩阵演示。"
                "基准：弦 c = 42 mm（120° 展开折扇，AR = 2.87）、h0/c = 0.83、U = 0.25 m/s。"
                "只画运动学，不含受力。",
                fontsize=9.2, color=MUTED, ha="left", va="bottom")
PROG = fig.text(0.975, 0.022, "", fontsize=9.6, color=INK2, ha="right", va="bottom")

state = {}


def build_scene(si):
    sc = SCENES[si]
    cases = sc["cases"]; n = len(cases)
    for a in state.get("art", []):
        try: a.remove()
        except Exception: pass
    axA.cla()
    for s in ("top", "right"): axA.spines[s].set_visible(False)
    for s in ("bottom", "left"): axA.spines[s].set_color(GRID)
    axA.grid(color=GRID, lw=0.8); axA.tick_params(colors=INK2, labelsize=9)
    axA.set_xlabel("时间 t [s]（各算例按各自频率实时拍动）", color=INK2, fontsize=9.5)
    axA.set_ylabel("攻角 α [°]", color=INK2, fontsize=9.5)
    axA.set_xlim(0, SEC); axA.axhline(0, color=INK2, lw=0.9)
    amx = max(np.degrees(np.abs(c.al).max()) for c in cases)
    axA.set_ylim(-amx*1.45, amx*1.45)

    art = []
    span = n*SLOTW + (n - 1)*GAP
    xoff = (WMAX - span)/2
    box = FancyBboxPatch((xoff - 0.014, -YH*1.04), span + 0.028, YH*2.08,
                         boxstyle="round,pad=0.005", facecolor=WATER, edgecolor="none", zorder=0)
    axM.add_patch(box); art.append(box)

    foils, pivs, trails, labs, aline, amark, x0s = [], [], [], [], [], [], []
    for k, cs in enumerate(cases):
        x0 = xoff + k*(SLOTW + GAP) + SLOTW*0.72   # 脚蹼所在的 x（尾迹向左拖）
        x0s.append(x0)
        p = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=cs.color, zorder=8)
        axM.add_patch(p); foils.append(p)
        pv, = axM.plot([], [], "o", color=GOLD, ms=7, mec=INK, mew=1.1, zorder=10)
        pivs.append(pv)
        tr, = axM.plot([], [], color=cs.color, lw=1.4, alpha=0.40, zorder=2)
        trails.append(tr)
        marg = SERVO_RATE/cs.rate_pk
        t1 = axM.text(x0 - SLOTW*0.22, YH*1.16, cs.label, ha="center", va="bottom",
                      fontsize=12.5, color=cs.color, weight="bold")
        t2 = axM.text(x0 - SLOTW*0.22, -YH*1.10,
                      f"f {cs.f:.2f} Hz · χ {cs.chi:.2f} · θ0 {np.degrees(cs.th0):.0f}°\n"
                      f"俯仰峰速 {cs.rate_pk:.0f} °/s · 舵机余量 {marg:.1f}×"
                      + ("　← 超出舵机" if marg < 1.0 else ""),
                      ha="center", va="top", fontsize=9.2,
                      color=RED if marg < 1.5 else INK2,
                      weight="bold" if marg < 1.5 else "normal")
        t3 = axM.text(x0 + SLOTW*0.30, YH*0.86, "", ha="right", va="center",
                      fontsize=10.5, color=INK)
        labs += [t1, t2, t3]
        ln, = axA.plot([], [], color=cs.color, lw=2.1, label=cs.label); aline.append(ln)
        mk, = axA.plot([], [], "o", color=cs.color, ms=6.5, mec="white", mew=1.2); amark.append(mk)
        art += [p, pv, tr, ln, mk]
    art += labs
    fa = FancyArrowPatch((0.018, -YH*1.45), (-0.004, -YH*1.45), arrowstyle="-|>",
                         mutation_scale=11, color="#9fc3e6", lw=1.3, zorder=3)
    axM.add_patch(fa)
    ft = axM.text(0.024, -YH*1.45, f"来流 U = {U} m/s（机体向右前进）", fontsize=9.5,
                  color="#4f86c6", va="center")
    art += [fa, ft]
    if n > 1:
        lg = axA.legend(frameon=False, fontsize=9.5, loc="upper right", ncol=n, labelcolor=INK2)
        art.append(lg)

    TIT.set_text(sc["title"])
    SUB.set_text(sc["why"])
    state.update(art=art, cases=cases, foils=foils, pivs=pivs, trails=trails, labs=labs,
                 aline=aline, amark=amark, x0s=x0s, si=si,
                 tx=[[] for _ in cases], ty=[[] for _ in cases],
                 at=[[] for _ in cases], av=[[] for _ in cases], circ=[])


def draw(fi):
    si, k = divmod(fi, NF_SCENE)
    if k == 0:
        build_scene(si)
    t = k*DT
    cases = state["cases"]
    for j, cs in enumerate(cases):
        tau = (t/cs.T) % 1.0
        s = cs.at(tau)
        x0 = state["x0s"][j]
        P = np.array([x0, s["h"]])
        poly, LE, TE = foil_poly(P, s["th"], cs.c, cs.pivot)
        state["foils"][j].set_xy(poly)
        state["pivs"][j].set_data([P[0]], [P[1]])
        # 尾迹：历史点以 U 向左漂移
        tx, ty = state["tx"][j], state["ty"][j]
        for i in range(len(tx)): tx[i] -= U*DT
        tx.append(x0); ty.append(s["h"])
        while tx and tx[0] < x0 - SLOTW*TRAIL_S:
            tx.pop(0); ty.pop(0)
        state["trails"][j].set_data(tx, ty)
        state["labs"][3*j + 2].set_text(f"α = {np.degrees(s['al']):+5.1f}°")
        at, av = state["at"][j], state["av"][j]
        at.append(t); av.append(np.degrees(s["al"]))
        state["aline"][j].set_data(at, av)
        state["amark"][j].set_data([t], [av[-1]])
    PROG.set_text(f"场景 {si+1} / {len(SCENES)}")
    return []


NF = NF_SCENE*len(SCENES)
anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000/FPS, blit=False)
out = os.path.join(HERE, "cfd_matrix_demo.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=5200,
          extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out, f"（{NF} 帧，{NF/FPS:.0f} s）")
