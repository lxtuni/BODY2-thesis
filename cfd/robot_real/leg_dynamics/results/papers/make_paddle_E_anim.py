#!/usr/bin/env python3
"""14 号动画：颈部扭簧铰铲桨一个周期（侧视 + 舵机/俯仰/攻角曲线）→ paddle_E_anim.mp4（放慢 3×）"""
import os, json, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation, font_manager
from matplotlib.patches import Polygon, Circle, Wedge
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
BLUE, RED, GOLD, GREEN, PUR, HULL, STEEL = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#c9c4b8", "#6b7280"
D = json.load(open(os.path.join(HERE, "paddle_E_design.json")))
O2 = np.array([63.4, 4.0]); WL = -10.0
L_RC, L_PAR, PSI, H0, F, DPHI, DROP, L_S, X_P, STOP = D["L_RC"], D["L_PAR"], D["PSI_FIX"], D["H0"], D["F"], D["DPHI"], D["DROP"], D["L_S"], D["X_P"], D["STOP"]
TH0, PSI_P = D["best"]["th0"], D["best"]["psi"]; U = 0.25
HEAD_L, THK = 42.0, 5.6
def rot(a): a = math.radians(a); return np.array([math.cos(a), math.sin(a)])
def fk(phi): C = O2 + L_RC * rot(PSI); return C, C + L_PAR * rot(PSI + phi)

SLOW, FPS, NCYC = 3, 30, 2
NF = int(NCYC / F * SLOW * FPS); tt = np.arange(NF) / FPS / SLOW
w = 2 * math.pi * F
phi_t = 90 + DPHI * np.sin(w * tt); z_t = H0 * np.sin(w * tt)
th_t = TH0 * np.sin(w * tt + math.radians(PSI_P))
al_t = np.degrees(np.arctan2(H0 * 1e-3 * w * np.cos(w * tt), U)) - th_t

fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], left=0.03, right=0.98, top=0.9, bottom=0.07, hspace=0.35, wspace=0.12)
axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[1, 1])
fig.text(0.03, 0.945, f"14 · 颈部扭簧铰铲桨：ψ 固定 −90°，φ = 90° ± {DPHI:.0f}°（{F} Hz，放慢 {SLOW}×），俯仰被动（k {D['K']*1e3:.0f} mN·m/rad）", fontsize=14, weight="bold", color=INK)
fig.text(0.03, 0.915, f"被动响应（UVLM 流固耦合）：θ0 ±{TH0:.0f}°，领先沉浮 {PSI_P:.0f}°，α_max {D['best']['alpha_max']:.0f}°，推力 {D['best']['T']*1e3:.0f} mN/桨，η {D['best']['eta']:.2f}。连杆尺寸为照片估计值。", fontsize=10, color=INK2)

def draw(i):
    axA.cla(); axA.set_facecolor(SURF); axA.set_aspect("equal"); axA.set_xticks([]); axA.set_yticks([])
    for s in axA.spines.values(): s.set_visible(False)
    axA.axhspan(-150, WL, color=WATER, zorder=0); axA.axhline(WL, color=BLUE, lw=1, ls="--"); axA.text(190, WL + 3, "水线", color=BLUE, fontsize=9)
    axA.add_patch(Polygon([(-60, 70), (200, 70), (215, 40), (200, 0), (170, WL - 6), (-40, WL - 6), (-60, 20)], closed=True, fc=HULL, ec=STEEL, lw=0.8, alpha=0.6, zorder=1))
    phi, th = phi_t[i], th_t[i]
    C, E = fk(phi); Dp = O2 + L_PAR * rot(PSI + phi)
    for a, b in ((O2, C), (O2, Dp), (C, E), (Dp, E)):
        axA.plot([a[0], b[0]], [a[1], b[1]], color=STEEL, lw=3, solid_capstyle="round", zorder=3)
    K = E + np.array([0, -DROP]); P = K + np.array([-L_S, 0])
    axA.plot([E[0], K[0], P[0]], [E[1], K[1], P[1]], color=GOLD, lw=4, solid_capstyle="round", zorder=4)
    u = -rot(th); n = np.array([-u[1], u[0]])
    poly = [P + X_P * u - n * THK / 2, P + (X_P + HEAD_L) * u - n * THK / 2, P + (X_P + HEAD_L) * u + n * THK / 2, P + X_P * u + n * THK / 2]
    axA.add_patch(Polygon(poly, closed=True, fc=RED, ec=RED, zorder=5))
    axA.add_patch(Wedge(P, HEAD_L + 6, 180 - STOP, 180 + STOP, fc=GREEN, alpha=0.1, ec=GREEN, lw=0.8, zorder=2))
    axA.add_patch(Circle(P, 2.4, fc="white", ec=INK, lw=1, zorder=6)); axA.add_patch(Circle(E, 3.4, fc=BLUE, ec="white", lw=1.2, zorder=7))
    axA.plot(*O2, "o", color=INK, ms=6, zorder=7)
    # 相对来流矢量 & 攻角
    vz = H0 * 1e-3 * w * math.cos(w * tt[i])
    M = P + (X_P + HEAD_L * 0.5) * u
    vv = np.array([-U, -vz]) / math.hypot(U, vz) * 40
    axA.annotate("", xy=M, xytext=M - vv, arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.8), zorder=8)
    axA.text(M[0] - 55, M[1] - 30, f"相对来流\nα = {al_t[i]:+.0f}°", fontsize=9, color=BLUE, ha="center")
    axA.text(E[0] + 8, E[1] + 6, "E", color=BLUE, fontsize=10, weight="bold"); axA.text(P[0] + 4, P[1] - 12, "颈部铰 P", fontsize=9, color=INK)
    axA.text(O2[0] + 6, O2[1] + 6, f"髋 ψ = {PSI:.0f}°（固定）", fontsize=9, color=INK)
    axA.text(Dp[0] - 4, Dp[1] + 8, f"φ = {phi:.0f}°", fontsize=9, color=INK2)
    axA.text(P[0] - 28, P[1] + 22, f"θ = {th:+.0f}°", fontsize=10, color=RED, weight="bold")
    axA.text(-60, -140, f"t = {tt[i]*F:.2f} T", fontsize=10, color=INK2)
    axA.set_xlim(-80, 230); axA.set_ylim(-150, 80)
    for ax, ys, labs, cols, ttl in ((axB, (phi_t - 90, z_t), ("φ − 90° [°]", "E 沉浮 z [mm]"), (STEEL, BLUE), "舵机指令与沉浮"),
                                    (axC, (th_t, al_t), ("被动俯仰 θ [°]", "头部攻角 α [°]"), (RED, GREEN), "俯仰响应与攻角")):
        ax.cla(); ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.6)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        for y, lab, col in zip(ys, labs, cols):
            ax.plot(tt * F, y, color=col, lw=1.8, label=lab); ax.plot(tt[i] * F, y[i], "o", color=col, ms=7)
        ax.axvline(tt[i] * F, color=INK2, lw=0.8); ax.axhline(0, color=INK, lw=0.5)
        ax.set_xlim(0, NCYC); ax.set_ylim(-45, 45); ax.legend(fontsize=9, frameon=False, loc="upper right", ncol=2); ax.set_title(ttl, fontsize=10, loc="left", color=INK)
        ax.tick_params(labelsize=8.5, colors=INK2)
    axC.set_xlabel("t / T")
    return []

anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000 / FPS, blit=False)
out = os.path.join(HERE, "paddle_E_anim.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=5000, extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
draw(int(NF * 0.12)); fig.savefig(os.path.join(HERE, "fig_paddle_E_frame.png"), dpi=130, facecolor=SURF)
print("→", out, NF, "frames")
