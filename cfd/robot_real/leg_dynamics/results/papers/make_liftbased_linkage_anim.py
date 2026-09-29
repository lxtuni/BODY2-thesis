#!/usr/bin/env python3
"""升力型运动 —— 用 BODY2 原来的双连杆实现（摇臂 + 平行四边形 + 脚踝俯仰）。
→ liftbased_linkage.mp4
"""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, animation
from matplotlib.patches import FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import liftbased_linkage as K

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
RED, BLUE, GOLD, GREEN, WATER, PUR = "#c2185b", "#2a78d6", "#d4a017", "#16a085", "#dbe9f6", "#7b3f9d"

T, F, U = K.T, K.F, K.U
tt = np.linspace(0, T, 721)
S  = [K.state(t) for t in tt]
Pz = np.array([s["P"][1] for s in S]); Px = np.array([s["P"][0] for s in S])
TH = np.degrees([s["th"] for s in S]); AL = np.degrees([s["al"] for s in S])
PS = np.degrees([s["psi"] for s in S]); FX = np.array([s["Fx"] for s in S])
ARM_Q = np.array([s["xcp"] - 0.25 for s in S])
ARM_M = np.array([s["xcp"] - 0.50 for s in S])
THR   = np.degrees(np.gradient(np.radians(TH), tt))
HEAVE = (Pz.max()-Pz.min()); SURGE = (Px.max()-Px.min()); ST = F*HEAVE/U
TH_MID = 0.5*(TH.max()+TH.min())

FPS, SLOW, CYC = 30, 3, 2
NF = int(CYC*T*FPS*SLOW)
fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(3, 2, width_ratios=[1.30, 1], height_ratios=[1, 1, 1],
                      left=0.045, right=0.972, top=0.795, bottom=0.065, wspace=0.17, hspace=0.55)
axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 1]); axD = fig.add_subplot(gs[2, 1])

def style(ax, water=True):
    ax.set_facecolor(WATER if water else SURF)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    if water:
        for s_ in ("bottom", "left"): ax.spines[s_].set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.grid(color=GRID, lw=0.8); ax.tick_params(colors=INK2, labelsize=8.5)
        for s_ in ("bottom", "left"): ax.spines[s_].set_color(GRID)
style(axA); [style(a, False) for a in (axB, axC, axD)]
axA.set_xlim(-0.062, 0.168); axA.set_ylim(-0.098, 0.078); axA.set_aspect("equal")
SC = 0.042/max(np.hypot(FX, [s["Fz"] for s in S]))

def arrow(ax, p0, v, c, lw=2.4, ms=15, z=9, ls="-"):
    a = FancyArrowPatch(p0, (p0[0]+v[0], p0[1]+v[1]), arrowstyle="-|>", mutation_scale=ms,
                        color=c, lw=lw, zorder=z, linestyle=ls, shrinkA=0, shrinkB=0)
    ax.add_patch(a); return a

for xx in np.linspace(-0.055, 0.145, 6):
    axA.add_patch(FancyArrowPatch((xx+0.018, 0.068), (xx, 0.068), arrowstyle="-|>",
                                  mutation_scale=9, color="#9fc3e6", lw=1.1, zorder=1))
axA.text(0.050, 0.073, f"来流（机体以 U = {U} m/s 向右前进）", ha="center", fontsize=9.5, color="#4f86c6")
axA.plot(Px, Pz, color=MUTED, lw=1.1, ls=(0, (3, 2)), alpha=0.8, zorder=2)

axA.set_title("A  原双连杆做升力型：摇臂摆动 → 桨板上下拍；脚踝只管把攻角控住",
              loc="left", fontsize=12.5, color=INK, pad=9)
axB.set_title("B  三个角度：髋摆 ψ、脚踝俯仰 θ、攻角 α", loc="left", fontsize=10.5, color=INK, pad=6)
axC.set_title("C  瞬时推力（准定常，非 CFD）", loc="left", fontsize=10.5, color=INK, pad=6)
axD.set_title("D  为什么枢轴不能放在中间：压心相对枢轴的位置", loc="left", fontsize=10.5, color=INK, pad=6)

axB.plot(tt/T, PS-PS.mean(), color=GREEN, lw=2.0, label=f"髋摆 ψ（±{np.degrees(K.PSI_A):.0f}°）")
axB.plot(tt/T, TH-TH_MID,    color=PUR,   lw=2.0, label=f"脚踝 θ（±{(TH.max()-TH.min())/2:.0f}°，已去零位偏置）")
axB.plot(tt/T, AL,           color=RED,   lw=2.0, label=f"攻角 α（±{np.abs(AL).max():.0f}°）")
axB.axhline(0, color=INK2, lw=0.9); axB.set_xlim(0, 1)
axB.set_xlabel("周期相位 τ", color=INK2, fontsize=9); axB.set_ylabel("角度 [°]", color=INK2, fontsize=9)
axB.legend(frameon=False, fontsize=8, loc="upper right", labelcolor=INK2, ncol=1)

axC.plot(tt/T, FX*1e3, color=BLUE, lw=2.1)
axC.axhline(0, color=INK2, lw=0.9)
axC.axhline(FX.mean()*1e3, color=BLUE, lw=1, ls=(0, (4, 3)))
axC.set_xlim(0, 1); axC.set_xlabel("周期相位 τ", color=INK2, fontsize=9)
axC.set_ylabel("推力 [mN]", color=INK2, fontsize=9)
axC.text(0.02, FX.mean()*1e3, f"  平均 {FX.mean()*1e3:+.1f} mN", color=BLUE, fontsize=8.5, va="bottom")
axC.text(0.98, axC.get_ylim()[1]*0.9, "两个半程都为正 —— 没有回程损失", fontsize=8.5,
         color=INK2, ha="right", va="top")

axD.axhline(0, color=INK, lw=1.2)
axD.plot(tt/T, ARM_Q*K.HEAD*1e3, color=GREEN, lw=2.2, label="枢轴在 1/4 弦：压心恒在枢轴之后 → 静稳定")
axD.plot(tt/T, ARM_M*K.HEAD*1e3, color=RED,   lw=2.2, label="枢轴在中间：压心恒在枢轴之前 → 静不稳定")
axD.set_xlim(0, 1); axD.set_xlabel("周期相位 τ", color=INK2, fontsize=9)
axD.set_ylabel("压心 − 枢轴 [mm]", color=INK2, fontsize=9)
axD.legend(frameon=False, fontsize=8, loc="center right", labelcolor=INK2)

fig.text(0.045, 0.958, "升力型运动，用 BODY2 原来的双连杆实现 —— 枢轴在 1/4 弦，髋有摆角",
         fontsize=15, color=INK, ha="left")
fig.text(0.045, 0.912,
    f"机构沿用原尺寸：摇臂 O2→C {K.L_ROCK_C*1e3:.1f} mm，从动杆 {K.L_PAR*1e3:.0f} mm，桨杆 {K.STEM*1e3:.0f} mm，"
    f"桨头 {K.HEAD*1e3:.0f} × {K.SPAN*1e3:.0f} mm（弦 × 展）　|　髋摆 ψ = {np.degrees(K.PSI_M):.0f}° ± {np.degrees(K.PSI_A):.0f}°\n"
    f"桨板枢轴走出的弧：垂直 {HEAVE*1e3:.0f} mm、前后仅 {SURGE*1e3:.0f} mm　→　St = {ST:.2f}（高效区 0.2–0.4 ✓）　|　"
    f"脚踝行程 ±{(TH.max()-TH.min())/2:.0f}°、峰值 {np.abs(THR).max():.0f} °/s（7465W 有 857 °/s，够用）\n"
    f"注意：这个构型要求摇臂工作在 {np.degrees(K.PSI_M-K.PSI_A):.0f}–{np.degrees(K.PSI_M+K.PSI_A):.0f}°，"
    f"而现在的四杆只走 76–136° —— 要换构型得重新设计四杆，不是改控制律就行",
    fontsize=9.2, color=MUTED, ha="left", va="top", linespacing=1.6)

# --- 动态
rockB, = axA.plot([], [], color="#6f6f6f", lw=5.5, solid_capstyle="round", zorder=5)
parOD, = axA.plot([], [], color="#a08b6a", lw=3.2, solid_capstyle="round", zorder=4)
parCE, = axA.plot([], [], color="#a08b6a", lw=3.2, solid_capstyle="round", zorder=4)
parED, = axA.plot([], [], color="#c9b493", lw=2.2, solid_capstyle="round", zorder=4)
stem,  = axA.plot([], [], color="#7b3f00", lw=3.6, solid_capstyle="round", zorder=6)
plate  = axA.add_patch(plt.Polygon([[0,0],[1,0],[1,1]], closed=True, color=BLUE, zorder=8))
joints,= axA.plot([], [], "o", color="#3f3f3f", ms=6.5, zorder=9)
hip,   = axA.plot([], [], "o", color=INK, ms=10, zorder=10)
piv,   = axA.plot([], [], "o", color=GOLD, ms=10, mec=INK, mew=1.4, zorder=12)
midm,  = axA.plot([], [], "x", color=RED, ms=8, mew=2, zorder=12)
mB, = axB.plot([], [], "o", color=INK, ms=5.5, zorder=10)
mC, = axC.plot([], [], "o", color=BLUE, ms=6, zorder=10)
dyn = []

def draw(fi):
    global dyn
    for a in dyn:
        try: a.remove()
        except Exception: pass
    dyn = []
    t = (fi/(FPS*SLOW)) % T
    s = K.state(t); lk = s["lk"]; P = s["P"]; d = s["d̂"]
    O2, B, C, D, E = K.O2, lk["B"], lk["C"], lk["D"], lk["E"]
    rockB.set_data([B[0], C[0]], [B[1], C[1]])
    parOD.set_data([O2[0], D[0]], [O2[1], D[1]])
    parCE.set_data([C[0], E[0]], [C[1], E[1]])
    parED.set_data([E[0], D[0]], [E[1], D[1]])
    stem.set_data([E[0], P[0]], [E[1], P[1]])
    LE = P + 0.25*K.HEAD*d; TE = P - 0.75*K.HEAD*d      # d = 机头方向（指向来流）
    nrm = np.array([-d[1], d[0]])*K.THK/2
    plate.set_xy([LE+nrm, TE+nrm, TE-nrm, LE-nrm])
    joints.set_data([B[0], C[0], D[0], E[0]], [B[1], C[1], D[1], E[1]])
    hip.set_data([O2[0]], [O2[1]]); piv.set_data([P[0]], [P[1]])
    mid = P - 0.25*K.HEAD*d
    midm.set_data([mid[0]], [mid[1]])

    Wv = s["ŵ"]*0.040
    a0 = arrow(axA, P - Wv, Wv, MUTED, lw=1.9, ms=12)
    a1 = arrow(axA, P, s["n̂"]*s["L"]*SC, BLUE, lw=2.8, ms=16)
    a2 = arrow(axA, P, s["ŵ"]*s["D"]*SC, GOLD, lw=1.9, ms=12)
    a3 = arrow(axA, [0.050, -0.086], [s["Fx"]*SC*2.6, 0], RED, lw=3.2, ms=18)
    t0 = axA.text(*(P - Wv*1.2), "来流 W", color=MUTED, fontsize=9, ha="right", va="center")
    t1 = axA.text(*(P + s["n̂"]*s["L"]*SC*1.15), "升力 L", color=BLUE, fontsize=10, weight="bold",
                  ha="center", va="center")
    t2 = axA.text(0.050 + s["Fx"]*SC*2.6 + 0.008, -0.086, f"推力 {s['Fx']*1e3:+.1f} mN",
                  color=RED, fontsize=10.5, weight="bold", ha="left", va="center")
    t3 = axA.text(P[0]+0.002, P[1]+0.010, "枢轴 1/4 弦", color=GOLD, fontsize=8.5, weight="bold")
    t4 = axA.text(mid[0]-0.002, mid[1]-0.011, "中点", ha="center", color=RED, fontsize=8)
    t5 = axA.text(-0.059, -0.093,
                  f"髋 ψ = {np.degrees(s['psi']):.0f}°　脚踝 θ = {np.degrees(s['th'])-TH_MID:+.0f}°　"
                  f"攻角 α = {np.degrees(s['al']):+.1f}°", color=INK2, fontsize=10)
    t6 = axA.text(O2[0]-0.004, O2[1]+0.007, "髋 O2", color=INK, fontsize=9, ha="right")
    dyn += [a0, a1, a2, a3, t0, t1, t2, t3, t4, t5, t6]
    mB.set_data([t/T], [np.degrees(s["al"])]); mC.set_data([t/T], [s["Fx"]*1e3])
    return []

anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000/FPS, blit=False)
out = os.path.join(HERE, "liftbased_linkage.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=4500,
          extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out)
