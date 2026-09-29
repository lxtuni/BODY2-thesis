#!/usr/bin/env python3
"""海豚尾鳍式运动轨迹（沉浮 heave + 俯仰 pitch），用展开态折扇脚蹼。
只考虑舵机限制，不含连杆。准定常受力仅作示意，非 CFD。
→ fluke_trajectory.mp4
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, animation
from matplotlib.patches import FancyArrowPatch, Wedge

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
RED, BLUE, GOLD, GREEN, WATER, PUR = "#c2185b", "#2a78d6", "#d4a017", "#16a085", "#dbe9f6", "#7b3f9d"

# ---------------- 脚蹼：展开态折扇（120°，半径 = 桨头长 42 mm）
RHO   = 1000.0
FAN   = np.radians(120.0)
R     = 0.042                       # 扇半径 = 弦长 chord
SPAN  = 2*R*np.sin(FAN/2)           # 最大展长 = 弦的两端
S     = 0.5*FAN*R**2                # 扇形面积
AR    = SPAN**2/S                   # 展弦比
THK   = 0.0056
PIVOT = 0.25                        # 俯仰轴距前缘的弦比例
CLA   = 2*np.pi*AR/(2+np.sqrt(AR**2+4))
AL_ST = np.radians(16.0)

# ---------------- 运动（海豚尾鳍：俯仰领先沉浮 90°）
U    = 0.25                         # 来流 [m/s]
F    = 1.0; T = 1/F; OM = 2*np.pi*F
H0   = 0.035                        # 沉浮半幅 → 峰峰 70 mm
CHI  = 0.70                         # 羽化系数 θ0 / 诱导角
SERVO_RATE = 820.0                  # 7465W @8.4V 带载角速度 [°/s]

IND0 = np.arctan(H0*OM/U)
TH0  = CHI*IND0
ST   = F*2*H0/U

def rot(a): return np.array([np.cos(a), np.sin(a)])

def state(t, u=U, f=F, h0=H0, chi=CHI):
    om = 2*np.pi*f
    ind0 = np.arctan(h0*om/max(u, 1e-6)); th0 = chi*ind0
    h    = h0*np.sin(om*t); hdot = h0*om*np.cos(om*t)
    th   = th0*np.cos(om*t)                      # 俯仰领先沉浮 90°
    W = np.array([-u, -hdot]); Wm = np.hypot(*W); ŵ = W/Wm
    al = (np.arctan2(W[1], W[0]) - th) % (2*np.pi) - np.pi
    ale = np.clip(al, -AL_ST, AL_ST)
    CL = CLA*np.sin(ale)
    CD = 0.02 + CL**2/(np.pi*AR*0.85) + 1.1*np.sin(al)**2
    q = 0.5*RHO*Wm**2*S
    n̂ = np.array([-ŵ[1], ŵ[0]])
    Fv = q*CL*n̂ + q*CD*ŵ
    return dict(h=h, hdot=hdot, th=th, th0=th0, al=al, W=W, Wm=Wm, ŵ=ŵ, n̂=n̂,
                L=q*CL, D=q*CD, Fx=Fv[0], Fz=Fv[1],
                Pin=-(Fv[1]*hdot))                # 输入功率（沉浮方向做的功）

tt = np.linspace(0, T, 721)
S_ = [state(t) for t in tt]
FX = np.array([s["Fx"] for s in S_]); PIN = np.array([s["Pin"] for s in S_])
AL = np.degrees([s["al"] for s in S_]); TH = np.degrees([s["th"] for s in S_])
HH = np.array([s["h"] for s in S_])
ETA = FX.mean()*U/max(PIN.mean(), 1e-9)
RATE_PK = np.degrees(TH0*OM)

print(f"脚蹼（展开折扇 120°, R={R*1e3:.0f} mm）：弦 {R*1e3:.0f} mm，展 {SPAN*1e3:.1f} mm，"
      f"面积 {S*1e4:.1f} cm²，展弦比 AR = {AR:.2f}（对比 42×44 矩形板 AR = 1.05）")
print(f"运动：U = {U} m/s，f = {F:.1f} Hz，沉浮峰峰 {2*H0*1e3:.0f} mm → St = {ST:.3f}")
print(f"诱导角 ±{np.degrees(IND0):.1f}°，俯仰 ±{np.degrees(TH0):.1f}°，攻角 ±{np.abs(AL).max():.1f}°")
print(f"俯仰峰值角速度 {RATE_PK:.0f} °/s（舵机 {SERVO_RATE:.0f} °/s → 余量 {SERVO_RATE/RATE_PK:.1f}×）")
print(f"准定常：平均推力 {FX.mean()*1e3:+.1f} mN，Froude 效率 {ETA*100:.0f} %（示意值）")

# ---------------- 画布
FPS, SLOW, CYC = 30, 3, 2
NF = int(CYC*T*FPS*SLOW)
fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(3, 2, width_ratios=[1.35, 1], height_ratios=[1, 1, 1],
                      left=0.045, right=0.972, top=0.815, bottom=0.065, wspace=0.18, hspace=0.62)
axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 1]); axD = fig.add_subplot(gs[2, 1])

def style(ax, water=True):
    ax.set_facecolor(WATER if water else SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if water:
        for s in ("bottom", "left"): ax.spines[s].set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.grid(color=GRID, lw=0.8); ax.tick_params(colors=INK2, labelsize=8.5)
        for s in ("bottom", "left"): ax.spines[s].set_color(GRID)
style(axA); style(axB); style(axC, False); style(axD, False)
axA.set_xlim(-0.085, 0.105); axA.set_ylim(-0.095, 0.085); axA.set_aspect("equal")
axB.set_xlim(-0.085, 0.085); axB.set_ylim(-0.052, 0.052); axB.set_aspect("equal")
SC = 0.040/max(np.hypot(FX, [s["Fz"] for s in S_]))

def arrow(ax, p0, v, c, lw=2.4, ms=15, z=9, ls="-"):
    a = FancyArrowPatch(p0, (p0[0]+v[0], p0[1]+v[1]), arrowstyle="-|>", mutation_scale=ms,
                        color=c, lw=lw, zorder=z, linestyle=ls, shrinkA=0, shrinkB=0)
    ax.add_patch(a); return a

for xx in np.linspace(-0.08, 0.075, 6):
    axA.add_patch(FancyArrowPatch((xx+0.016, 0.074), (xx, 0.074), arrowstyle="-|>",
                                  mutation_scale=9, color="#9fc3e6", lw=1.1, zorder=1))
axA.text(0.010, 0.079, f"来流 U = {U} m/s（机体向右前进）", ha="center", fontsize=9.5, color="#4f86c6")
axA.plot(np.zeros(60), np.linspace(-H0, H0, 60), color=MUTED, lw=1.2, ls=(0, (3, 2)), zorder=2)
axA.text(-0.080, 0.052, "淡蓝曲线 = 水坐标系里脚蹼走出的波形轨迹", fontsize=9, color="#7fb0dd")
axA.set_title("A  海豚尾鳍式运动：沉浮 + 俯仰（俯仰领先 90°），推力 = 升力的前向分量",
              loc="left", fontsize=12.5, color=INK, pad=9)

# ---- B 平面形状（俯视）
axB.add_patch(Wedge((-R/2, 0), R, -60, 60, color=BLUE, alpha=0.85, zorder=4))
axB.plot([-R/2, -R/2+R*np.cos(np.radians(60))], [0, R*np.sin(np.radians(60))],
         color=INK, lw=0.8, ls=(0, (2, 2)), zorder=5)
axB.annotate("", (-R/2, -0.036), (-R/2+R, -0.036), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.2))
axB.text(-R/2+R/2, -0.040, f"弦 chord = {R*1e3:.0f} mm", ha="center", va="top", fontsize=9, color=INK2)
axB.annotate("", (R/2, -SPAN/2), (R/2, SPAN/2), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.2))
axB.text(R/2+0.005, 0, f"展 span\n{SPAN*1e3:.0f} mm", fontsize=9, color=INK2, va="center")
axB.plot(-R/2 + PIVOT*R, 0, "o", color=GOLD, ms=10, mec=INK, mew=1.5, zorder=8)

axB.text(0.012, 0.96, f"120° 展开态折扇 · 面积 {S*1e4:.1f} cm² · AR = {AR:.2f} · 俯仰轴 1/4 弦（金点）",
         transform=axB.transAxes, fontsize=9.0, color=INK, weight="bold", va="top")
axB.set_title("B  脚蹼平面形状（俯视，尖端朝前 = 三角翼）", loc="left", fontsize=10.5, color=INK, pad=6)

# ---- C 角度与速率
axC.plot(tt/T, np.degrees(np.arctan2(-np.array([s['hdot'] for s in S_]), -U))*0 +
         np.degrees([np.arctan(s['hdot']/U) for s in S_]), color=MUTED, lw=1.6, ls=(0, (4, 2)),
         label=f"诱导角 ±{np.degrees(IND0):.0f}°")
axC.plot(tt/T, TH, color=PUR, lw=2.1, label=f"俯仰 θ ±{np.degrees(TH0):.0f}°")
axC.plot(tt/T, AL, color=RED, lw=2.1, label=f"攻角 α ±{np.abs(AL).max():.0f}°")
axC.plot(tt/T, HH*1e3/2, color=GREEN, lw=1.6, label=f"沉浮 h（÷2, mm）±{H0*1e3:.0f}")
axC.axhline(0, color=INK2, lw=0.9); axC.set_xlim(0, 1)
axC.set_xlabel("周期相位 τ", color=INK2, fontsize=9); axC.set_ylabel("角度 [°]", color=INK2, fontsize=9)
axC.legend(frameon=False, fontsize=7.6, loc="lower left", labelcolor=INK2, ncol=2)
axC.set_title("C  运动角度", loc="left", fontsize=10.5, color=INK, pad=6)

# ---- D 舵机可行域
ff = np.linspace(0.3, 4.0, 220)
for h0_, c_, ls_ in ((0.025, "#9ec5e8", "-"), (0.035, BLUE, "-"), (0.050, "#16558f", "-")):
    ind = np.arctan(h0_*2*np.pi*ff/U); rate = np.degrees(CHI*ind*2*np.pi*ff)
    axD.plot(ff, rate, color=c_, lw=2.0, ls=ls_, label=f"沉浮半幅 {h0_*1e3:.0f} mm")
axD.axhline(SERVO_RATE, color=RED, lw=1.8, ls=(0, (4, 3)))
axD.text(3.95, SERVO_RATE*1.06, f"7465W 带载 {SERVO_RATE:.0f} °/s", color=RED, fontsize=8.5, ha="right")
f_lo, f_hi = 0.2*U/(2*H0), 0.4*U/(2*H0)
axD.axvspan(f_lo, f_hi, color="#eaf3fb", lw=0, zorder=0)
axD.text((f_lo+f_hi)/2, 470, "St 0.2–0.4\n高效区", fontsize=8, color="#4f86c6", ha="center", va="top")
axD.plot([F], [RATE_PK], "o", color=INK, ms=8, zorder=10)
axD.text(F+0.12, RATE_PK, f"本工况\n{RATE_PK:.0f} °/s", fontsize=8.5, color=INK, va="center")
axD.set_xlim(0.3, 4.0); axD.set_ylim(0, 1250)
axD.set_xlabel("拍动频率 f [Hz]", color=INK2, fontsize=9)
axD.set_ylabel("俯仰峰值角速度 [°/s]", color=INK2, fontsize=9)
axD.legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK2)
axD.set_title("D  舵机是不是瓶颈？（结论：在高效区里不是）", loc="left", fontsize=10.5, color=INK, pad=6)

fig.text(0.045, 0.962, "海豚尾鳍式运动轨迹 —— 展开态折扇脚蹼，只看舵机限制，不含连杆",
         fontsize=15, color=INK, ha="left")
fig.text(0.045, 0.918,
    f"脚蹼 = 120° 展开折扇：弦 {R*1e3:.0f} mm、展 {SPAN*1e3:.1f} mm、面积 {S*1e4:.1f} cm²、"
    f"AR = {AR:.2f}（矩形板只有 1.05，这是换用展开态最实在的好处）\n"
    f"工况 U = {U} m/s，f = {F:.1f} Hz，沉浮峰峰 {2*H0*1e3:.0f} mm → St = {ST:.2f}（高效区 0.2–0.4 ✓）；"
    f"俯仰领先沉浮 90°，羽化系数 χ = {CHI}，攻角 ±{np.abs(AL).max():.0f}°（未失速）\n"
    f"俯仰峰值角速度仅 {RATE_PK:.0f} °/s，7465W 有 {SERVO_RATE:.0f} °/s —— 余量 {SERVO_RATE/RATE_PK:.1f} 倍。"
    f"受力为准定常示意，不是 CFD 结果。",
    fontsize=9.3, color=MUTED, ha="left", va="top", linespacing=1.65)

trail, = axA.plot([], [], color=BLUE, lw=1.0, alpha=0.35, zorder=2)
foil   = axA.add_patch(plt.Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=BLUE, zorder=8))
piv,   = axA.plot([], [], "o", color=GOLD, ms=10, mec=INK, mew=1.4, zorder=12)
mC,    = axC.plot([], [], "o", color=INK, ms=5.5, zorder=10)
dyn = []
TRX, TRZ = [], []

def draw(fi):
    global dyn, TRX, TRZ
    for a in dyn:
        try: a.remove()
        except Exception: pass
    dyn = []
    t = (fi/(FPS*SLOW)) % T
    s = state(t)
    P = np.array([0.0, s["h"]])
    d = rot(s["th"])                                   # 机头方向（指向来流一侧）
    LE = P + PIVOT*R*d; TE = P - (1-PIVOT)*R*d
    n = np.array([-d[1], d[0]])*THK/2
    foil.set_xy([LE+n, TE+n, TE-n, LE-n])
    piv.set_data([P[0]], [P[1]])
    # 水坐标系里的波形轨迹（机体前进 → 轨迹被拉成正弦）
    TRX.append(-U*t if fi % (FPS*SLOW*CYC) else 0); TRZ.append(s["h"])
    tw = np.linspace(0, T, 200)
    trail.set_data(-U*(tw - t)*0.55, [state(x)["h"] for x in tw])

    Wv = s["ŵ"]*0.042
    a0 = arrow(axA, P - Wv, Wv, MUTED, lw=1.9, ms=12)
    a1 = arrow(axA, P, s["n̂"]*s["L"]*SC, BLUE, lw=2.8, ms=16)
    a2 = arrow(axA, P, s["ŵ"]*s["D"]*SC, GOLD, lw=1.9, ms=12)
    a3 = arrow(axA, [0.0, -0.082], [s["Fx"]*SC*2.6, 0], RED, lw=3.2, ms=18)
    t1 = axA.text(*(P - Wv*1.18), "相对来流 W", color=MUTED, fontsize=9, ha="right", va="center")
    t2 = axA.text(*(P + s["n̂"]*s["L"]*SC*1.14), "升力 L", color=BLUE, fontsize=10.5,
                  weight="bold", ha="center", va="center")
    t3 = axA.text(s["Fx"]*SC*2.6 + 0.008, -0.082, f"推力 {s['Fx']*1e3:+.0f} mN",
                  color=RED, fontsize=11, weight="bold", ha="left", va="center")
    t4 = axA.text(-0.082, -0.090, f"沉浮 h = {s['h']*1e3:+.0f} mm　俯仰 θ = {np.degrees(s['th']):+.0f}°　"
                                  f"攻角 α = {np.degrees(s['al']):+.1f}°", color=INK2, fontsize=10)
    dyn += [a0, a1, a2, a3, t1, t2, t3, t4]
    mC.set_data([t/T], [np.degrees(s["al"])])
    return []

anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000/FPS, blit=False)
out = os.path.join(HERE, "fluke_trajectory.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=4500,
          extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out)
