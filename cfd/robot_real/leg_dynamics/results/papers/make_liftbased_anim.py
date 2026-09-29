#!/usr/bin/env python3
"""升力型（lift-based）关节运动概念动画 vs 阻力型（drag-based）对照。
准定常（quasi-steady）受力模型，仅作机理示意，不是 CFD 结果。
→ liftbased_vs_dragbased.mp4
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, animation
from matplotlib.patches import FancyArrowPatch

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))

INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
RED, BLUE, GOLD, GREEN, WATER = "#c2185b", "#2a78d6", "#d4a017", "#16a085", "#dbe9f6"

# ---------------------------------------------------------------- 参数
RHO  = 1000.0
U    = 0.25          # 动画工况的前进速度 [m/s]
T    = 1.0           # 周期 [s]
F    = 1.0 / T
OM   = 2 * np.pi * F
H0   = 0.042         # 垂直扫动半幅 [m] → 峰峰 84 mm = Lu 的行程
C    = 0.042         # 弦长 = 板长
SPAN = 0.044         # 展长 = 板宽
THK  = 0.0056
S    = C * SPAN                      # 18.5 cm²
AR   = SPAN / C                      # ≈ 1.05
S_FE = THK * SPAN                    # 羽化后迎流面积
CHI  = 0.70                          # 羽化系数 θ0 / 诱导角
CLA  = 2*np.pi*AR / (2 + np.sqrt(AR**2 + 4))    # 低展弦比升力线斜率 (Helmbold)
CD_N = 3.30      # 平板法向力系数：按 Lu 自己 CFD 反推的非定常有效值（准定常只有 1.17）
DUTY, RAMP = 0.416, 0.15
STROKE = 2 * H0
AL_STALL = np.radians(16.0)

def lift_based(t, u=U):
    """升力型：沉浮 h(t) + 俯仰 θ(t)，俯仰领先 90°"""
    ind0 = np.arctan(H0 * OM / max(u, 1e-6))
    th0  = CHI * ind0
    h    = H0 * np.sin(OM * t)
    hdot = H0 * OM * np.cos(OM * t)
    th   = th0 * np.cos(OM * t)
    Wx, Wz = -u, -hdot                            # 板感受到的相对来流
    Wm = np.hypot(Wx, Wz); wx, wz = Wx/Wm, Wz/Wm
    al = (np.arctan2(Wz, Wx) - th) % (2*np.pi) - np.pi   # 攻角 α
    ale = np.clip(al, -AL_STALL, AL_STALL)
    CL = CLA * np.sin(ale)
    CD = 0.02 + CL**2/(np.pi*AR*0.8) + CD_N*np.sin(al)**2
    q  = 0.5*RHO*Wm**2*S
    nx, nz = -wz, wx                              # 升力方向 ⊥ W
    return dict(h=h, hdot=hdot, th=th, th0=th0, ind0=ind0, W=(Wx, Wz), Wm=Wm, al=al,
                L=q*CL, D=q*CD, n=(nx, nz), w=(wx, wz),
                Fx=q*CL*nx + q*CD*wx, Fz=q*CL*nz + q*CD*wz)

def trap_v(x, ramp):
    x = np.clip(x, 0, 1)
    if x <= ramp:      return x/ramp
    if x >= 1 - ramp:  return (1-x)/ramp
    return 1.0

def drag_based(t, u=U):
    """阻力型：板正对后划 + 羽化 90° 前收。相对水的速度才是产生力的量。"""
    tau = (t/T) % 1.0
    if tau < DUTY:                                    # 划水 power phase
        vpk = STROKE/(DUTY*T*(1-RAMP)); vp = trap_v(tau/DUTY, RAMP)*vpk
        vrel = vp - u                                 # >0 才有推力：桨必须追得上船
        Fx = 0.5*RHO*CD_N*S*vrel*abs(vrel)
        return dict(phase="power", Fx=Fx, beta=0.0, v=vp, vrel=vrel,
                    x=STROKE/2 - STROKE*np.clip(tau/DUTY, 0, 1))
    tr = (tau-DUTY)/(1-DUTY)                          # 回程 recovery phase
    vpk = STROKE/((1-DUTY)*T*(1-RAMP)); vr = trap_v(tr, RAMP)*vpk
    vrel = vr + u                                     # 前收 + 船前进 → 相对速度相加
    Fx = -0.5*RHO*CD_N*S_FE*vrel**2
    return dict(phase="recovery", Fx=Fx, beta=90.0, v=vr, vrel=vrel,
                x=-STROKE/2 + STROKE*tr)

def cycle_mean(fn, u, n=720):
    return np.mean([fn(t, u)["Fx"] for t in np.linspace(0, T, n, endpoint=False)])

# ---------------------------------------------------------------- 预算
NT = 600
tt = np.linspace(0, T, NT, endpoint=False)
LB = [lift_based(t) for t in tt];  DB = [drag_based(t) for t in tt]
FxL = np.array([d["Fx"] for d in LB]); FxD = np.array([d["Fx"] for d in DB])
ST  = F*2*H0/U
IND0, TH0 = LB[0]["ind0"], LB[0]["th0"]
UU  = np.linspace(0.02, 0.50, 80)
ML  = np.array([cycle_mean(lift_based, u) for u in UU])
MD  = np.array([cycle_mean(drag_based, u) for u in UU])
cross = UU[np.argmin(np.abs(ML - MD))]
u_band = (F*2*H0/0.4, F*2*H0/0.2)                 # St 0.2–0.4 对应的 U 区间

print(f"St = {ST:.3f}  诱导角 {np.degrees(IND0):.1f}°  θ0 {np.degrees(TH0):.1f}°  "
      f"α幅值 {np.degrees(max(abs(d['al']) for d in LB)):.1f}°  AR {AR:.2f}")
print(f"@U={U}: 升力型 {FxL.mean()*1e3:+.1f} mN   阻力型 {FxD.mean()*1e3:+.1f} mN")
print(f"阻力型桨速峰值 {STROKE/(DUTY*T*(1-RAMP)):.3f} m/s → 速度天花板")
print(f"交叉点 U ≈ {cross:.3f} m/s；St 0.2–0.4 对应 U = {u_band[0]:.2f}–{u_band[1]:.2f} m/s")
print(f"脚踝角速度：升力型 {np.degrees(TH0*OM):.0f} °/s vs 阻力型羽化 {90/0.11:.0f} °/s")

# ---------------------------------------------------------------- 画布
FPS, SLOW, CYC = 30, 3, 2
NF = int(CYC*T*FPS*SLOW)
fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(3, 2, width_ratios=[1.28, 1], height_ratios=[1.05, 0.62, 0.62],
                      left=0.045, right=0.972, top=0.805, bottom=0.065, wspace=0.16, hspace=0.42)
axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 1]); axD = fig.add_subplot(gs[2, 1])

def style(ax, water=True, box=False):
    ax.set_facecolor(WATER if water else SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if water:
        for s in ("bottom", "left"): ax.spines[s].set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.grid(color=GRID, lw=0.8); ax.tick_params(colors=INK2, labelsize=8.5)
        for s in ("bottom", "left"): ax.spines[s].set_color(GRID)
style(axA); style(axB); style(axC, water=False); style(axD, water=False)
axA.set_xlim(-0.118, 0.128); axA.set_ylim(-0.148, 0.098); axA.set_aspect("equal")
axB.set_xlim(-0.105, 0.125); axB.set_ylim(-0.098, 0.050); axB.set_aspect("equal")

SC  = 0.052/max(np.hypot(FxL, [d["Fz"] for d in LB]))
SCB = 0.048/max(np.abs(FxD))

def plate_poly(cx, cz, ang, c=C, th=THK):
    pts = np.array([[-c/2,-th/2],[c/2,-th/2],[c/2,th/2],[-c/2,th/2]])
    R = np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]])
    return pts @ R.T + np.array([cx, cz])

def arrow(ax, p0, vec, color, lw=2.4, ms=15, z=9, ls="-"):
    a = FancyArrowPatch(p0, (p0[0]+vec[0], p0[1]+vec[1]), arrowstyle="-|>", mutation_scale=ms,
                        color=color, lw=lw, zorder=z, linestyle=ls, shrinkA=0, shrinkB=0)
    ax.add_patch(a); return a

# ---- 静态装饰
for ax, y0, x1 in ((axA, 0.086, 0.105), (axB, 0.041, 0.105)):
    for xx in np.linspace(ax.get_xlim()[0]+0.015, x1, 6):
        ax.add_patch(FancyArrowPatch((xx+0.020, y0), (xx, y0), arrowstyle="-|>",
                                     mutation_scale=9, color="#9fc3e6", lw=1.1, zorder=1))
axA.text(0.005, 0.094, f"水相对机体的来流（机体以 U = {U} m/s 向右前进）",
         ha="center", fontsize=10, color="#4f86c6")
axA.set_title("A  升力型 lift-based：板上下扫动 + 始终保持攻角，推力 = 升力的前向分量",
              loc="left", fontsize=12.5, color=INK, pad=9)
axB.set_title("B  阻力型 drag-based（你现在的 D2）：板正对后划 + 羽化前收",
              loc="left", fontsize=11, color=INK, pad=7)
axC.set_title("C  瞬时推力波形（各自归一化，看形状不看大小）", loc="left", fontsize=10.5, color=INK, pad=6)
axD.set_title("D  平均推力随前进速度 U 的走势（归一化；绝对值两者不可直接比）",
              loc="left", fontsize=10.5, color=INK, pad=6)

nL, nD = np.max(np.abs(FxL)), np.max(np.abs(FxD))
axC.plot(tt/T, FxL/nL, color=BLUE, lw=2.1, label="升力型：两个半程都产生推力")
axC.plot(tt/T, FxD/nD, color=RED,  lw=2.1, label="阻力型：一推一亏")
axC.axhline(0, color=INK2, lw=0.9); axC.set_xlim(0, 1); axC.set_ylim(-0.75, 1.15)
axC.set_xlabel("周期相位 τ", color=INK2, fontsize=9)
axC.set_ylabel("推力 / 各自峰值", color=INK2, fontsize=9)
axC.legend(frameon=False, fontsize=8.5, loc="lower left", labelcolor=INK2)

axD.axvspan(u_band[0], u_band[1], color="#eaf3fb", lw=0, zorder=0)
axD.text(np.mean(u_band), axD.get_ylim()[1], "", fontsize=8)
axD.plot(UU, ML/ML[0], color=BLUE, lw=2.2, label="升力型")
axD.plot(UU, MD/MD[0], color=RED,  lw=2.2, label="阻力型")
axD.axhline(0, color=INK2, lw=0.9)
axD.axvline(U, color=INK, lw=1.0, ls=(0, (3, 2)))
axD.set_xlim(UU[0], UU[-1]); axD.set_xlabel("前进速度 U [m/s]", color=INK2, fontsize=9)
axD.set_ylabel("平均推力 / 各自低速值", color=INK2, fontsize=9)
axD.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK2)
yl = axD.get_ylim()
axD.text(u_band[1]-0.005, 1.52, "St 0.2–0.4 升力型高效区", fontsize=8, color="#4f86c6",
         ha="right", va="top")
axD.text(U-0.008, -1.95, "动画工况", fontsize=8, color=INK, ha="right")
vpk_ = STROKE/(DUTY*T*(1-RAMP))
axD.annotate(f"天花板 U ≈ 桨速 {vpk_:.2f} m/s\n阻力型推力归零", (vpk_*0.86, 0), xytext=(0.30, -1.4),
             fontsize=8, color=RED, arrowprops=dict(arrowstyle="->", color=RED, lw=1))
axD.set_ylim(-2.2, 1.6)

fig.text(0.045, 0.962, "升力型 vs 阻力型：同一条腿、同样的 84 mm 行程与周期，只改运动方式",
         fontsize=15, color=INK, ha="left")
fig.text(0.045, 0.918,
   f"参数 U = {U} m/s，f = {F:.1f} Hz，扫动峰峰值 {2*H0*1e3:.0f} mm，板 {C*1e3:.0f} × {SPAN*1e3:.0f} mm　|　"
   f"Strouhal St = {ST:.2f}　展弦比 AR = {AR:.2f}（低，对升力型不利）\n"
   f"诱导角 ±{np.degrees(IND0):.0f}°，俯仰幅值 theta0 = ±{np.degrees(TH0):.0f}°（χ = {CHI}），攻角 α ≈ ±{np.degrees(IND0-TH0):.0f}°　|　"
   f"脚踝角速度：升力型 {np.degrees(TH0*OM):.0f} °/s vs 阻力型羽化 {90/0.11:.0f} °/s —— 升力型对舵机反而更轻松\n"
   "模型：准定常受力，阻力型的法向力系数已按 Lu 的 CFD 标定，升力型用低展弦比升力线理论 —— 两者绝对值不可直接比较，图 C/D 只看形状与趋势",
   fontsize=9.2, color=MUTED, ha="left", va="top", linespacing=1.6)

# ---- 动态元素
railA, = axA.plot([], [], color="#b9b6b0", lw=3, solid_capstyle="round", zorder=2)
trailA,= axA.plot([], [], color=BLUE, lw=1.1, ls=(0, (3, 2)), alpha=0.6, zorder=2)
strutA,= axA.plot([], [], color="#7b7b7b", lw=5, solid_capstyle="round", zorder=5)
plateA = axA.add_patch(plt.Polygon(plate_poly(0, 0, 0), closed=True, color=BLUE, zorder=8))
hingeA,= axA.plot([], [], "o", color=GOLD, ms=9, mec=INK, zorder=10)
railB, = axB.plot([], [], color="#b9b6b0", lw=3, solid_capstyle="round", zorder=2)
strutB,= axB.plot([], [], color="#7b7b7b", lw=5, solid_capstyle="round", zorder=5)
plateB = axB.add_patch(plt.Polygon(plate_poly(0, 0, 0), closed=True, color=RED, zorder=8))
hingeB,= axB.plot([], [], "o", color=GOLD, ms=8, mec=INK, zorder=10)
mC,  = axC.plot([], [], "o", color=BLUE, ms=6.5, zorder=10)
mC2, = axC.plot([], [], "o", color=RED,  ms=6.5, zorder=10)
STRUT = 0.050
dyn = []

def draw(fi):
    global dyn
    for a in dyn:
        try: a.remove()
        except Exception: pass
    dyn = []
    t  = (fi/(FPS*SLOW)) % T
    d  = lift_based(t); db = drag_based(t)

    # ---------- A
    p = np.array([0.0, d["h"] - STRUT])
    railA.set_data([0, 0], [-H0, H0])
    trailA.set_data(np.zeros(60), np.linspace(-H0, H0, 60) - STRUT)
    strutA.set_data([0, 0], [d["h"], p[1]])
    plateA.set_xy(plate_poly(p[0], p[1], d["th"]))
    hingeA.set_data([p[0]], [p[1]])

    Wv = np.array(d["W"])/d["Wm"]*0.050
    a0 = arrow(axA, p - Wv, Wv, MUTED, lw=2.0, ms=13)
    t0 = axA.text(*(p - Wv*1.18), "相对来流 W", color=MUTED, fontsize=9.5, ha="right", va="center")
    Lv = np.array(d["n"])*d["L"]*SC; Dv = np.array(d["w"])*d["D"]*SC
    Fv = np.array([d["Fx"], d["Fz"]])*SC
    a1 = arrow(axA, p, Lv, BLUE, lw=2.9, ms=17)
    a2 = arrow(axA, p, Dv, GOLD, lw=2.0, ms=13)
    a3 = arrow(axA, p, Fv, GREEN, lw=1.6, ms=12, ls=(0, (3, 2)))
    a4 = arrow(axA, [0.0, -0.122], [d["Fx"]*SC*2.4, 0], RED, lw=3.4, ms=19)
    ch = np.array([np.cos(d["th"]), np.sin(d["th"])])*C*0.85
    a5 = axA.plot([p[0]-ch[0], p[0]+ch[0]], [p[1]-ch[1], p[1]+ch[1]],
                  color=INK, lw=0.9, ls=(0, (2, 2)), zorder=7)[0]
    t1 = axA.text(p[0]+Lv[0]*1.14, p[1]+Lv[1]*1.14, "升力 L", color=BLUE, fontsize=10.5,
                  ha="center", va="center", weight="bold")
    t2 = axA.text(p[0]+Dv[0]*1.7, p[1]+Dv[1]*1.7, "阻力 D", color=GOLD, fontsize=9, ha="center", va="center")
    t3 = axA.text(p[0]+Fv[0]*1.12, p[1]+Fv[1]*1.12, "合力", color=GREEN, fontsize=8.5, ha="center", va="center")
    t4 = axA.text(d["Fx"]*SC*2.4 + (0.010 if d["Fx"] > 0 else -0.010), -0.122,
                  f"推力 {d['Fx']*1e3:+.0f} mN", color=RED, fontsize=11, weight="bold",
                  ha="left" if d["Fx"] > 0 else "right", va="center")
    t5 = axA.text(-0.112, -0.140,
                  f"攻角 α = {np.degrees(d['al']):+.1f}°　俯仰 θ = {np.degrees(d['th']):+.1f}°　"
                  f"扫动速度 dh/dt = {d['hdot']:+.2f} m/s", color=INK2, fontsize=10, ha="left")
    dyn += [a0, t0, a1, a2, a3, a4, a5, t1, t2, t3, t4, t5]

    # ---------- B
    bx, bz = db["x"], -0.030
    railB.set_data([-STROKE/2, STROKE/2], [0.012, 0.012])
    strutB.set_data([bx, bx], [0.012, bz])
    plateB.set_xy(plate_poly(bx, bz, np.radians(90.0 - db["beta"])))
    plateB.set_color(RED if db["phase"] == "power" else BLUE)
    hingeB.set_data([bx], [bz])
    b1 = arrow(axB, [0.0, -0.066], [db["Fx"]*SCB, 0], RED, lw=3.0, ms=17)
    b2 = axB.text(0.01, -0.089,
                  f"{'划水（板正对）' if db['phase']=='power' else '回程（羽化 90°）'}"
                  f"　桨相对水 {db['vrel']:+.2f} m/s　推力 {db['Fx']*1e3:+.0f} mN",
                  color=INK2, fontsize=9.5, ha="center")
    dyn += [b1, b2]

    mC.set_data([t/T], [d["Fx"]*1e3]); mC2.set_data([t/T], [db["Fx"]*1e3])
    return []

anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000/FPS, blit=False)
out = os.path.join(HERE, "liftbased_vs_dragbased.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=4500,
          extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out)
