#!/usr/bin/env python3
"""两栖海豚式尾鳍腿（髋 + 膝 + 脚蹼俯仰，3 舵机，全部矢状面 y 轴）示意动画 → amphib_leg_anim.mp4
游泳：髋把大腿定成向后下 30° 的"支柱"，膝上下摆 ±27° 做沉浮（小腿向后拖，平均下倾 10°），S3 经小腿内推杆定脚蹼绝对俯仰（ψ 90°）。
行走：髋 + 膝走标准步态，S3 把脚蹼放平当脚掌。
整机：同侧前后腿反相 = 对角同相（trot），垂向力 / 纵摇 / 横摇抵消。
"""
import os, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, animation
from matplotlib.patches import Polygon, Ellipse, FancyBboxPatch, Circle, FancyArrowPatch, Arc
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER, GROUND = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa", "#efe6d2"
BLUE, RED, GOLD, GREEN, PUR, HULL, SERVO = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#c9c4b8", "#333"

# ---------------- 参数（mm，船体坐标：x 向前，z 向上，水线 z = -10）
F, H0, TH0 = 2.0, 50.0, 44.0
L1, L2 = 90.0, 110.0
TH1_SWIM = math.radians(30.0)            # 大腿：从竖直向后转 30°
G_MEAN = math.radians(10.0)              # 小腿：向后、平均下倾 10°
G0 = math.asin(H0 / L2)                  # 膝摆幅 ±27°
CM = 28.7; PIV = 0.33 * CM
LWL, BWL, DRAFT, WL = 233.0, 99.0, 38.0, -10.0
HIP_Z = -22.0
sec = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_section_NACA0021.csv"), delimiter=",", skiprows=1)
OM = 2 * math.pi * F; T = 1 / F; U = 0.25

def swim(t, ph=0.0):
    """返回 (大腿角 th1[从竖直向后], 小腿角 g[从向后水平向下], 脚蹼绝对俯仰 θ[前缘上抬为正])"""
    g = G_MEAN - G0 * math.sin(OM * t + ph)          # 先向上（z 随 g 减小而升高）
    th = math.radians(TH0) * math.sin(OM * t + ph + math.pi / 2)
    return TH1_SWIM, g, th

def leg_points(hip, th1, g):
    K = hip + L1 * np.array([-math.sin(th1), -math.cos(th1)])
    P = K + L2 * np.array([-math.cos(g), -math.sin(g)])
    return K, P

def foil_poly(P, th, scale_t=1.6):
    d = np.array([math.cos(th), math.sin(th)]); n = np.array([-d[1], d[0]])
    return np.array([P + (PIV - x * CM) * d + y * CM * scale_t * n for x, y in sec])

def alpha_deg(t, ph=0.0):
    _, g, th = swim(t, ph)
    zdot = L2 * math.cos(g) * G0 * OM * math.cos(OM * t + ph) * 1e-3     # 脚蹼竖向速度 [m/s]（向上为正）
    return math.degrees(math.atan2(zdot, U) - th)

# 行走：脚在地面的直线 + 摆动相抬腿；膝向后弯（与游泳同一弯曲方向）
GROUND_Z = HIP_Z - 165.0
def walk_foot(tau):
    """tau ∈ [0,1)：0–0.6 支撑相（脚相对身体向后走），0.6–1 摆动相（向前抬）"""
    xs, xe = 45.0, -45.0
    if tau < 0.6:
        s = tau / 0.6; return np.array([xs + (xe - xs) * s, GROUND_Z])
    s = (tau - 0.6) / 0.4
    return np.array([xe + (xs - xe) * (0.5 - 0.5 * math.cos(math.pi * s)), GROUND_Z + 32 * math.sin(math.pi * s)])

def ik(hip, P):
    """两连杆逆解，膝向后（-x）"""
    d = P - hip; r = np.linalg.norm(d); r = min(r, L1 + L2 - 1e-6)
    c2 = (r**2 - L1**2 - L2**2) / (2 * L1 * L2); a2 = math.acos(max(-1, min(1, c2)))
    a1 = math.atan2(d[1], d[0]) - math.atan2(L2 * math.sin(a2), L1 + L2 * math.cos(a2))
    K = hip + L1 * np.array([math.cos(a1), math.sin(a1)])
    # 选膝在后的分支
    K2 = hip + L1 * np.array([math.cos(math.atan2(d[1], d[0]) + math.atan2(L2 * math.sin(a2), L1 + L2 * math.cos(a2))), math.sin(math.atan2(d[1], d[0]) + math.atan2(L2 * math.sin(a2), L1 + L2 * math.cos(a2)))])
    return K if K[0] < K2[0] else K2

def servo(ax, xy, w=16, h=11, col=SERVO):
    ax.add_patch(FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h, boxstyle="round,pad=1.2", fc=col, ec=INK, lw=1.0, zorder=9))
    ax.add_patch(Circle(xy, 2.6, fc=GOLD, ec=INK, lw=0.8, zorder=10))

def hull(ax, x0=0.0, deck=True):
    ax.add_patch(Ellipse((x0, WL), LWL, DRAFT * 2, fc=HULL, ec=INK2, lw=1, zorder=2))
    if deck: ax.add_patch(FancyBboxPatch((x0 - LWL / 2 + 10, WL), LWL - 20, 28, boxstyle="round,pad=2", fc="#e4dfd3", ec=INK2, lw=1, zorder=3))

FPS, SEC, SLOW = 30, 10.0, 2.0
NF = int(SEC * FPS); DT = 1 / FPS / SLOW
fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(2, 3, height_ratios=[2.5, 1], width_ratios=[1.15, 1, 1.05], left=0.03, right=0.985, top=0.80, bottom=0.07, hspace=0.30, wspace=0.10)
axS = fig.add_subplot(gs[0, 0]); axW = fig.add_subplot(gs[0, 1]); axB = fig.add_subplot(gs[0, 2]); axC = fig.add_subplot(gs[1, :])
for ax in (axS, axW, axB):
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
axS.set_facecolor(WATER); axB.set_facecolor(WATER); axW.set_facecolor("#f6f4ef")
for s in ("top", "right"): axC.spines[s].set_visible(False)
axC.grid(color=GRID, lw=0.8); axC.tick_params(colors=INK2, labelsize=9); axC.axhline(0, color=INK2, lw=0.8)
axC.set_xlim(0, SEC / SLOW); axC.set_ylim(-75, 75); axC.set_xlabel("真实时间 t [s]（视频放慢 2×）", color=INK2); axC.set_ylabel("舵机角 [°]（游泳模式）", color=INK2)

fig.text(0.03, 0.965, "两栖海豚式尾鳍腿：髋 S1 + 膝 S2 + 脚蹼俯仰 S3（三舵机，全部 y 轴，矢状面内）—— 同一条腿既走路又当尾鳍", fontsize=16, color=INK, va="top")
fig.text(0.03, 0.92, "游泳：S1 把大腿（90 mm）固定成向后下 30° 的流线支柱，把膝放深；S2 让小腿（110 mm，向后拖）上下摆 ±27° → 脚蹼沉浮 ±50 mm；"
         "S3 经小腿内推杆定脚蹼绝对俯仰 ±44°、领先 90°。脚蹼顶端离水面 ≈ 2 弦。\n"
         "行走：S1 + S2 走标准两连杆步态（膝向后弯，与游泳同向），S3 把脚蹼放平当脚掌（前缘朝前，尾缘加胶条）。切换 = 髋转 60° + 膝折 100°，不换零件。",
         fontsize=9.6, color=MUTED, va="top", linespacing=1.6)

# ============ 面板 1：游泳（单条后腿，大图）
HIP_S = np.array([-40.0, HIP_Z])
axS.set_xlim(-270, 80); axS.set_ylim(-205, 45)
hull(axS)
axS.axhline(WL, color="#4f86c6", lw=1.2, ls=(0, (5, 3))); axS.text(85, WL + 3, "水线", fontsize=8.5, color="#4f86c6")
servo(axS, HIP_S); axS.text(HIP_S[0] + 12, HIP_S[1] + 8, "S1 髋（船内）：游泳时固定大腿角", fontsize=8.5, color=INK, weight="bold")
servo(axS, HIP_S + np.array([-28, 0])); axS.text(HIP_S[0] - 40, HIP_S[1] - 14, "S2 膝（船内，平行四边形下传）", fontsize=8.5, color=INK, weight="bold", ha="right")
axS.add_patch(FancyArrowPatch((90, -120), (60, -120), arrowstyle="-|>", mutation_scale=12, color="#4f86c6", lw=1.3)); axS.text(92, -120, "来流", fontsize=8.5, color="#4f86c6", va="center")
axS.set_title("游泳模式 —— 侧视图（后腿）：膝上下摆 = 沉浮，S3 = 俯仰", loc="left", fontsize=11, color=INK)
thighS, = axS.plot([], [], color=INK, lw=7, solid_capstyle="round", zorder=5)
shankS, = axS.plot([], [], color=INK, lw=6, solid_capstyle="round", zorder=5)
parS, = axS.plot([], [], color=GOLD, lw=2, zorder=6)          # 膝平行四边形
rodS, = axS.plot([], [], color=RED, lw=2, zorder=6)           # S3 推杆
kneeS, = axS.plot([], [], "o", color=GOLD, ms=8, mec=INK, zorder=9)
s3S = FancyBboxPatch((0, 0), 12, 9, boxstyle="round,pad=1.2", fc=SERVO, ec=INK, lw=1, zorder=9); axS.add_patch(s3S)
s3St = axS.text(0, 0, "S3 俯仰（膝部，推杆经小腿到脚蹼）", fontsize=8.5, color=INK, weight="bold")
foilS = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=BLUE, zorder=8); axS.add_patch(foilS)
pivS, = axS.plot([], [], "o", color=GOLD, ms=6, mec=INK, zorder=10)
trailS, = axS.plot([], [], color=BLUE, lw=1.3, alpha=0.4, zorder=3)
arcS = Arc((0, 0), 2 * L2, 2 * L2, color=GOLD, lw=1, ls=(0, (3, 3)), zorder=3); axS.add_patch(arcS)
txtS = axS.text(-265, -198, "", fontsize=9, color=INK2)
depS = axS.text(-265, -182, "", fontsize=9, color=RED)
trS = []

# ============ 面板 2：行走
HIP_W = np.array([-40.0, HIP_Z])
axW.set_xlim(-230, 150); axW.set_ylim(GROUND_Z - 25, 45)
hull(axW)
axW.axhspan(GROUND_Z - 25, GROUND_Z, color=GROUND, zorder=1); axW.axhline(GROUND_Z, color="#a89a78", lw=1.5)
servo(axW, HIP_W); servo(axW, HIP_W + np.array([-28, 0]))
axW.text(HIP_W[0] + 12, HIP_W[1] + 8, "S1 髋 + S2 膝：标准两连杆步态", fontsize=8.5, color=INK, weight="bold")
axW.set_title("行走模式 —— 同一条腿：脚蹼放平当脚掌", loc="left", fontsize=11, color=INK)
thighW, = axW.plot([], [], color=INK, lw=7, solid_capstyle="round", zorder=5)
shankW, = axW.plot([], [], color=INK, lw=6, solid_capstyle="round", zorder=5)
kneeW, = axW.plot([], [], "o", color=GOLD, ms=8, mec=INK, zorder=9)
foilW = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=BLUE, zorder=8); axW.add_patch(foilW)
pathW, = axW.plot([], [], color=GOLD, lw=1, ls=(0, (3, 3)), zorder=3)
txtW = axW.text(-225, GROUND_Z - 20, "", fontsize=9, color=INK2)
gtick = [axW.plot([], [], color="#a89a78", lw=1)[0] for _ in range(12)]

# ============ 面板 3：整机（前后腿反相）
axB.set_xlim(-300, 220); axB.set_ylim(-240, 45)
hull(axB)
axB.axhline(WL, color="#4f86c6", lw=1.2, ls=(0, (5, 3)))
HIP_R, HIP_F_ = np.array([-60.0, HIP_Z]), np.array([80.0, HIP_Z])
axB.set_title("整机 —— 同侧前后腿反相 = 对角同相（trot）：垂向力、纵摇、横摇抵消", loc="left", fontsize=11, color=INK)
legsB = []
for hip, ph, col in ((HIP_R, 0.0, INK), (HIP_F_, math.pi, "#666")):
    th_, = axB.plot([], [], color=col, lw=6, solid_capstyle="round", zorder=5); sh_, = axB.plot([], [], color=col, lw=5, solid_capstyle="round", zorder=5)
    fo_ = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=BLUE, zorder=8); axB.add_patch(fo_)
    tr_, = axB.plot([], [], color=BLUE, lw=1.2, alpha=0.35, zorder=3)
    servo(axB, hip, w=12, h=8)
    legsB.append(dict(hip=hip, ph=ph, th=th_, sh=sh_, fo=fo_, tr=tr_, pts=[]))
axB.text(-295, -232, "前腿脚蹼落在后腿前方 ≈ 5 弦：后腿在前腿尾迹里，相位要随速度调（Muscutt 2017）", fontsize=8.3, color=INK2)
txtB = axB.text(-295, -215, "", fontsize=9, color=INK2)

# ============ 曲线
l1, = axC.plot([], [], color=INK, lw=1.6, ls=":", label="S1 髋（固定 30°）")
l2, = axC.plot([], [], color=GOLD, lw=2, label="S2 膝摆角（相对大腿）")
l3, = axC.plot([], [], color=RED, lw=2, label="S3 俯仰（相对小腿）")
l4, = axC.plot([], [], color=BLUE, lw=1.8, ls="--", label="脚蹼绝对俯仰 θ")
l5, = axC.plot([], [], color=GREEN, lw=1.8, ls="-.", label="攻角 α")
axC.legend(frameon=False, fontsize=8.5, ncol=5, labelcolor=INK2, loc="upper right")
ts, y1, y2, y3, y4, y5 = [], [], [], [], [], []

def draw(fi):
    t = fi * DT
    # ---- 游泳
    th1, g, th = swim(t)
    K, P = leg_points(HIP_S, th1, g)
    thighS.set_data([HIP_S[0], K[0]], [HIP_S[1], K[1]]); shankS.set_data([K[0], P[0]], [K[1], P[1]])
    kneeS.set_data([K[0]], [K[1]]); pivS.set_data([P[0]], [P[1]])
    # 膝平行四边形：与大腿平行、偏置 9 mm，曲柄角 = 小腿角
    cr = 9.0; ang_c = -g + math.pi          # 曲柄方向随小腿转
    cH = HIP_S + cr * np.array([math.cos(ang_c + math.pi / 2), math.sin(ang_c + math.pi / 2)])
    cK = K + cr * np.array([math.cos(ang_c + math.pi / 2), math.sin(ang_c + math.pi / 2)])
    parS.set_data([HIP_S[0], cH[0], cK[0], K[0]], [HIP_S[1], cH[1], cK[1], K[1]])
    # S3 推杆：从膝沿小腿到脚蹼曲柄
    cF = P + 8 * np.array([-math.sin(th), math.cos(th)])
    s3c = K + 22 * np.array([-math.cos(g), -math.sin(g)])
    rodS.set_data([s3c[0], cF[0], P[0]], [s3c[1], cF[1], P[1]])
    s3S.set_bounds(s3c[0] - 6, s3c[1] - 4.5, 12, 9); s3St.set_position((s3c[0] - 30, s3c[1] + 14))
    foilS.set_xy(foil_poly(P, th, 2.0))
    arcS.center = (K[0], K[1]); arcS.theta1 = 180 + math.degrees(G_MEAN - G0); arcS.theta2 = 180 + math.degrees(G_MEAN + G0)
    for i in range(len(trS)): trS[i][0] -= U * 1e3 * DT
    trS.append([P[0] - CM, P[1]])
    while trS and trS[0][0] < -265: trS.pop(0)
    trailS.set_data([p[0] for p in trS], [p[1] for p in trS])
    al = alpha_deg(t)
    txtS.set_text(f"膝 g = {math.degrees(g):+5.1f}°　沉浮 {-(P[1] - (K[1] - L2*math.sin(G_MEAN))):+5.1f} mm　θ = {math.degrees(th):+5.1f}°　α = {al:+5.1f}°")
    top = P[1] + 0.33 * CM * math.sin(th) if th > 0 else P[1]
    depS.set_text(f"脚蹼离水面 {WL - (P[1] + 6):.0f} mm ≈ {(WL - (P[1] + 6))/CM:.1f} 弦（冲程最高点 ≈ 2.0 弦，自由面效应可忽略）")
    # ---- 行走
    tau = (t * 1.2) % 1.0
    Pw = walk_foot(tau); Kw = ik(HIP_W, Pw)
    thighW.set_data([HIP_W[0], Kw[0]], [HIP_W[1], Kw[1]]); shankW.set_data([Kw[0], Pw[0]], [Kw[1], Pw[1]]); kneeW.set_data([Kw[0]], [Kw[1]])
    foilW.set_xy(foil_poly(Pw + np.array([0, 5]), 0.0, scale_t=1.6))
    pw = np.array([walk_foot(s) for s in np.linspace(0, 1, 80)]); pathW.set_data(pw[:, 0], pw[:, 1])
    shift = (-t * 1.2 * 90.0 / 0.6) % 40
    for i, gl in enumerate(gtick):
        x = -230 + i * 40 + shift; gl.set_data([x, x + 12], [GROUND_Z - 6, GROUND_Z - 6])
    txtW.set_text("支撑相 0–60%（脚相对身体向后）· 摆动相抬 32 mm；脚掌 = 放平的脚蹼" if tau < 0.6 else "摆动相：膝先折再伸，脚蹼保持水平")
    # ---- 整机
    for L_ in legsB:
        th1b, gb, thb = swim(t, L_["ph"]); Kb, Pb = leg_points(L_["hip"], th1b, gb)
        L_["th"].set_data([L_["hip"][0], Kb[0]], [L_["hip"][1], Kb[1]]); L_["sh"].set_data([Kb[0], Pb[0]], [Kb[1], Pb[1]])
        L_["fo"].set_xy(foil_poly(Pb, thb, 1.3))
        for i in range(len(L_["pts"])): L_["pts"][i][0] -= U * 1e3 * DT
        L_["pts"].append([Pb[0] - CM, Pb[1]])
        while L_["pts"] and L_["pts"][0][0] < -295: L_["pts"].pop(0)
        L_["tr"].set_data([p[0] for p in L_["pts"]], [p[1] for p in L_["pts"]])
    txtB.set_text("前腿（浅色）与后腿（深色）反相；另一侧镜像 → 对角两腿同相")
    # ---- 曲线
    ts.append(t); y1.append(30.0); y2.append(math.degrees(g - G_MEAN)); y3.append(math.degrees(th - g)); y4.append(math.degrees(th)); y5.append(al)
    for l_, y_ in ((l1, y1), (l2, y2), (l3, y3), (l4, y4), (l5, y5)): l_.set_data(ts, y_)
    return []

anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000 / FPS, blit=False)
out = os.path.join(HERE, "amphib_leg_anim.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=5000, extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out)
