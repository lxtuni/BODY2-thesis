#!/usr/bin/env python3
"""尾鳍腿的机械结构示意动画（含舵机位置）→ mechanism_anim.mp4
方案 F（鱼尾式，推荐）：S1 髋滚转舵机（x 轴）摆动竖直流线腿 → 脚蹼横向沉浮；S2 装在髋支架上，输出轴沿腿，
    通过腿内扭转轴让脚蹼绕腿轴俯仰。脚蹼竖放（展向 z），深度恒定。
方案 D（海豚式）：龙骨固定支柱把髋放深；S1 髋俯仰舵机（y 轴）摆动向后拖的腿 → 脚蹼上下沉浮；
    S2 通过平行四边形连杆控制脚蹼绝对俯仰角。
运动学 = 可信最优：f 2 Hz、h0 ±50 mm、θ0 ±44°、ψ 90°；腿长 110 mm → 摆角 ±27°。
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
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
BLUE, RED, GOLD, GREEN, PUR, HULL, SERVO = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#c9c4b8", "#333"

# ---------------- 参数（mm）
F, H0, TH0, PSI = 2.0, 50.0, 44.0, 90.0
L = 110.0; PHI0 = math.degrees(math.asin(H0 / L))
CM, B = 28.7, 86.0                      # 展中弦、展长
LWL, BWL, DRAFT, WL = 233.0, 99.0, 38.0, -10.0
pl = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_planform_mm.csv"), delimiter=",", skiprows=1)   # x 向后、y 展向
sec = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_section_NACA0021.csv"), delimiter=",", skiprows=1)
PIV = 0.33 * CM
OM = 2 * math.pi * F; T = 1 / F
U = 0.25

def phi(t): return math.radians(PHI0) * math.sin(OM * t)            # 髋摆角
def theta(t): return math.radians(TH0) * math.sin(OM * t + math.radians(PSI))   # 脚蹼绝对俯仰
def alpha(t):
    hd = H0 * 1e-3 * OM * math.cos(OM * t)
    return math.degrees(math.atan2(hd, U) - theta(t))

def rot(p, a): c, s = math.cos(a), math.sin(a); return np.array([[c, -s], [s, c]]) @ p
def servo(ax, xy, w=16, h=11, lab="", axis="o", col=SERVO):
    ax.add_patch(FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h, boxstyle="round,pad=1.2", fc=col, ec=INK, lw=1.0, zorder=9))
    if axis == "o": ax.add_patch(Circle(xy, 2.6, fc=GOLD, ec=INK, lw=0.8, zorder=10))
    ax.text(xy[0], xy[1] - h / 2 - 3, lab, ha="center", va="top", fontsize=8.5, color=INK, weight="bold", zorder=10)

FPS, SEC, SLOW = 30, 8.0, 2.0
NF = int(SEC * FPS); DT = 1 / FPS / SLOW
fig = plt.figure(figsize=(16, 9), facecolor=SURF)
gs = fig.add_gridspec(2, 3, height_ratios=[2.4, 1], left=0.03, right=0.985, top=0.82, bottom=0.07, hspace=0.32, wspace=0.12)
axR = fig.add_subplot(gs[0, 0]); axT = fig.add_subplot(gs[0, 1]); axD = fig.add_subplot(gs[0, 2]); axS = fig.add_subplot(gs[1, :])
for ax in (axR, axT, axD):
    ax.set_facecolor(WATER); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
for s in ("top", "right"): axS.spines[s].set_visible(False)
axS.grid(color=GRID, lw=0.8); axS.tick_params(colors=INK2, labelsize=9); axS.axhline(0, color=INK2, lw=0.8)
axS.set_xlim(0, SEC / SLOW); axS.set_ylim(-60, 60); axS.set_xlabel("真实时间 t [s]（视频放慢 2×）", color=INK2); axS.set_ylabel("角度 [°]", color=INK2)

fig.text(0.03, 0.965, "尾鳍腿的机械结构：两种“双连杆”—— 第二连杆就是脚蹼本身（髋摆动 + 脚蹼俯仰，2 个舵机，不需要膝）", fontsize=16.5, color=INK, va="top")
fig.text(0.03, 0.92, "运动学 = 可信最优：f 2 Hz · 沉浮 ±50 mm · 俯仰 ±44° · ψ 90°；腿长 110 mm → 髋摆 ±27°。摇臂弧线 vs 直线沉浮：UVLM 推力 +6–8%、η 不变 → 膝关节没有必要。\n"
         "S1 = 髋舵机（产生沉浮）；S2 = 俯仰舵机（脚蹼攻角）。方案 F 的 S2 装在髋支架上、输出轴沿腿，通过腿内扭转轴转脚蹼；方案 D 的 S2 在船体内，经平行四边形连杆把绝对俯仰角传到脚蹼。",
         fontsize=9.6, color=MUTED, va="top", linespacing=1.6)

# ============ 方案 F 后视图（y–z）
HIP_F = np.array([BWL / 2 + 8, -18.0])            # 髋在船侧、水线下 8 mm
axR.set_xlim(-30, 190); axR.set_ylim(-200, 40)
axR.add_patch(Ellipse((0, WL), BWL, DRAFT * 2, fc=HULL, ec=INK2, lw=1, zorder=2))
axR.add_patch(FancyBboxPatch((-BWL / 2, WL), BWL, 30, boxstyle="round,pad=2", fc="#e4dfd3", ec=INK2, lw=1, zorder=3))
axR.text(0, WL - DRAFT / 2, "船体（后视）", ha="center", va="center", fontsize=8.5, color=INK2, zorder=4)
axR.axhline(WL, color="#4f86c6", lw=1.2, ls=(0, (5, 3))); axR.text(-28, WL + 3, "水线", fontsize=8.5, color="#4f86c6")
servo(axR, HIP_F, lab="")
axR.text(HIP_F[0] + 12, HIP_F[1] + 4, "S1 髋滚转舵机（轴沿 x）", fontsize=8.5, color=INK, weight="bold")
s2R = FancyBboxPatch((0, 0), 10, 15, boxstyle="round,pad=1.2", fc=SERVO, ec=INK, lw=1.0, zorder=9); axR.add_patch(s2R)
s2Rt = axR.text(0, 0, "S2 扭转舵机\n（随腿摆，轴沿腿）", fontsize=8.5, color=INK, weight="bold")
axR.add_patch(Arc(HIP_F, 2 * L, 2 * L, theta1=250, theta2=290, color=GOLD, lw=1, ls=(0, (3, 3)), zorder=3))
axR.set_title("方案 F（推荐）鱼尾式 —— 后视图：S1 横向摆腿 = 沉浮", loc="left", fontsize=11, color=INK)
legR, = axR.plot([], [], color=INK, lw=6, solid_capstyle="round", zorder=5)
shaftR, = axR.plot([], [], color=GOLD, lw=1.6, zorder=6)
foilR = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=BLUE, zorder=7); axR.add_patch(foilR)
txtR = axR.text(60, -190, "", fontsize=9, color=INK2)

# ============ 方案 F 俯视图（x–y）：x 向前为右
axT.set_xlim(-190, 160); axT.set_ylim(-130, 130)
axT.add_patch(Ellipse((0, 0), LWL, BWL, fc=HULL, ec=INK2, lw=1, zorder=2, alpha=0.6))
axT.text(0, 0, "船体（俯视）", ha="center", va="center", fontsize=9, color=INK2)
axT.add_patch(FancyArrowPatch((120, 110), (150, 110), arrowstyle="-|>", mutation_scale=12, color="#4f86c6", lw=1.3)); axT.text(118, 118, "前进方向", fontsize=8.5, color="#4f86c6")
XLEG = -60.0                                       # 后腿 x
servo(axT, (XLEG, BWL / 2 - 2), lab="S1+S2（髋部）"); servo(axT, (XLEG, -BWL / 2 + 2), lab="")
axT.set_title("方案 F —— 俯视图：S2 绕腿轴转脚蹼 = 俯仰；脚蹼深度恒定", loc="left", fontsize=11, color=INK)
foilT, = axT.plot([], [], color=BLUE, lw=7, solid_capstyle="round", zorder=7); foilT2, = axT.plot([], [], color=BLUE, lw=7, solid_capstyle="round", zorder=7)
leT, = axT.plot([], [], "o", color=GOLD, ms=5, zorder=8); leT2, = axT.plot([], [], "o", color=GOLD, ms=5, zorder=8)
axT.text(-185, 118, "脚蹼竖放 → 俯视只见弦线（粗蓝线），金点 = 前缘", fontsize=8.5, color=INK2)
legT, = axT.plot([], [], "o", color=INK, ms=7, zorder=8); legT2, = axT.plot([], [], "o", color=INK, ms=7, zorder=8)
trailT, = axT.plot([], [], color=BLUE, lw=1.2, alpha=0.4, zorder=3); trailT2, = axT.plot([], [], color=BLUE, lw=1.2, alpha=0.4, zorder=3)
arrT = [None]; txtT = axT.text(-185, -122, "", fontsize=9, color=INK2)

# ============ 方案 D 侧视图（x–z）：x 向前为右
axD.set_xlim(-230, 150); axD.set_ylim(-200, 40)
axD.add_patch(Ellipse((0, WL), LWL, DRAFT * 2, fc=HULL, ec=INK2, lw=1, zorder=2))
axD.add_patch(FancyBboxPatch((-LWL / 2 + 10, WL), LWL - 20, 30, boxstyle="round,pad=2", fc="#e4dfd3", ec=INK2, lw=1, zorder=3))
axD.text(40, WL - 18, "船体（侧视）", ha="center", fontsize=8.5, color=INK2, zorder=4)
axD.axhline(WL, color="#4f86c6", lw=1.2, ls=(0, (5, 3))); axD.text(120, WL + 3, "水线", fontsize=8.5, color="#4f86c6")
HIP_D = np.array([-70.0, -75.0])                   # 龙骨支柱底部的髋
axD.add_patch(Polygon([[-78, -46], [-62, -46], [-66, -75], [-74, -75]], closed=True, fc="#888", ec=INK, lw=0.8, zorder=4)); axD.text(-58, -62, "固定流线支柱", fontsize=8, color=INK2)
servo(axD, (-70, -26), lab=""); axD.text(-56, -30, "S1 髋俯仰舵机\n（轴沿 y，经支柱内轴驱动）", fontsize=8.2, color=INK, weight="bold", va="center")
servo(axD, (-100, -26), lab=""); axD.text(-112, -26, "S2 俯仰舵机", fontsize=8.2, color=INK, weight="bold", ha="right", va="center")
axD.text(-225, -150, "金色 = 平行四边形连杆：\n把 S2 的角度原样传到脚蹼，\n与髋摆角无关", fontsize=8.5, color="#8a6d00", va="top")
axD.add_patch(Circle(HIP_D, 3, fc=GOLD, ec=INK, zorder=9))
axD.add_patch(Arc(HIP_D, 2 * L, 2 * L, theta1=180 + 15 - PHI0, theta2=180 + 15 + PHI0, color=GOLD, lw=1, ls=(0, (3, 3)), zorder=3))
axD.set_title("方案 D 海豚式 —— 侧视图：腿向后拖，S1 上下摆 = 沉浮", loc="left", fontsize=11, color=INK)
legD, = axD.plot([], [], color=INK, lw=6, solid_capstyle="round", zorder=5)
rodD, = axD.plot([], [], color=GOLD, lw=2, zorder=6)
crankD, = axD.plot([], [], color=GOLD, lw=2.5, zorder=6)
foilD = Polygon([[0, 0], [1, 0], [1, 1]], closed=True, color=BLUE, zorder=7); axD.add_patch(foilD)
trailD, = axD.plot([], [], color=BLUE, lw=1.2, alpha=0.4, zorder=3)
txtD = axD.text(-100, -192, "", fontsize=9, color=INK2)
PHI_MEAN_D = math.radians(15.0)                    # 腿平均下倾 15°

# ============ 指令曲线
lnP, = axS.plot([], [], color=INK, lw=2, label="S1 髋摆角 φ")
lnT, = axS.plot([], [], color=GOLD, lw=2, label="S2 俯仰角 θ（脚蹼绝对角）")
lnA, = axS.plot([], [], color=RED, lw=2, ls="--", label="攻角 α")
axS.legend(frameon=False, fontsize=9, ncol=3, labelcolor=INK2, loc="upper right")
ts, ps, th_s, al_s = [], [], [], []
trT, trT2, trD = [], [], []

def sec_poly(center, ang, c, flipx=False):
    """截面多边形：弦向 (cos ang, sin ang)，俯仰轴在 0.33c"""
    d = np.array([math.cos(ang), math.sin(ang)]); n = np.array([-d[1], d[0]])
    pts = [center + (0.33 * c - x * c) * d * (-1 if flipx else 1) * -1 + y * c * n for x, y in sec]
    return np.array(pts)

def planform_poly(center, ang):
    """俯视平面形：x 向后 → 前进方向为 +x 画面右，所以脚蹼弦向指向 -x（尾缘在左）。"""
    d = np.array([-math.cos(ang), -math.sin(ang)]); n = np.array([math.sin(ang), -math.cos(ang)])
    return np.array([center + (x - PIV) * d + y * n for x, y in pl])

def draw(fi):
    t = fi * DT
    ph, th, al = phi(t), theta(t), alpha(t)
    # --- F 后视图：腿从髋向下，横向摆 ph
    tip = HIP_F + L * np.array([math.sin(ph), -math.cos(ph)])
    legR.set_data([HIP_F[0], tip[0]], [HIP_F[1], tip[1]])
    shaftR.set_data([HIP_F[0], tip[0]], [HIP_F[1], tip[1]])
    dl = np.array([math.sin(ph), -math.cos(ph)]); s2c = HIP_F + 26 * dl
    s2R.set_bounds(s2c[0] - 5, s2c[1] - 7.5, 10, 15); s2Rt.set_position((s2c[0] + 14, s2c[1] - 2))
    # 脚蹼竖放：展向沿腿方向，看到的是弦的投影宽度 |cos θ|·c
    d = np.array([math.sin(ph), -math.cos(ph)]); n = np.array([math.cos(ph), math.sin(ph)])
    w = CM * abs(math.cos(th)) * 0.5 + 1.5
    foilR.set_xy([tip + B / 2 * d + w * n, tip + B / 2 * d - w * n, tip - B / 2 * d - w * n, tip - B / 2 * d + w * n])
    txtR.set_text(f"髋摆 φ = {math.degrees(ph):+5.1f}°　脚蹼横向位移 {L*math.sin(ph):+5.1f} mm")
    # --- F 俯视图：左右腿反相
    for k, (fp, le_, lt, trl, tr, sgn) in enumerate(((foilT, leT, legT, trailT, trT, +1), (foilT2, leT2, legT2, trailT2, trT2, -1))):
        y = sgn * (BWL / 2 + 8 + L * math.sin(ph))          # 左右镜像（反相）
        cen = np.array([XLEG, y]); ang = th * sgn
        dd = np.array([math.cos(ang), math.sin(ang)])         # 弦向：前缘朝 +x
        LE = cen + PIV * dd; TE = cen - (CM - PIV) * dd
        fp.set_data([LE[0], TE[0]], [LE[1], TE[1]]); le_.set_data([LE[0]], [LE[1]]); lt.set_data([cen[0]], [cen[1]])
        for i in range(len(tr)): tr[i][0] -= U * 1e3 * DT
        tr.append([cen[0] - 0.7 * CM, cen[1]])
        while tr and tr[0][0] < -185: tr.pop(0)
        trl.set_data([p[0] for p in tr], [p[1] for p in tr])
    txtT.set_text(f"S2 俯仰 θ = {math.degrees(th):+5.1f}°　攻角 α = {al:+5.1f}°（左右腿反相，侧向力抵消）")
    # --- D 侧视图：腿向后拖（-x），下倾 15° + 摆 ph
    a = PHI_MEAN_D + ph
    tipD = HIP_D + L * np.array([-math.cos(a), -math.sin(a)])
    legD.set_data([HIP_D[0], tipD[0]], [HIP_D[1], tipD[1]])
    # 平行四边形：平行杆偏置 9 mm，曲柄在髋和脚蹼各一段
    off = 9.0; nrm = np.array([math.sin(a), -math.cos(a)])
    crank_h = HIP_D + off * rot(np.array([1.0, 0.0]), th + math.pi / 2)
    crank_f = tipD + off * rot(np.array([1.0, 0.0]), th + math.pi / 2)
    rodD.set_data([crank_h[0], crank_f[0]], [crank_h[1], crank_f[1]])
    crankD.set_data([HIP_D[0], crank_h[0], np.nan, tipD[0], crank_f[0]], [HIP_D[1], crank_h[1], np.nan, tipD[1], crank_f[1]])
    # 脚蹼侧视 = 截面，弦向指向 -x（前缘朝前 +x）
    dd = np.array([math.cos(th), math.sin(th)]); nn = np.array([-dd[1], dd[0]])
    pts = np.array([tipD + (PIV - x * CM) * dd * -1 * -1 + 0 for x in [0]])  # placeholder
    poly = np.array([tipD + (x * CM - PIV) * (-dd) + y * CM * 1.6 * nn for x, y in sec])   # 厚度放大 1.6× 便于看清
    foilD.set_xy(poly)
    for i in range(len(trD)): trD[i][0] -= U * 1e3 * DT
    trD.append([tipD[0] - CM, tipD[1]])
    while trD and trD[0][0] < -225: trD.pop(0)
    trailD.set_data([p[0] for p in trD], [p[1] for p in trD])
    txtD.set_text(f"髋摆 φ = {math.degrees(ph):+5.1f}°　沉浮 {-(L*math.sin(a)-L*math.sin(PHI_MEAN_D)):+5.1f} mm　θ = {math.degrees(th):+5.1f}°　α = {al:+5.1f}°")
    # --- 曲线
    ts.append(t); ps.append(math.degrees(ph)); th_s.append(math.degrees(th)); al_s.append(al)
    lnP.set_data(ts, ps); lnT.set_data(ts, th_s); lnA.set_data(ts, al_s)
    return []

anim = animation.FuncAnimation(fig, draw, frames=NF, interval=1000 / FPS, blit=False)
out = os.path.join(HERE, "mechanism_anim.mp4")
anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=5000, extra_args=["-pix_fmt", "yuv420p", "-vf", "scale=1600:900"]))
print("→", out)
