#!/usr/bin/env python3
"""双连杆 vs 三连杆：两种候选腿机构的对比图（结构无关，杆长为示意比例）
→ fig_linkage_2vs3.png / .pdf
"""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Arc

for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
RED, BLUE, GOLD, GREEN, WATER, PUR = "#c2185b", "#2a78d6", "#d4a017", "#16a085", "#eef5fb", "#7b3f9d"
LINK1, LINK2, LINK3 = "#5f6a72", "#8d7a5e", "#b8926a"
GHOST, GHOST2 = "#a8b0b6", "#c7bda8"

A1, A2      = 0.085, 0.075            # 双连杆
B1, B2, B3  = 0.075, 0.062, 0.032     # 三连杆
CHORD, THK  = 0.042, 0.0058

def rot(a): return np.array([np.cos(a), np.sin(a)])
def ik2(P, l1, l2, elbow=+1):
    r2 = P @ P; r = np.sqrt(r2)
    if r > l1+l2 or r < abs(l1-l2): return None
    c2 = np.clip((r2-l1**2-l2**2)/(2*l1*l2), -1, 1); q2 = elbow*np.arccos(c2)
    q1 = np.arctan2(P[1], P[0]) - np.arctan2(l2*np.sin(q2), l1+l2*np.cos(q2))
    return q1, q2
def fk2(q1, q2, l1, l2):
    J = l1*rot(q1); return J, J + l2*rot(q1+q2)
def ik3(P, th, l1, l2, l3, elbow=+1):
    s = ik2(P - l3*rot(th), l1, l2, elbow)
    return None if s is None else (s[0], s[1], th - s[0] - s[1])
def fk3(q1, q2, q3, l1, l2, l3):
    J1 = l1*rot(q1); J2 = J1 + l2*rot(q1+q2); return J1, J2, J2 + l3*rot(q1+q2+q3)
def plate(P, th, frac=0.25):
    d = rot(th); n = np.array([-d[1], d[0]])*THK/2
    LE = P + frac*CHORD*d; TE = P - (1-frac)*CHORD*d
    return np.array([LE+n, TE+n, TE-n, LE-n])
def arm(ax, pts, cols, lw=6.5, z=5, a=1.0):
    for i in range(len(pts)-1):
        ax.plot([pts[i][0], pts[i+1][0]], [pts[i][1], pts[i+1][1]], color=cols[i],
                lw=lw, solid_capstyle="round", zorder=z, alpha=a)
def jt(ax, pts, ms=9, z=8):
    for p in pts:
        ax.plot(*p, "o", color="#2b2b2b", ms=ms, zorder=z, mfc="white", mew=2.0)
def ang(ax, c, a0, a1, r, col, lab, fs=10):
    ax.add_patch(Arc(c, 2*r, 2*r, theta1=np.degrees(min(a0, a1)), theta2=np.degrees(max(a0, a1)),
                     color=col, lw=1.6, zorder=7))
    am = (a0+a1)/2
    ax.text(c[0]+r*1.5*np.cos(am), c[1]+r*1.5*np.sin(am), lab, color=col, fontsize=fs,
            ha="center", va="center", zorder=9, weight="bold")

# =======================================================================
fig = plt.figure(figsize=(15.5, 13.2), facecolor=SURF)
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.86], left=0.045, right=0.975,
                      top=0.885, bottom=0.028, wspace=0.10, hspace=0.26)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axC1 = fig.add_subplot(gs[1, 0]); axC2 = fig.add_subplot(gs[1, 1])
axD = fig.add_subplot(gs[2, :]); axD.axis("off")

def style(ax, xl, yl):
    ax.set_facecolor(WATER)
    for s in ("top", "right", "bottom", "left"): ax.spines[s].set_visible(False)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal")
style(axA, (-0.062, 0.168), (-0.205, 0.034))
style(axB, (-0.062, 0.168), (-0.205, 0.034))
for a in (axC1, axC2): style(a, (-0.062, 0.168), (-0.222, 0.034))

# ---------------------------------------------------------------- A 双连杆
q1a, q2a = np.deg2rad(-52), np.deg2rad(-56)
J, P = fk2(q1a, q2a, A1, A2)
arm(axA, [np.zeros(2), J, P], [LINK1, LINK2])
axA.add_patch(plt.Polygon(plate(P, q1a+q2a), closed=True, color=BLUE, zorder=6))
jt(axA, [np.zeros(2), J]); axA.plot(*P, "o", color=GOLD, ms=11, mec=INK, mew=1.6, zorder=9)
axA.plot([0, 0.060], [0, 0], color=MUTED, lw=1, ls=(0, (3, 2)), zorder=3)
axA.plot([J[0], J[0]+0.052], [J[1], J[1]], color=MUTED, lw=1, ls=(0, (3, 2)), zorder=3)
ang(axA, np.zeros(2), q1a, 0, 0.034, GREEN, "q1")
ang(axA, J, q1a+q2a, 0, 0.030, PUR, "q2")
axA.text(*(0.5*J+np.array([-0.023, 0.005])), "L1", color=LINK1, fontsize=12, weight="bold")
axA.text(*(0.5*(J+P)+np.array([0.008, 0.006])), "L2", color=LINK2, fontsize=12, weight="bold")
axA.text(0.004, 0.014, "髋 Hüfte", fontsize=10, color=INK)
axA.text(J[0]+0.008, J[1]+0.011, "膝 Knie", fontsize=10, color=INK)
axA.text(P[0]+0.012, P[1]+0.004, "桨板（刚性固定在 L2 上）", fontsize=10, color=BLUE, weight="bold")
axA.text(-0.057, -0.163, "自由度 = 2　→　只能独立控制位置 (x, z)\n"
                          "桨板姿态 θ = q1 + q2，被杆件位形锁死",
         fontsize=11.5, color=INK, va="top", linespacing=1.7)
axA.set_title("A  双连杆 Zweigelenk（2 自由度）", loc="left", fontsize=14, color=INK, pad=9)

# ---------------------------------------------------------------- B 三连杆
Pb = np.array([0.076, -0.108])
for th_g, al_ in ((np.deg2rad(-78), 0.30), (np.deg2rad(-38), 0.30)):
    s2 = ik3(Pb, th_g, B1, B2, B3, elbow=-1)
    if s2:
        g = fk3(*s2, B1, B2, B3)
        arm(axB, [np.zeros(2), *g], [GHOST, GHOST2, GHOST2], lw=5.0, z=2, a=al_)
        axB.add_patch(plt.Polygon(plate(g[2], th_g), closed=True, color=BLUE, alpha=0.22, zorder=2))
thb = np.deg2rad(-6)
q1b, q2b, q3b = ik3(Pb, thb, B1, B2, B3, elbow=-1)
J1, J2, Pb_ = fk3(q1b, q2b, q3b, B1, B2, B3)
arm(axB, [np.zeros(2), J1, J2, Pb_], [LINK1, LINK2, LINK3])
axB.add_patch(plt.Polygon(plate(Pb_, thb), closed=True, color=BLUE, zorder=6))
jt(axB, [np.zeros(2), J1, J2]); axB.plot(*Pb_, "o", color=GOLD, ms=11, mec=INK, mew=1.6, zorder=9)
axB.plot([0, 0.055], [0, 0], color=MUTED, lw=1, ls=(0, (3, 2)), zorder=3)
ang(axB, np.zeros(2), q1b, 0, 0.032, GREEN, "q1")
ang(axB, J1, q1b+q2b, q1b, 0.028, PUR, "q2")
ang(axB, J2, q1b+q2b+q3b, q1b+q2b, 0.020, GOLD, "q3")
axB.text(*(0.5*J1+np.array([-0.023, 0.004])), "L1", color=LINK1, fontsize=12, weight="bold")
axB.text(*(0.5*(J1+J2)+np.array([-0.018, -0.002])), "L2", color=LINK2, fontsize=12, weight="bold")
axB.text(*(0.5*(J2+Pb_)+np.array([0.004, 0.016])), "L3", color=LINK3, fontsize=12, weight="bold")
axB.text(0.004, 0.014, "髋", fontsize=10, color=INK)
axB.text(J1[0]-0.019, J1[1]+0.002, "膝", fontsize=10, color=INK)
axB.text(J2[0]-0.014, J2[1]-0.020, "踝 Knöchel", fontsize=10, color=INK, ha="right")
axB.text(-0.057, -0.163, "自由度 = 3　→　位置 (x, z) 与姿态 θ 可同时独立设定\n"
                          "淡色虚影：桨板位置不变，姿态换了三个值，都做得到",
         fontsize=11.5, color=INK, va="top", linespacing=1.7)
axB.set_title("B  三连杆 Dreigelenk（3 自由度）", loc="left", fontsize=14, color=INK, pad=9)

# ---------------------------------------------------------------- C 决定性差别
PC = np.array([0.110, -0.111])          # 选在双连杆的两组解都不等于任一需求姿态的位置
for ax, th_need, lab, col in ((axC1, np.deg2rad(-90), "需求①  划水：板正对来流", RED),
                              (axC2, np.deg2rad(0.0), "需求②  回程：板羽化侧立", GREEN)):
    for elbow, dyl in ((+1, 0.020), (-1, -0.022)):
        s_ = ik2(PC, A1, A2, elbow)
        if s_ is None: continue
        g1, g2 = s_; Jg, Pg = fk2(g1, g2, A1, A2)
        arm(ax, [np.zeros(2), Jg, Pg], [GHOST, GHOST2], lw=5.0, z=3)
        ax.add_patch(plt.Polygon(plate(Pg, g1+g2), closed=True, color="#8d969c", zorder=5))
        th_deg = (np.degrees(g1+g2)+180) % 360 - 180
        ax.text(Pg[0]-0.030, Pg[1]+dyl, f"θ = {th_deg:.0f}°", fontsize=10, color="#5f686e",
                ha="right", va="center")
    s3 = ik3(PC, th_need, B1, B2, B3, elbow=-1)
    g = fk3(*s3, B1, B2, B3)
    arm(ax, [np.zeros(2), g[0], g[1], g[2]], [LINK1, LINK2, LINK3], lw=5.5, z=6)
    ax.add_patch(plt.Polygon(plate(g[2], th_need), closed=True, color=col, zorder=8))
    jt(ax, [np.zeros(2), g[0], g[1]], ms=8, z=9)
    ax.plot(*PC, "o", color=GOLD, ms=13, mec=INK, mew=1.8, zorder=12)
    ax.text(PC[0]+0.016, PC[1]+0.030, "指定的桨板位置\n（两图相同）", fontsize=9.5, color=INK2, va="center")
    ax.text(-0.057, -0.172, f"三连杆：θ = {np.degrees(th_need):.0f}°　做得到",
            fontsize=11.5, color=col, weight="bold")
    ax.text(-0.057, -0.196, "双连杆：只有灰色那两个姿态，都不是要的",
            fontsize=11, color="#5f686e")
    ax.set_title(f"C{'①' if col == RED else '②'}  {lab}", loc="left", fontsize=13, color=INK, pad=8)

# ---------------------------------------------------------------- D 对比表
rows = [
    ("平面内自由度",        "2",                                "3"),
    ("可独立控制的量",      "位置 (x, z)",                      "位置 (x, z) + 姿态 θ"),
    ("桨板姿态",            "= q1 + q2，不可选",                "任意设定，与位置解耦"),
    ("羽化 / 攻角控制",     "做不到（除非另加被动铰）",         "做得到，且可连续调"),
    ("适配的推进方式",      "只能阻力型，回程靠减面积",         "阻力型 + 升力型都行"),
    ("每腿舵机 / 整机",     "2 / 8",                            "3 / 12"),
    ("逆运动学",            "闭式，2 组解（肘上 / 肘下）",      "闭式（先退桨板求腕点），2 组解"),
    ("自由度与任务",        "欠驱动：3 个量只控得住 2 个",      "恰定：3 个自由度控 3 个量"),
    ("水下代价",            "密封点少，末端惯量小",             "最远端多一个密封点，惯量与漏水风险最大"),
]
axD.set_xlim(0, 1); axD.set_ylim(0, 1)
X0, X1, X2 = 0.0, 0.30, 0.635
y0, dy = 0.905, 0.087
axD.add_patch(plt.Rectangle((X1-0.012, y0-len(rows)*dy+0.006), X2-X1-0.008, len(rows)*dy+0.052,
                            color="#f0efec", lw=0, zorder=0))
axD.add_patch(plt.Rectangle((X2-0.012, y0-len(rows)*dy+0.006), 1.012-X2, len(rows)*dy+0.052,
                            color="#e8f3ee", lw=0, zorder=0))
axD.text(X1, y0+0.020, "双连杆（2 DOF）", fontsize=12.5, color=INK, weight="bold")
axD.text(X2, y0+0.020, "三连杆（3 DOF）", fontsize=12.5, color=GREEN, weight="bold")
axD.plot([0, 1], [y0-0.002, y0-0.002], color=INK, lw=1.3)
for i, (k, a, b) in enumerate(rows):
    y = y0 - (i+0.62)*dy
    axD.text(X0, y, k, fontsize=11, color=INK2)
    axD.text(X1, y, a, fontsize=11, color=INK)
    axD.text(X2, y, b, fontsize=11, color=INK, weight="bold" if i in (1, 3, 4) else "normal")
    if i < len(rows)-1:
        axD.plot([0, 1], [y-0.030, y-0.030], color=GRID, lw=0.8)
axD.plot([0, 1], [y0-len(rows)*dy+0.010, y0-len(rows)*dy+0.010], color=INK, lw=1.3)
axD.text(X0, y0-len(rows)*dy-0.048,
         "折中方案：双连杆 + 一个被动铰（自由铰 + 限位），2 个舵机，靠水流自己把板翻过来。省一个水下舵机，"
         "但翻转时刻由流场决定、不可指定，需要实测标定。",
         fontsize=10.5, color=MUTED, va="top", linespacing=1.7)
axD.set_title("D  对比", loc="left", fontsize=14, color=INK, pad=6)

fig.text(0.045, 0.962, "腿部机构选型：双连杆 vs 三连杆", fontsize=19, color=INK, ha="left")
fig.text(0.045, 0.930,
    "结构待重新设计，这里只比机构本身，杆长为示意比例。核心判据只有一条：平面内桨板的位姿是 3 个量（x, z, θ）。"
    "双连杆 2 个自由度，位置到了姿态就被锁死；三连杆 3 个自由度，位置与姿态才能同时管住。"
    "羽化和攻角控制都是对 θ 的要求，所以这一条就决定了选型。",
    fontsize=11.5, color=INK2, ha="left", va="top", linespacing=1.7)

for ext in ("png", "pdf"):
    fn = os.path.join(HERE, f"fig_linkage_2vs3.{ext}")
    plt.savefig(fn, dpi=140 if ext == "png" else None, facecolor=SURF)
    print("→", fn)
