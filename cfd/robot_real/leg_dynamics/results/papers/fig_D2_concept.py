#!/usr/bin/env python3
"""第二套设计（刚性板 + 脚踝羽化关节）概念草图 + 仿真方案示意 → fig_D2_concept.png"""
import os, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"; RED, BLUE, GOLD, GREEN = "#c2185b", "#2a78d6", "#d4a017", "#16a085"
D = np.genfromtxt(os.path.join(HERE, "optimal_trajectory.csv"), delimiter=",", names=True, dtype=None, encoding="utf-8"); EV = json.load(open(os.path.join(HERE, "optimal_trajectory_events.json")))
T = EV["T"]; tau = D["tau"]; X = D["x_abs_m"] * 1e3; Z = D["z_abs_m"] * 1e3; ph = D["phase"]
def style(ax):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=9)
fig = plt.figure(figsize=(17, 10), facecolor=SURF); gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1, 1], height_ratios=[1.35, 1], left=0.04, right=0.985, top=0.88, bottom=0.07, wspace=0.22, hspace=0.32)
# ---------- A 侧视：轨迹 + 腿 + 板的两种状态
ax = fig.add_subplot(gs[0, :2]); style(ax); ax.set_aspect("equal")
ax.axhspan(-130, -10, color="#dbe9f6", alpha=0.5, lw=0); ax.axhline(-10, color="#4f86c6", lw=1.2, ls=(0, (5, 3))); ax.text(-25, -7, "水线 z = −10 mm", color="#4f86c6", fontsize=9)
pw = ph == "power"; ax.plot(X[pw], Z[pw], color=RED, lw=3, zorder=3, label="划水：直线 84 mm，桨板正对（β = 0），梯形速度律")
ax.plot(X[~pw], Z[~pw], color=BLUE, lw=2, ls=(0, (4, 2)), zorder=3, label="回收：抬升 30 mm，桨板边缘迎流（β = 90°）")
HIP = np.array([63.4, 4.0]); ax.plot(*HIP, "o", color=INK, ms=9, zorder=6); ax.text(HIP[0] + 12, HIP[1] + 4, "髋关节 Hüfte（形式任意：转动 / 直线 / 2 自由度）", fontsize=9.5, color=INK, ha="left")
for k, colr in ((int(0.08 * len(D)), RED), (int(0.72 * len(D)), BLUE)):
    P = np.array([X[k], Z[k]]); Et = P + np.array([0, 71 if colr == RED else 45]); tip = P - np.array([0, 21])
    ax.plot([HIP[0], Et[0]], [HIP[1], Et[1]], color="#7b7b7b", lw=4, solid_capstyle="round", zorder=4)
    ax.plot(*Et, "o", color=GOLD, ms=8, mec=INK, zorder=7)
    ax.plot([Et[0], Et[0]], [Et[1], tip[1] + 42], color="#7b3f00", lw=3, zorder=5)
    if colr == RED: ax.add_patch(Rectangle((Et[0] - 2.8, tip[1]), 5.6, 42, color=RED, zorder=6))
    else: ax.add_patch(Rectangle((Et[0] - 22, tip[1]), 44, 42, color=BLUE, alpha=0.3, zorder=6)); ax.plot([Et[0] - 22, Et[0] + 22], [tip[1] + 21] * 2, color=BLUE, lw=1.2, zorder=7)
kR, kB = int(0.08 * len(D)), int(0.72 * len(D))
ax.annotate("脚踝 Knöchel = 羽化轴（沿桨杆）\n主动舵机，只有 0° / 90° 两个位置", (X[kR], Z[kR] + 71), xytext=(125, -20), fontsize=9.5, color=INK, arrowprops=dict(arrowstyle="->", color=INK2))
ax.annotate("划水：板面垂直于 x，侧视是 5.6 mm 的一条线", (X[kR] + 2.8, Z[kR] - 5), xytext=(125, -75), fontsize=9, color=RED, arrowprops=dict(arrowstyle="->", color=RED))
ax.annotate("回程：板面与 x–z 面共面\n对面内任何运动都是最小截面", (X[kB] - 22, Z[kB] - 21), xytext=(-28, -128), fontsize=9, color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
ax.set_xlim(-30, 200); ax.set_ylim(-135, 45); ax.set_xlabel("x [mm] → 前进方向", color=INK2); ax.set_ylabel("z [mm]", color=INK2)
ax.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK2); ax.set_title("A  第二套：刚性矩形板 + 脚踝羽化关节，走结构无关的目标轨迹（04：T 0.68 s，DUTY 0.42，v_max 0.35 m/s）", loc="left", fontsize=11, color=INK, pad=8)
# ---------- B 俯视：两种状态
ax = fig.add_subplot(gs[0, 2]); style(ax); ax.set_aspect("equal")
ax.add_patch(Rectangle((-2.8, -22), 5.6, 44, color=RED, zorder=4)); ax.text(0, 27, "划水 β = 0\n迎流面 44 × 42 = 18.5 cm²\n（= 120° 扇的面积）", ha="center", fontsize=9, color=RED)
ax.add_patch(Rectangle((-22, -70), 44, 5.6, color=BLUE, zorder=4)); ax.text(0, -80, "回程 β = 90°\n迎流 5.6 × 42 = 2.4 cm²（13 %）\n折扇收拢 12 mm 是 5.0 cm²（27 %）", ha="center", va="top", fontsize=9, color=BLUE)
ax.plot(0, 0, "o", color=GOLD, ms=9, mec=INK, zorder=6); ax.plot(0, -67, "o", color=GOLD, ms=9, mec=INK, zorder=6)
ax.add_patch(FancyArrowPatch((14, -8), (-8, -40), connectionstyle="arc3,rad=-0.5", arrowstyle="->", mutation_scale=16, color=INK2, lw=1.4)); ax.text(20, -30, "绕桨杆轴转 90°\n（舵机或被动）", fontsize=9, color=INK2)
ax.set_xlim(-45, 60); ax.set_ylim(-95, 45); ax.set_xlabel("x [mm]（前进方向 →）", color=INK2); ax.set_ylabel("y [mm]（横向）", color=INK2); ax.set_title("B  俯视：两个位置，中间无需停留", loc="left", fontsize=11, color=INK, pad=8)
ax.text(0, 8, "桨杆（羽化轴）", ha="center", fontsize=8, color=GOLD)
# ---------- C β(t) 与舵机速度
ax = fig.add_subplot(gs[1, :2]); style(ax)
tt = np.linspace(0, 1, 1000); fc, fo = EV["events_tau"]["close_cmd"], EV["events_tau"]["open_cmd"] - 1
def beta(tau_, tf):
    tfr = tf / T; dc = (tau_ - fc) % 1; do = (tau_ - fo) % 1; L = (fc - fo) % 1
    r = lambda x: 0.5 - 0.5 * np.cos(np.pi * np.clip(x, 0, 1))
    return np.where(do < L, 90 * (1 - r(do / tfr)), 90 * r(dc / tfr))
ax.axvspan(0, EV["DUTY"], color=RED, alpha=0.05, lw=0)
for tf, col, lab in ((0.04, RED, "40 ms → 2250 °/s（D2b）"), (0.10, GOLD, "100 ms → 900 °/s（D2c100）"), (0.20, BLUE, "200 ms → 450 °/s（D2c200，普通舵机）")):
    ax.plot(tt, beta(tt, tf), color=col, lw=2, label=f"羽化过渡 {lab}")
v = D["v_m_s"]; ax2 = ax.twinx(); ax2.plot(tau, v, color=MUTED, lw=1.2, ls=(0, (3, 2))); ax2.set_ylabel("桨头速度 [m/s]", color=MUTED); ax2.set_ylim(0, 0.9)
for s_ in ("top", "right"): ax2.spines[s_].set_visible(False)
for te, nm in ((0, "划水开始"), (fc, "收/羽化指令\n= 减速起点"), (EV["DUTY"], "划水结束"), (fo + 1, "展开指令\n= 回收减速起点")):
    ax.axvline(te, color=INK2, lw=0.8, ls=(0, (2, 2))); ax.text(te + 0.005, 95, nm, fontsize=8, color=INK2, va="top")
ax.set_xlim(0, 1); ax.set_ylim(0, 100); ax.set_xlabel("周期相位 τ", color=INK2); ax.set_ylabel("羽化角 β [°]", color=INK2); ax.legend(frameon=False, fontsize=9, loc="center right", labelcolor=INK2)
ax.set_title("C  D2‑2 扫描：舵机越慢，板在划水减速段 / 回程里越长时间是斜的 → 推力损失曲线 = 舵机选型依据", loc="left", fontsize=11, color=INK, pad=8)
# ---------- D 舵机要求估算
ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
txt = ("脚踝舵机要求（估算，板 44 × 42 × 5.6 mm，绕中轴）\n\n"
       "速度：90° / 40 ms = 2250 °/s（0.027 s/60°）\n          90° / 100 ms = 900 °/s（0.067 s/60°）\n          普通 9 g 舵机 ≈ 500 °/s → 180 ms\n\n"
       "力矩：划水时板平衡（压心过轴）≈ 0\n          旋转中水动力矩 ≈ 15 mN·m\n          + 附加质量惯量 ≈ 25 mN·m\n          → 峰值 ≈ 40 mN·m = 0.4 kg·cm（任何舵机都够）\n\n"
       "髋关节（沿用 04 的包络）：\n          v_max 0.35 m/s，斜坡 45 ms，峰值横向力 ≈ 0.5 N\n\n"
       "结论：瓶颈是脚踝的「速度」，不是力矩。\n若买不到 ≥ 1000 °/s 的微型舵机 → 被动风向标式羽化\n（自由铰 + 限位，减速段松开靠水流翻转）")
ax.text(0.0, 1.0, txt, va="top", ha="left", fontsize=9.5, color=INK, family=plt.rcParams["font.family"], linespacing=1.45, bbox=dict(fc="white", ec=GRID, boxstyle="round,pad=0.6"))
ax.set_title("D  驱动器要求", loc="left", fontsize=11, color=INK, pad=8)
fig.suptitle("第二套设计概念与仿真方案：刚性板 + 脚踝羽化 —— 六个算例回答「羽化值多少、板要多快转、板能做多大」", fontsize=13.5, x=0.04, ha="left", y=0.965, color=INK)
fig.text(0.04, 0.925, "CFD：桨板绕自身轴的旋转直接写进 6DoF 表（一次算例含整周期，不再两态拼接）；等面积板 44 mm 与第一套的 120° 扇可直接对比；D2w73 把板加宽到 73 mm 看面积换推力的代价", fontsize=9.5, color=MUTED)
fn = os.path.join(HERE, "fig_D2_concept.png"); plt.savefig(fn, dpi=110); print("→", fn)
