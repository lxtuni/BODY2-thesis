#!/usr/bin/env python3
"""被动俯仰尾鳍脚：完整设计图（侧视总图 / 脚部剖面 / 俯视 / 弹簧参数曲线 / 弹簧规格）→ fig_passive_foot_design.png/.pdf"""
import os, json, math, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Arc, FancyBboxPatch, Rectangle, Ellipse
import sys; sys.path.insert(0, os.path.join(HERE_, "..", "..")) if False else None
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
import leg_dynamics as ld
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
BLUE, RED, GOLD, GREEN, PUR, HULL, STEEL = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#c9c4b8", "#6b7280"

# ---------------- 设计参数
X_P = 9.0            # 铰轴距展中前缘 [mm]
X_AC = 18.6          # UVLM 识别的气动中心（展中前缘后）[mm]
C_MID, C_TIP, B = 28.7, 14.3, 86.0
L_STRUT = 105.0      # 支杆：E → 铰轴 [mm]
H0 = 30.0; F = 1.5
STOP = 45.0
sec = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_section_NACA0021.csv"), delimiter=",", skiprows=1)
pl = np.loadtxt(os.path.join(HERE, "fluke_AR4_t0.5_s35_planform_mm.csv"), delimiter=",", skiprows=1)
res = pd.DataFrame([r for r in json.load(open(os.path.join(HERE, "passive_pitch_sweep.json"))) if "error" not in r])
sel = res[(res.f == F) & (res.x_p == X_P)].sort_values("k")
c_pick = 0.15e-3
sub = sel[np.isclose(sel.c, c_pick)]
# 推荐 k：α_max 在 25–32° 且推力最大
K_PICK = 3e-3
best = sub.iloc[(sub.k - K_PICK).abs().argsort().iloc[0]]
b2 = res[(res.f == 2.0) & (res.x_p == X_P) & np.isclose(res.c, c_pick)]; b2 = b2.iloc[(b2.k - K_PICK).abs().argsort().iloc[0]]

fig = plt.figure(figsize=(18, 12.5), facecolor=SURF)
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], height_ratios=[1.25, 1], left=0.035, right=0.985, top=0.845, bottom=0.05, hspace=0.32, wspace=0.16)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, 0]); axE = fig.add_subplot(gs[1, 1]); axF = fig.add_subplot(gs[1, 2])
for ax in (axA, axB, axC, axF):
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
for ax in (axD, axE):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(color=GRID, lw=0.8); ax.tick_params(colors=INK2, labelsize=9)

def dim(ax, p0, p1, text, off=0, fs=8, color=INK2, rot=0):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0; n = np.array([-d[1], d[0]]) / max(np.linalg.norm(d), 1e-9) * off
    ax.annotate("", p1 + n, p0 + n, arrowprops=dict(arrowstyle="<->", color=color, lw=0.9, shrinkA=0, shrinkB=0))
    m = (p0 + p1) / 2 + n; ax.text(m[0], m[1], text, fontsize=fs, color=color, ha="center", va="bottom", rotation=rot)

# ============ A 侧视总图（游泳姿态，ψ = 90°）
axA.set_facecolor(WATER)
O2 = ld.O2 * 1e3
psi = math.radians(90)
hull_y = -10
axA.add_patch(Ellipse((O2[0] - 60, hull_y), 233, 76, fc=HULL, ec=INK2, lw=1, zorder=2))
axA.add_patch(FancyBboxPatch((O2[0] - 60 - 105, hull_y), 210, 26, boxstyle="round,pad=2", fc="#e4dfd3", ec=INK2, lw=1, zorder=3))
axA.axhline(hull_y, color="#4f86c6", lw=1.2, ls=(0, (5, 3))); axA.text(O2[0] - 165, hull_y + 3, "水线 z = −10", fontsize=8, color="#4f86c6")
axA.plot([O2[0]], [O2[1]], "o", color=GOLD, ms=8, mec=INK, zorder=9); axA.text(O2[0] + 5, O2[1] + 6, "O2 髋轴（ψ 摇臂，游泳时固定 90°）", fontsize=8, color=INK, weight="bold")
for ph, a in ((40, 0.25), (90, 1.0), (140, 0.25)):
    p = ld.pose(psi, math.radians(ph)); C, D, E = p["C"] * 1e3, p["D"] * 1e3, p["E"] * 1e3
    axA.plot([O2[0], C[0]], [O2[1], C[1]], color=INK, lw=4, alpha=a, zorder=5)
    axA.plot([O2[0], D[0], E[0], C[0]], [O2[1], D[1], E[1], C[1]], color=GOLD, lw=2, alpha=a, zorder=6)
    foot = E + np.array([0, -L_STRUT])
    axA.plot([E[0], foot[0]], [E[1], foot[1]], color=INK, lw=3.5, alpha=a, zorder=5)
    th = math.radians(best.th0) * (1 if ph == 40 else (-1 if ph == 140 else 0))   # 端点姿态示意
    d = np.array([math.cos(th), math.sin(th)]); n = np.array([-d[1], d[0]])
    poly = np.array([foot + (X_P - x * C_MID) * d + y * C_MID * 1.4 * n for x, y in sec])
    axA.add_patch(Polygon(poly, closed=True, color=BLUE, alpha=0.35 + 0.6 * a, zorder=7))
    axA.plot([foot[0]], [foot[1]], "o", color=RED, ms=5, zorder=9)
# E 的半圆轨迹
th2 = np.radians(np.linspace(102, 258, 80)); Cc = (ld.O2 + ld.L_ROCK_C * ld.rot(psi)) * 1e3
Ep = Cc[None, :] + 39 * np.column_stack([np.cos(th2), np.sin(th2)])
axA.plot(Ep[:, 0], Ep[:, 1], color=GOLD, lw=1, ls=(0, (3, 3)))
axA.plot(Ep[:, 0], Ep[:, 1] - L_STRUT, color=RED, lw=1, ls=(0, (3, 3)))
pE = ld.pose(psi, math.radians(90))["E"] * 1e3
dim(axA, (pE[0] + 55, pE[1] - L_STRUT - H0), (pE[0] + 55, pE[1] - L_STRUT + H0), f"沉浮 ±{H0:.0f}\n(φ 40°→140°)", off=0, fs=8, color=RED)
dim(axA, (pE[0] - 22, pE[1]), (pE[0] - 22, pE[1] - L_STRUT), f"支杆 {L_STRUT:.0f}", off=0, fs=8, rot=90)
axA.text(pE[0] - 20, pE[1] - L_STRUT - H0 - 22, f"铰轴中心深度 {-(pE[1]-L_STRUT)+hull_y:.0f} mm，顶端离水面 {-(pE[1]-L_STRUT+H0)+hull_y:.0f} mm ≈ {(-(pE[1]-L_STRUT+H0)+hull_y)/C_MID:.1f} 弦", fontsize=8, color=INK2)
axA.text(Cc[0] + 44, Cc[1] + 10, "φ 平行四边形：E 走 R39 半圆\n= 竖直平移，不转动", fontsize=8, color="#8a6d00")
axA.set_xlim(O2[0] - 190, O2[0] + 130); axA.set_ylim(pE[1] - L_STRUT - H0 - 40, 40)
axA.set_title("A  侧视总图（游泳姿态）：摇臂固定竖直，平行四边形做沉浮，脚部被动俯仰", loc="left", fontsize=11.5, color=INK)

# ============ B 脚部剖面（铰 + 扭簧 + 限位）
axB.set_facecolor(SURF)
sc = 3.0   # 放大倍数（画布单位 = mm/… 只是比例）
foil = np.array([((x * C_MID - X_P), y * C_MID) for x, y in sec]) * sc
axB.add_patch(Polygon(-foil * np.array([1, -1]) * np.array([1, 1]), closed=True, fc="#9fc3e6", ec=INK, lw=1.2, zorder=4))  # 前缘朝左
# 上面一行把 x 取反：前缘在左（来流从左），铰轴在原点
axB.add_patch(Circle((0, 0), 2.2 * sc, fc=GOLD, ec=INK, lw=1.2, zorder=8)); axB.text(4 * sc, 3.5 * sc, "铰轴 Ø2 mm 钢针", ha="left", fontsize=8, color=INK)
axB.plot([-(X_AC - X_P) * sc], [0], "x", color=RED, ms=9, mew=2, zorder=9); axB.text(-(X_AC - X_P) * sc - 4 * sc, 4.5 * sc, f"× 气动中心 AC（UVLM 识别，前缘后 {X_AC:.1f} mm）", ha="right", va="bottom", fontsize=7.8, color=RED)
# 支杆
axB.add_patch(Rectangle((-1.5 * sc, 2 * sc), 3 * sc, 30 * sc, fc=INK, ec=INK, zorder=3)); axB.text(2.5 * sc, 20 * sc, "支杆 3×5.6 → 改流线 4×8\n（沿自身轴运动，阻力≈0）", fontsize=8, color=INK2, va="center")
# 扭簧示意
ang = np.linspace(0, 8 * math.pi, 400); axB.plot(-7.5 * sc + 1.2 * sc * np.cos(ang), 2.2 * sc + ang / (8 * math.pi) * 6 * sc + 0.6 * sc * np.sin(ang), color=STEEL, lw=1.5, zorder=7)
axB.text(-9.5 * sc, 9.5 * sc, "扭簧 k（琴钢丝，\n或 TPU 柔性铰替代）", ha="right", fontsize=8, color=STEEL, va="center")
# 限位
for s_ in (1, -1):
    a = math.radians(STOP) * s_
    axB.plot([0, 12 * sc * math.cos(a)], [0, 12 * sc * math.sin(a)], color=RED, lw=1, ls=(0, (3, 2)), zorder=6)
axB.add_patch(Arc((0, 0), 20 * sc, 20 * sc, theta1=-STOP, theta2=STOP, color=RED, lw=1.2, zorder=6))
axB.text(13 * sc, 0, f"限位 ±{STOP:.0f}°\n（α 削成平台型）", fontsize=8, color=RED, va="center")
dim(axB, (X_P * sc, -9 * sc), (-(C_MID - X_P) * sc, -9 * sc), f"展中弦 {C_MID:.1f}", off=0, fs=8)
dim(axB, (X_P * sc, -13.5 * sc), (0, -13.5 * sc), f"x_p = {X_P:.0f} mm（0.31 c）", off=0, fs=8, color=GOLD)
axB.add_patch(FancyArrowPatch((28 * sc, 14 * sc), (18 * sc, 14 * sc), arrowstyle="-|>", mutation_scale=12, color="#4f86c6", lw=1.3)); axB.text(28.5 * sc, 14 * sc, "来流", fontsize=8.5, color="#4f86c6", va="center")
axB.text(-31 * sc, 31 * sc, "阻尼：铰轴涂硅脂或小旋转阻尼器，c ≥ 0.1 mN·m·s/rad\n（x_p < 12 mm 时流体阻尼为负，无阻尼会自激颤振）", fontsize=8, color=INK2, va="top")
axB.set_xlim(-32 * sc, 36 * sc); axB.set_ylim(-17 * sc, 34 * sc)
axB.set_title("B  脚部剖面（展中）：铰轴在气动中心之前 → 风向标稳定", loc="left", fontsize=11.5, color=INK)

# ============ C 俯视
axC.set_facecolor(SURF)
axC.add_patch(Polygon(np.column_stack([pl[:, 1], -(pl[:, 0] - X_P)]), closed=True, fc="#9fc3e6", ec=INK, lw=1.2, zorder=3))
axC.axhline(0, color=GOLD, lw=1.6, ls=(0, (4, 2)), zorder=5); axC.text(-B / 2 - 3, 0, "铰轴", fontsize=8.5, color=GOLD, va="center", ha="right")
axC.add_patch(Rectangle((-2, -4), 4, 8, fc=INK, zorder=6)); axC.text(6, 5, "支杆（展中，Ø4×8 流线）", ha="left", fontsize=8, color=INK)
axC.plot([0], [-(X_AC - X_P)], "x", color=RED, ms=8, mew=2, zorder=7)
axC.add_patch(FancyArrowPatch((0, 28), (0, 18), arrowstyle="-|>", mutation_scale=12, color="#4f86c6", lw=1.3)); axC.text(4, 24, "来流", fontsize=8.5, color="#4f86c6")
dim(axC, (-B / 2, -30), (B / 2, -30), f"展 {B:.0f}", off=0, fs=8)
dim(axC, (B / 2 + 10, X_P), (B / 2 + 10, X_P - C_MID), "", off=0, fs=8); axC.text(B / 2 + 13, X_P - C_MID / 2, "展中弦\n28.7", fontsize=8, color=INK2, va="center")
axC.text(-B / 2, -40, "月牙形 AR 4 · 后掠 35° · 翼尖弦 14.3 · S 18.5 cm²\n截面 NACA 0021（t_max 6.0 mm @ 8.6 mm，前缘 R 1.4，尾缘尖）\n面积不变、形状换新；打印件 ≈ 5 g", fontsize=8, color=INK2, va="top")
axC.set_xlim(-62, 72); axC.set_ylim(-62, 34)
axC.set_title("C  俯视：脚蹼平面形，铰轴沿展向穿过展中", loc="left", fontsize=11.5, color=INK)

# ============ D 弹簧刚度扫描：θ0、ψ、α_max
for c_, ls, lab in ((0.05e-3, ":", "c 0.05"), (0.15e-3, "-", "c 0.15"), (0.3e-3, "--", "c 0.30")):
    s_ = sel[np.isclose(sel.c, c_)]
    if not len(s_): continue
    axD.plot(s_.k * 1e3, s_.th0, ls, color=PUR, lw=1.8, marker="o", ms=4, label=f"θ0 [°]（{lab} mN·m·s/rad）")
    axD.plot(s_.k * 1e3, s_.psi, ls, color=GOLD, lw=1.8, marker="s", ms=4, label=f"相位 ψ [°]" if c_ == c_pick else None)
    axD.plot(s_.k * 1e3, s_.alpha_max, ls, color=RED, lw=1.8, marker="^", ms=4, label=f"α_max [°]" if c_ == c_pick else None)
axD.axhspan(15, 30, color="#fde8ee", zorder=0); axD.text(sel.k.max() * 1e3 * 0.98, 16, "α_max 目标区 15–30°", fontsize=8, color=RED, ha="right")
axD.axhline(90, color=GOLD, lw=0.8, ls=(0, (2, 2))); axD.text(sel.k.min() * 1e3, 92, "理想相位 90°", fontsize=8, color=GOLD)
axD.axvline(K_PICK * 1e3, color=INK, lw=1, ls=(0, (4, 3))); axD.text(K_PICK * 1e3 * 1.04, 60, f"选 k = {K_PICK*1e3:.0f} mN·m/rad", fontsize=8.5, color=INK, weight="bold")
axD.set_xscale("log"); axD.set_xlabel("扭簧刚度 k [mN·m/rad]", color=INK2); axD.set_ylabel("角度 [°]", color=INK2)
axD.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper right", ncol=1)
axD.set_title(f"D  被动响应 vs 弹簧刚度（f {F:g} Hz, h0 {H0:.0f} mm, x_p {X_P:.0f} mm）", loc="left", fontsize=11.5, color=INK)

# ============ E 推力 / 功率 / 效率
for c_, ls, lab in ((0.05e-3, ":", "0.05"), (0.15e-3, "-", "0.15"), (0.3e-3, "--", "0.30")):
    s_ = sel[np.isclose(sel.c, c_)]
    if not len(s_): continue
    axE.plot(s_.k * 1e3, s_["T"] * 1e3, ls, color=BLUE, lw=2, marker="o", ms=4, label=f"推力 T [mN]（c {lab}）")
    axE.plot(s_.k * 1e3, s_.P * 1e3, ls, color=GREEN, lw=1.6, marker="s", ms=4, label="沉浮舵机功率 P [mW]" if c_ == c_pick else None)
axE.axvline(K_PICK * 1e3, color=INK, lw=1, ls=(0, (4, 3)))
axE.text(K_PICK * 1e3 * 1.04, best["T"] * 1e3 * 0.55, f"T = {best['T']*1e3:.0f} mN\nP = {best.P*1e3:.0f} mW\nη = {min(best.eta, 0.99):.2f}（无黏，偏高）", fontsize=8.5, color=INK, va="top")
axE.set_xscale("log"); axE.set_xlabel("扭簧刚度 k [mN·m/rad]", color=INK2); axE.set_ylabel("mN / mW", color=INK2)
axE.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper left")
axE.set_title("E  推力与功率 vs 刚度（UVLM 耦合，含 C_D 0.03）", loc="left", fontsize=11.5, color=INK)

# ============ F 弹簧规格表
E_MW = 200e9
rows = []
for D_, N_ in ((5.0, 3), (6.0, 4), (8.0, 5)):
    d = (K_PICK * 64 * D_ * 1e-3 * N_ / E_MW) ** 0.25 * 1e3
    rows.append(("琴钢丝扭簧", f"d {d:.2f} mm · D {D_:.0f} mm · {N_} 圈（双向：两只对置或预压中立）"))
E_TPU, E_PETG = 30e6, 2.0e9
for lab, E_, t_, L_ in (("TPU 95A 柔性铰", E_TPU, 1.2, 8.0), ("TPU 95A 柔性铰", E_TPU, 1.0, 8.0), ("PETG 柔性铰", E_PETG, 0.35, 8.0)):
    k_ = E_ * 20e-3 * (t_ * 1e-3) ** 3 / (12 * L_ * 1e-3)
    rows.append((lab, f"宽 20 × 厚 {t_} × 长 {L_:.0f} mm → k ≈ {k_*1e3:.1f} mN·m/rad"))
d_ss = (K_PICK * 12 * 8e-3 / (200e9 * 20e-3)) ** (1 / 3) * 1e3
rows.append(("弹簧钢片 0.05 级", f"宽 20 × 长 8 × 厚 {d_ss:.3f} mm（薄，易疲劳，不推荐）"))
axF.set_xlim(0, 1); axF.set_ylim(0, 1)
axF.add_patch(FancyBboxPatch((0.01, 0.02), 0.98, 0.96, boxstyle="round,pad=0.01", fc="#f5f4f1", ec=GRID, lw=1))
axF.text(0.04, 0.95, f"F  弹簧与阻尼规格（目标 k = {K_PICK*1e3:.0f} mN·m/rad = {K_PICK*1e3*math.pi/180:.3f} mN·m/°）", fontsize=11, color=INK, weight="bold", va="top")
y = 0.85
for lab, txt in rows:
    axF.text(0.04, y, f"• {lab}：{txt}", fontsize=8.0, color=INK2, va="top"); y -= 0.07
y -= 0.02
axF.text(0.04, y, f"• 阻尼 c ≥ 0.10–0.15 mN·m·s/rad：流体阻尼 C_f ≈ −0.05…−0.09（x_p 6–9 mm，\n   负值 = 无阻尼会自激）；铰轴硅脂 + 0.1 mm 轴向间隙，实测调到不颤为止\n"
                     f"• 限位 ±45°；中立位置 θ = 0（弦 ⟂ 支杆）；预紧 0\n"
                     f"• 铰轴 Ø2 mm 不锈钢针，PTFE/POM 衬套；I ≈ {best.I*1e6:.2f}e-6 kg·m²，\n   含附加质量固有频率 ≈ 4–5 Hz > 2 Hz（共振以下）\n"
                     f"• 流体刚度 K_f ≈ {best.Kf*1e3:.1f} mN·m/rad（1.5 Hz），χ = K_f/(k+K_f) ≈ {best.chi:.2f}",
          fontsize=8.0, color=INK2, va="top", linespacing=1.55)

fig.text(0.035, 0.965, "BODY2 · 被动俯仰尾鳍脚 —— 在现有髋 ψ + 平行四边形 φ 上只换脚：完整设计图与弹簧参数", fontsize=16.5, color=INK, va="top")
fig.text(0.035, 0.925,
    f"运动：ψ 固定 90°（摇臂竖直）；φ = 90° ± asin(h/39) 做正弦沉浮 h0 = {H0:.0f} mm（φ 峰值角速度 1.5 Hz 415 °/s、2 Hz 554 °/s；7465W 余量 2.0× / 1.5×）。"
    f"俯仰完全被动：扭簧 + 限位 + 阻尼，无第三个舵机。\n"
    f"被动响应（UVLM 流固耦合，基波相量法，x_p {X_P:.0f} mm，k = {K_PICK*1e3:.0f} mN·m/rad，c 0.15）：巡航 1.5 Hz → θ0 ±{best.th0:.0f}°、相位 {best.psi:.0f}°、α_max {best.alpha_max:.0f}°、推力 {best['T']*1e3:.0f} mN/桨、沉浮功率 {best.P*1e3:.0f} mW；"
    f"冲刺 2 Hz → θ0 ±{b2.th0:.0f}°、α_max {b2.alpha_max:.0f}°（会部分失速）、推力 ≤ {b2['T']*1e3:.0f} mN、{b2.P*1e3:.0f} mW。\n"
    f"UVLM 无黏、无失速：绝对值上偏 ~20%，趋势可信。连杆尺寸按 leg_dynamics.py 的 CAD 值（平行四边形 39 mm）；照片上的腿更长，规律不变、数值按实际杆长代入。",
    fontsize=9.4, color=MUTED, va="top", linespacing=1.65)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"fig_passive_foot_design.{ext}"), dpi=160 if ext == "png" else None, facecolor=SURF)
json.dump(dict(best=best.to_dict(), K_PICK=K_PICK, X_P=X_P, X_AC=X_AC, L_STRUT=L_STRUT, H0=H0, F=F, STOP=STOP), open(os.path.join(HERE, "passive_foot_design.json"), "w"), indent=1, default=float)
print("best", best[["k", "c", "th0", "psi", "chi", "alpha_max", "T", "P", "eta", "resid", "Kf", "Cf"]].to_dict())
