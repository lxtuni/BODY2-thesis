#!/usr/bin/env python3
"""
阻力式划水的“最优行程”准定常理论（结构无关）—— 四张图
    python3 optimal_stroke_theory.py [--out .]
A  折叠比 a = A_min/A_max 决定的最优占空比 DUTY*（划水速率给定、同弧长、U = 0）：x* = v_r/v_p = sqrt(1+1/a) − 1，DUTY* = x*/(1+x*)
B  周期平均推力随回收/划水速度比 x 的变化（几个 a），标出最优
C  同峰值速率下，余弦律 vs 梯形律（匀速划水）的速度曲线与有效推力窗口（v > U 才产生推力）
D  结构无关的“理论最优单腿行程”时间线：速度 / 面积 / 深度的开合时序
输出 fig_optimal_stroke_theory.png
"""
import os, argparse, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__))); a = ap.parse_args()
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
BLUE, ORANGE, PURPLE, RED = "#2a78d6", "#eb6834", "#7a5cc0", "#e34948"
RAMP = ["#8fb9ea", "#2a78d6", "#16416f"]

def duty_opt(a): x = np.sqrt(1 + 1 / a) - 1; return x / (1 + x), x
fig = plt.figure(figsize=(15.4, 10.4), facecolor=SURF)
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.27, left=0.06, right=0.985, top=0.875, bottom=0.06)
def style(ax):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9.5)
def title2(ax, t1, t2): ax.set_title(t1 + "\n" + t2, loc="left", fontsize=11, color=INK, pad=9)

# ---------------- A：DUTY*(a) ----------------
ax = fig.add_subplot(gs[0, 0]); style(ax)
aa = np.linspace(0.05, 1.0, 400); d, x = duty_opt(aa)
ax.plot(aa, d, color=INK, lw=2.2, zorder=3, label="无速率限制：DUTY* = x*/(1+x*)")
ax.plot(aa, np.minimum(d, 0.5), color=BLUE, lw=2.2, ls=(0, (5, 3)), zorder=3, label="回收速率 ≤ 划水速率时：min(DUTY*, 0.5)")
pts = [(0.445, "麝鼠 Fish 1984\n实测 DUTY 0.45", ORANGE, (-12, -26), "right", ""), (0.67, "BODY2 几何比\n（收拢 30 mm）", BLUE, (8, 14), "left", ""),
       (0.22, "BODY2 CFD 有效比\n（含蜷缩抬起）", BLUE, (12, 6), "left", "，速率受限 → 0.50")]
for av, nm, c, off, ha, extra in pts:
    dv, _ = duty_opt(av); ax.scatter([av], [dv], s=95, color=c, zorder=5, edgecolors=SURF, linewidths=1.6)
    ax.annotate(f"{nm}\n→ DUTY* {dv:.2f}{extra}", (av, dv), xytext=off, textcoords="offset points", fontsize=9, color=c, ha=ha, va="top" if off[1] < 0 else "bottom")
ax.scatter([0.445], [0.45], s=60, marker="x", color=ORANGE, zorder=6, linewidths=2)
ax.set_xlim(0, 1.0); ax.set_ylim(0.2, 0.8); ax.set_xlabel("回收相 / 划水相迎流面积比 a = A_min / A_max", color=INK2, fontsize=10)
ax.set_ylabel("周期平均推力最大的占空比 DUTY*", color=INK2, fontsize=10); ax.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK2)
title2(ax, "A  最优占空比由折叠比决定（准定常，U = 0，划水速率给定）", "     麝鼠 a = 0.445 → 理论 0.445，实测 0.45；折叠越好，回收越该快")

# ---------------- B：T̄ vs x ----------------
ax = fig.add_subplot(gs[0, 1]); style(ax)
xx = np.linspace(0.05, 3.0, 500)
for k, (av, lab) in enumerate([(0.67, "a = 0.67（BODY2 几何比）"), (0.445, "a = 0.45（麝鼠）"), (0.22, "a = 0.22（BODY2 CFD 有效比）")]):
    f = (1 - av * xx) * xx / (1 + xx); f = f / f.max() if False else f
    ax.plot(xx, f, color=RAMP[2 - k], lw=2.2, zorder=3, label=lab)
    dv, xs = duty_opt(av); fs = (1 - av * xs) * xs / (1 + xs)
    ax.scatter([xs], [fs], s=70, color=RAMP[2 - k], zorder=5, edgecolors=SURF, linewidths=1.5)
    off, ha, va = [((0, -12), "center", "top"), ((44, 2), "left", "bottom"), ((8, 6), "left", "bottom")][k]
    ax.annotate(f"x* = {xs:.2f}\nDUTY* = {dv:.2f}", (xs, fs), xytext=off, textcoords="offset points", fontsize=8.8, color=RAMP[2 - k], ha=ha, va=va)
ax.axvline(1.0, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2); ax.text(1.02, 0.02, "x = 1：回收与划水同速\n（舵机速率上限时的可达边界）", fontsize=8.5, color=INK2, va="bottom")
ax.set_xlim(0, 3.0); ax.set_ylim(0, 0.55); ax.set_xlabel("回收速度 / 划水速度  x = v_r / v_p", color=INK2, fontsize=10)
ax.set_ylabel("周期平均推力 ∝ (1 − a·x)·x / (1 + x)", color=INK2, fontsize=10); ax.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK2)
title2(ax, "B  周期平均推力随回收速度的变化", "     折叠越差（a 大），回收越该慢；BODY2 的有效 a ≈ 0.22 → 回收尽量快")

# ---------------- C：余弦 vs 梯形 ----------------
ax = fig.add_subplot(gs[1, 0]); style(ax)
tt = np.linspace(0, 1, 600); vc = np.sin(np.pi * tt); r = 0.12
vt = np.clip(np.minimum(tt / r, (1 - tt) / r), 0, 1)
U = 0.3
ax.fill_between(tt, 0, U, color=RED, alpha=0.08, lw=0, zorder=1); ax.axhline(U, color=RED, lw=1, ls=(0, (4, 3)), zorder=2)
ax.text(0.5, U - 0.03, "v < U：桨张开也只产生阻力（U / v_max = 0.3 示例）", ha="center", va="top", fontsize=8.8, color=RED)
ax.plot(tt, vc, color=ORANGE, lw=2.3, zorder=3, label="余弦律（现有全部算例）")
ax.plot(tt, vt, color=BLUE, lw=2.3, zorder=3, label="梯形律：12 % 加减速 + 匀速")
Ic, It = np.trapezoid(vc ** 2, tt), np.trapezoid(vt ** 2, tt); Ec, Et = np.trapezoid(vc ** 3, tt), np.trapezoid(vt ** 3, tt)
uc = np.trapezoid(np.clip(vc - U, 0, None) ** 2, tt) - np.trapezoid(np.clip(U - vc, 0, None) ** 2 * (vc < U), tt)
ut = np.trapezoid(np.clip(vt - U, 0, None) ** 2, tt) - np.trapezoid(np.clip(U - vt, 0, None) ** 2 * (vt < U), tt)
ax.text(0.02, 1.16, f"同峰值速率、同划水时长：冲量 ∫v²dt 梯形/余弦 = {It / Ic:.2f}，能量 ∫v³dt = {Et / Ec:.2f}；\n"
                    f"U/v_max = 0.3 时净推进冲量 ∫(v−U)²·sgn dt 之比 = {ut / uc:.2f}；张开窗口：余弦 {1 - 2 * np.arcsin(U) / np.pi:.0%} vs 梯形 {1 - 2 * r * U:.0%} 的划水相",
        fontsize=8.8, color=INK2, va="top")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.2); ax.set_xlabel("划水相内的相位 t / τ_p", color=INK2, fontsize=10); ax.set_ylabel("桨速 / 峰值速率 v / v_max", color=INK2, fontsize=10)
ax.legend(frameon=False, fontsize=9, loc="lower center", labelcolor=INK2)
title2(ax, "C  同一舵机速率上限下，匀速划水的冲量≈余弦律的两倍", "     余弦律把速率预算浪费在起止段；U > 0 时起止段还在产生阻力")

# ---------------- D：理论最优行程时间线 ----------------
ax = fig.add_subplot(gs[1, 1]); style(ax); ax.grid(False)
for s in ("left",): ax.spines[s].set_visible(False)
tau = np.linspace(0, 1, 1000); D = 0.5; rp = 0.06; pause = 0.04
def trap(t0, t1, ramp, sign):
    seg = (tau >= t0) & (tau < t1); u = (tau - t0) / (t1 - t0); v = np.clip(np.minimum(u / ramp, (1 - u) / ramp), 0, 1); return np.where(seg, sign * v, 0.0)
v = trap(0, D, rp / D, 1.0) + trap(D + pause, 1.0, rp / (1 - D - pause), -1.0)
area = np.where((tau >= 0.0) & (tau < D), 1.0, 0.3); area = np.where((tau >= D) & (tau < D + pause), 0.3 + 0.7 * (1 - (tau - D) / pause), area)
depth = np.where(tau < D, 1.0, 0.0); depth = np.where((tau >= D) & (tau < D + 0.12), 1 - (tau - D) / 0.12, depth); depth = np.where(tau > 0.88, (tau - 0.88) / 0.12, depth)
base = {"v": 2.6, "A": 1.4, "z": 0.2}; H = 0.9
for key, y0, sig, col, lab in (("v", base["v"], v, BLUE, "桨速 v（+ 向后划）"), ("A", base["A"], area, ORANGE, "迎流面积 A / A_max"), ("z", base["z"], depth, PURPLE, "入水深度\n1 = 最深，0 = 抬出")):
    ax.plot(tau, y0 + H * (sig if key != "v" else 0.5 + 0.5 * sig), color=col, lw=2.2, zorder=3)
    ax.axhline(y0 + (H / 2 if key == "v" else 0), color=GRID, lw=0.8, zorder=1)
    ax.text(-0.01, y0 + H * (0.5 if key == "v" else 0.5), lab, ha="right", va="center", fontsize=9, color=INK2)
ax.axvspan(0, D, color=BLUE, alpha=0.06, lw=0); ax.axvspan(D, D + pause, color=RED, alpha=0.12, lw=0); ax.axvspan(D + pause, 1, color=ORANGE, alpha=0.06, lw=0)
ax.text(D / 2, 3.62, "划水相：匀速（速率/力矩上限）、桨面垂直运动方向、全开、最深", ha="center", fontsize=9, color=BLUE)
ax.text(D + pause / 2, 0.04, "停顿·收扇", ha="center", va="bottom", fontsize=8.5, color=RED)
ax.text(D + pause + (1 - D - pause) / 2, 3.62, "回收相：同样快、已收拢、抬高", ha="center", fontsize=9, color=ORANGE)
for tx, lab, dy in ((0.0, "开扇：划水起步时\n（v 刚超过 U）", 0.06), (D, "收扇：v 降到 U 时\n（减速段内）", 0.06)):
    ax.annotate(lab, (tx, base["A"] + H * 1.0 + dy), xytext=(0, 14), textcoords="offset points", fontsize=8.3, color=INK2, ha="left" if tx == 0 else "right",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xlim(-0.02, 1.0); ax.set_ylim(0.0, 3.9); ax.set_yticks([]); ax.set_xticks([0, D, D + pause, 1]); ax.set_xticklabels(["0", "DUTY ≈ 0.5", "", "1"])
ax.set_xlabel("周期相位 τ", color=INK2, fontsize=10)
title2(ax, "D  结构无关的“理论最优单腿行程”（准定常 + 速率/力矩上限）", "     再用 CFD 校正非定常效应：起动附加质量、停止涡、尾流再捕获、自由面")

fig.suptitle("阻力式划水的最优行程：准定常理论能先回答什么", fontsize=13.5, x=0.06, ha="left", y=0.972, color=INK)
fig.text(0.06, 0.93, "模型：F = ½ρC_d A (v − U)²，划水相 A = A_max、回收相 A = A_min，同弧长 L；A：给定 v_p 求使 T̄ = (I_p − I_r)/T 最大的 v_r；C/D 的“最优”指速率或力矩受限下的最大冲量，效率最优（巡航）要求更大面积、更慢划水（η = U/v_p）",
         fontsize=9, color=MUTED, ha="left")
fn = os.path.join(a.out, "fig_optimal_stroke_theory.png"); plt.savefig(fn, dpi=120, facecolor=SURF); print("→", fn)
print(f"C: I_trap/I_cos = {It / Ic:.3f}, E_trap/E_cos = {Et / Ec:.3f}, net(U=0.3) = {ut / uc:.3f}")
for av in (0.67, 0.445, 0.27, 0.22):
    dv, xs = duty_opt(av); print(f"a = {av:.3f}: x* = {xs:.3f}, DUTY* = {dv:.3f}")
