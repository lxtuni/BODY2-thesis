#!/usr/bin/env python3
"""文献对照图：周期相位占比（BODY2 vs 狗 vs 方案 A）+ CFD 单腿推力随髋角。  python3 make_lit_figure.py [results_dir]"""
import sys, os, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_leg")
def _cjk():
    import glob as g
    c = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + g.glob("/mnt/c/Windows/Fonts/msyh.ttc")
    for f in c:
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return
        except Exception: pass
_cjk(); plt.rcParams["axes.unicode_minus"] = False
INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"
T = 1.25; DUTY = 0.5; PSI0 = 136; SWEEP = 60
d = np.loadtxt(os.path.join(R, "force_open.dat"), comments="#"); t = d[:, 0]; Fx = d[:, 1] * 1e3
dc = np.loadtxt(os.path.join(R, "force_closed.dat"), comments="#"); tc = dc[:, 0]; Fxc = dc[:, 1] * 1e3
def smooth(tt, y, w=0.02):
    dt = np.median(np.diff(tt)); n = max(1, int(round(w / dt))); return np.convolve(y, np.ones(n) / n, mode="same")
def cyc(tt, y, n=3):
    m = (tt >= (n - 1) * T) & (tt < (n - 1) * T + DUTY * T); tau = (tt[m] - (n - 1) * T) / T
    return tau, PSI0 - SWEEP * (1 - np.cos(np.pi * tau / DUTY)) / 2, smooth(tt, y)[m]
tau, psi, fo = cyc(t, Fx); tauc, psic, fc = cyc(tc, Fxc)
i0 = np.where((fo[:-1] > 0) & (fo[1:] <= 0) & (tau[:-1] > 0.3))[0]; psi_zero = psi[i0[0]]; ip = np.argmax(fo)
psi_peak = PSI0 - SWEEP * (1 - np.cos(np.pi * 0.254 / DUTY)) / 2          # 与 analyze_stroke.py 的 peak_tau 一致
fig, axs = plt.subplots(1, 2, figsize=(13, 4.3), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.05, 1]))
ax = axs[0]; ax.set_facecolor(SURF)
rows = [("BODY2 现步态 (CFD 基准)", [(0.50, INK, "划水 50 %"), (0.15, ORANGE, "收拢切换"), (0.20, GRAY, "蜷缩前扫"), (0.15, ORANGE, "伸展")]),
        ("狗 前肢 (Fish 2020)", [(0.37, INK, "划水 34–39 %"), (0.63, GRAY, "回收（屈曲、内收、慢）")]),
        ("狗 后肢 (Fish 2020)", [(0.30, INK, "划水 26–34 %"), (0.70, GRAY, "回收")]),
        ("方案 A：狗式占空比 0.35，T 不变", [(0.35, GREEN, "划水 35 %（速率 ×1.43）"), (0.65, GRAY, "慢蜷缩回收 65 %")])]
for i, (lab, segs) in enumerate(rows):
    x = 0; y = len(rows) - 1 - i
    for w_, col, txt in segs:
        ax.barh(y, w_, left=x, height=0.58, color=col, alpha=0.9 if col != GRAY else 0.55, edgecolor=SURF, linewidth=1.5)
        if w_ > 0.12: ax.text(x + w_ / 2, y, txt, ha="center", va="center", fontsize=8.6, color="white" if col in (INK, GREEN, BLUE) else INK)
        x += w_
    ax.text(-0.02, y, lab, ha="right", va="center", fontsize=9.5, color=INK)
ax.set_xlim(0, 1); ax.set_ylim(-0.6, len(rows) - 0.4); ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(["0", "25 %", "50 %", "75 %", "100 %"])
ax.set_xlabel("周期相位 τ"); ax.set_title("划水/回收占比：狗用 ~1/3 划水、~2/3 回收；BODY2 现为 1/2", loc="left", fontsize=10.5)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.grid(axis="x", color=GRID, lw=0.8)
ax = axs[1]; ax.set_facecolor(SURF)
ax.plot(psi, fo, color=INK, lw=2.0, label="open 桨（划水态）CFD 第 3 周期"); ax.plot(psic, fc, color=BLUE, lw=1.4, ls="--", label="closed 桨（收拢态）")
ax.axhline(0, color=GRAY, lw=0.8); ax.axvline(90, color=RED, lw=1.2, ls=":")
ax.text(91, 38, "狗：肢体与躯干垂直处\n结束划水 (ψ = 90°)", fontsize=8.6, color=RED, ha="right", va="top")
ax.axvspan(76, psi_zero, color=RED, alpha=0.08, lw=0); ax.text(84.5, -30, f"F_x < 0：{psi_zero:.0f}° → 76°", ha="center", fontsize=8.6, color=RED)
ax.annotate(f"峰值 {fo[ip]:.0f} mN @ ψ ≈ {psi_peak:.0f}°", (psi[ip], fo[ip]), xytext=(8, 4), textcoords="offset points", fontsize=8.8, color=INK)
ax.annotate(f"过零 ψ ≈ {psi_zero:.0f}°", (psi_zero, 0), xytext=(-70, 26), textcoords="offset points", fontsize=8.8, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.set_xlim(140, 72); ax.set_ylim(-40, 130); ax.set_xlabel("髋角 ψ [°]（136° 前伸 → 76° 后扫；90° = 桨竖直）"); ax.set_ylabel("单腿推力 F_x [mN]")
ax.set_title(f"推力随髋角：峰值在 {psi_peak:.0f}°，过零在 {psi_zero:.0f}°，与狗停在 90° 一致", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.8, loc="upper left")
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.grid(color=GRID, lw=0.8)
plt.tight_layout(); plt.savefig(os.path.join(R, "fig_lit_compare.png"), dpi=120, facecolor=SURF); print("→", os.path.join(R, "fig_lit_compare.png"))
