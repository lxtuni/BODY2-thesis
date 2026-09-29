#!/usr/bin/env python3
"""三视图线图（全部光滑）：水下 = hull_design（解析），水线以上 = 0826 外表面（B 样条截面）。输出 cad_M3_full/fig_M3_full_lines.png"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hull_design import Design, WL, DECK
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.interpolate import UnivariateSpline
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "cad_M3_full") if os.path.isdir(os.path.join(HERE, "cad_M3_full")) else (os.path.dirname(HERE) if os.path.basename(HERE) == "scripts" else HERE)
D = Design(); S = D.summary(); O = D.O
def outline(y):
    P = np.array(O[str(float(y))]); c0 = P.mean(0); ang = np.arctan2(P[:, 1] - c0[1], P[:, 0] - c0[0]); P = P[np.argsort(ang)]; return np.vstack([P, P[:1]])
def half_at(y):
    P = np.array(O[str(float(y))]); P = P[P[:, 1] > 1e-6]; P = P[np.argsort(P[:, 0])]; xs, idx = np.unique(np.round(P[:, 0], 4), return_index=True)
    return xs, P[idx, 1]
BLUE, RED, GREY, GREEN = "#1f77b4", "#c0392b", "#7f8c8d", "#2e8b57"
fig, axs = plt.subplots(3, 1, figsize=(16, 17), gridspec_kw=dict(height_ratios=[1.25, 0.65, 1.2]))
# ---- 俯视 ----
ax = axs[0]
for y in (15.0, 10.0, 5.0, 0.0, -5.0):
    P = outline(y); ax.plot(P[:, 0], P[:, 1], "-", color=GREY, lw=0.9, alpha=0.8)
P = outline(-10.0); ax.plot(P[:, 0], P[:, 1], "-", color=RED, lw=2.4, label="对接面 / 水线 y = −10")
for y in np.arange(-15, -46, -5):
    x, b = D.waterline(float(y), 600); m = b > 0
    if m.sum() > 5: ax.plot(np.concatenate([x[m], x[m][::-1]]), np.concatenate([b[m], -b[m][::-1]]), "-", color=BLUE, lw=1.0)
xd, yd = D.deepest(); ax.axvline(xd, color=GREEN, ls=":", lw=1); ax.text(xd, 54, f"最深站位 x = {xd:.1f}\n（浮心 LCB = {S['LCB']:.1f} 与原船相同）", ha="center", va="bottom", color=GREEN, fontsize=9)
ax.annotate(f"尖端 spitz\nx = {D.xmin:.1f}（水线）/ −136（顶缘）", (D.xmin, 0), xytext=(-165, 0), textcoords="data", fontsize=9, ha="center", va="center", arrowprops=dict(arrowstyle="->"))
ax.annotate(f"钝端 stumpf\nx = {D.xmax:.1f}（水线）/ +126（顶缘）", (D.xmax, 0), xytext=(165, 0), textcoords="data", fontsize=9, ha="center", va="center", arrowprops=dict(arrowstyle="->"))
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_xlim(-190, 190); ax.set_ylim(-65, 68); ax.set_ylabel("z [mm]（横向）")
ax.set_title(f"俯视 Draufsicht — 红 = 水线 y = −10（LWL {S['LWL']:.1f}，BWL {S['BWL']:.1f}）；蓝 = 水线以下每 5 mm 的水线；灰 = 水线以上原船体轮廓 y = −5…+15（顶缘 262 × 116）", fontsize=11)
ax.legend(loc="lower right", fontsize=9)
# ---- 侧视 ----
ax = axs[1]
x = np.linspace(D.xmin, D.xmax, 800); ax.plot(x, D.keel(x), "-", color=BLUE, lw=2.2, label="龙骨线（新）")
ys = sorted([float(k) for k in O.keys()]); xl = [np.array(O[str(y)])[:, 0].min() for y in ys]; xr = [np.array(O[str(y)])[:, 0].max() for y in ys]
ax.plot(xl, ys, "-", color=GREY, lw=1.6); ax.plot(xr, ys, "-", color=GREY, lw=1.6, label="首尾轮廓（原船体上部）")
ax.plot([xl[-1], xr[-1]], [DECK, DECK], "-", color=GREY, lw=1.6); ax.axhline(WL, color=RED, ls="--", lw=1.2, label="水线 y = −10")
ax.plot([xd], [yd], "o", color=GREEN); ax.text(xd, yd - 2.5, f"最深 y = {yd:.2f}（吃水 {WL - yd:.1f} mm）", ha="center", va="top", color=GREEN, fontsize=9)
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_xlim(-150, 140); ax.set_ylim(-56, 22); ax.set_ylabel("y [mm]（向上）")
ax.set_title(f"侧视 Seitenansicht — 龙骨 D(x) = T3(1−ξ)^0.8(1+ξ)^0.85(1+aξ)，T3 = {S['T3']:.2f} mm，a = {S['a']:.3f}（首尾龙骨倾角与上部轮廓相接）；排水量 {S['V_cm3']:.1f} cm³（= 原船 385.1）；Cb = {S['Cb']:.3f}", fontsize=11)
ax.legend(loc="lower right", fontsize=9)
# ---- 横剖线图 ----
ax = axs[2]
stations = np.linspace(D.xmin + 4, D.xmax - 4, 13)
hs = {y: half_at(y) for y in ys}
for i, xs_ in enumerate(stations):
    z, y = D.section(float(xs_), 60); ax.plot(z, y, "-", color=BLUE, lw=1.0)
    # 水线以上：原船体半宽随高度
    bb = []; yy = []
    for yv in ys:
        xx, hh = hs[yv]
        if xx.min() <= xs_ <= xx.max(): bb.append(np.interp(xs_, xx, hh)); yy.append(yv)
    if len(bb) > 2:
        bb, yy = np.array(bb), np.array(yy); ax.plot(bb, yy, "-", color=GREY, lw=0.9); ax.plot(-bb, yy, "-", color=GREY, lw=0.9)
ax.axhline(WL, color=RED, ls="--", lw=1.2); ax.axhline(DECK, color=GREY, ls=":", lw=1)
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_xlim(-64, 64); ax.set_ylim(-52, 20); ax.set_xlabel("z [mm]"); ax.set_ylabel("y [mm]")
ax.set_title("横剖线图 Spantenriss — 13 站位：蓝 = 新水下剖面 b(t) = b_w(1−kt)√(1−t²)（水线处斜率与上部外飘连续，无折角）；灰 = 原船体上部（含外飘，腿槽处以外表面计）", fontsize=11)
plt.tight_layout(); fn = os.path.join(OUT, "fig_M3_full_lines.png"); plt.savefig(fn, dpi=110); print("→", fn)
