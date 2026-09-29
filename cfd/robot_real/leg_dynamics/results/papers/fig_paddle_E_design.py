#!/usr/bin/env python3
"""14 号：现有铲形桨 + E 点铰的三种做法对比与颈部铰设计图 → fig_paddle_E_design.png/.pdf, paddle_E_design.json"""
import os, json, math, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon, Circle, Wedge, FancyBboxPatch, Rectangle, FancyArrowPatch
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, SURF, GRID, WATER = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1", "#eaf2fa"
BLUE, RED, GOLD, GREEN, PUR, HULL, STEEL = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#c9c4b8", "#6b7280"

# ---------------- 几何（照片估计值，mm；机体坐标 x 向前 z 向上）
O2 = np.array([63.4, 4.0]); WL = -10.0
L_RC, L_PAR = 40.0, 60.0          # 摇臂短端、平行四边形长边（照片估计，请实测）
PSI_FIX = -90.0                    # 摇臂竖直向下
H0 = 30.0; F = 1.5                 # 设计沉浮：φ = 90° ± 30° → h0 = 60·sin30° = 30 mm
DROP, L_S = 45.0, 40.0             # 颈件：E 下折 45、向后 40 到铰
X_P = 5.0                          # 铰在铲头尖端之后 5 mm（尖端太窄，销孔放 5 mm 处）
HEAD_L, HEAD_W, THK = 42.0, 44.0, 5.6
STOP = 40.0
K_PICK, C_PICK = 3e-3, 0.1e-3

def rot(a): a = math.radians(a); return np.array([math.cos(a), math.sin(a)])
def fk(phi): C = O2 + L_RC * rot(PSI_FIX); return C, C + L_PAR * rot(PSI_FIX + phi)
DPHI = math.degrees(math.asin(H0 / L_PAR))

# ---------------- 数据
pend = pd.DataFrame([r for r in json.load(open(os.path.join(HERE, "pendulum_paddle_sweep.json"))) if "error" not in r])
act = pd.DataFrame([r for r in json.load(open(os.path.join(HERE, "paddle_active_E_A.json"))) if "error" not in r])
actB = pd.DataFrame(json.load(open(os.path.join(HERE, "paddle_active_E_bigh0.json"))))
neck = pd.DataFrame([r for r in json.load(open(os.path.join(HERE, "passive_spade_sweep.json"))) if "error" not in r])
sel = neck[(neck.f == F) & (neck.x_p == X_P) & np.isclose(neck.h0, H0 * 1e-3)].sort_values("k")
best = sel.iloc[(sel.k - K_PICK).abs().argsort().iloc[0]]
sel2 = neck[(neck.f == 2.0) & (neck.x_p == X_P) & np.isclose(neck.h0, H0 * 1e-3)]
best2 = sel2.iloc[(sel2.k - K_PICK).abs().argsort().iloc[0]]
PURE = {1.5: 11.6, 2.0: 21.8}   # 刚性纯沉浮 h0 20（UVLM）

fig = plt.figure(figsize=(18, 12.5), facecolor=SURF)
gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1], height_ratios=[1.2, 1], left=0.035, right=0.985, top=0.86, bottom=0.085, hspace=0.32, wspace=0.17)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, 0]); axE = fig.add_subplot(gs[1, 1]); axF = fig.add_subplot(gs[1, 2])
for ax in (axA, axB, axC, axD, axE, axF):
    ax.set_facecolor(SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8.5); ax.grid(color=GRID, lw=0.6)

def dim(ax, p0, p1, text, off=0, fs=8, color=INK2):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float); d = p1 - p0; n = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-9)
    a, b = p0 + n * off, p1 + n * off
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="<->", color=color, lw=0.8, shrinkA=0, shrinkB=0))
    ax.text(*((a + b) / 2 + n * 3), text, fontsize=fs, color=color, ha="center", va="center")

# ============ A 侧视总图
ax = axA; ax.set_aspect("equal"); ax.grid(False)
ax.axhspan(-140, WL, color=WATER, zorder=0); ax.axhline(WL, color=BLUE, lw=1, ls="--"); ax.text(180, WL + 3, "水线", color=BLUE, fontsize=8)
# 船体（示意）
hull = Polygon([(-60, 38), (200, 38), (215, 22), (200, 0), (170, WL - 6), (-40, WL - 6), (-60, 12)], closed=True, fc=HULL, ec=STEEL, lw=0.8, alpha=0.6, zorder=1)
ax.add_patch(hull); ax.text(0, 22, "船体（示意）", fontsize=8.5, color=INK2, ha="center")
def draw_leg(ax, phi, th_deg, alpha=1.0, lw=2.2, label=False):
    C, E = fk(phi)
    D = O2 + L_PAR * rot(PSI_FIX + phi)
    for a, b in ((O2, C), (O2, D), (C, E), (D, E)):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=STEEL, lw=lw, alpha=alpha, solid_capstyle="round", zorder=3)
    # 颈件：E → 下折 DROP → 向后 L_S → 铰 P
    K = E + np.array([0, -DROP]); P = K + np.array([-L_S, 0])
    ax.plot([E[0], K[0], P[0]], [E[1], K[1], P[1]], color=GOLD, lw=lw + 1.2, alpha=alpha, solid_capstyle="round", zorder=4)
    # 铲头（侧视：厚 5.6 的板，绕 P 转 θ，头向后拖）
    u = -rot(th_deg)            # 弦向（向后），θ>0 = 前缘抬头 → 尾缘向下
    n = np.array([-u[1], u[0]])
    poly = [P + X_P * u - n * THK / 2, P + (X_P + HEAD_L) * u - n * THK / 2, P + (X_P + HEAD_L) * u + n * THK / 2, P + X_P * u + n * THK / 2]
    ax.add_patch(Polygon(poly, closed=True, fc=RED, ec=RED, alpha=alpha * 0.9, zorder=5))
    ax.add_patch(Circle(P, 2.2, fc="white", ec=INK, lw=1, alpha=alpha, zorder=6))
    if label:
        ax.add_patch(Circle(E, 3.2, fc=BLUE, ec="white", lw=1.2, zorder=7)); ax.text(E[0] - 2, E[1] + 7, "E（蓝点）", color=BLUE, fontsize=9, weight="bold", ha="right")
        ax.plot(*O2, "o", color=INK, ms=5, zorder=7); ax.text(O2[0] - 6, O2[1] + 6, "髋 O2", fontsize=8.5, color=INK, ha="right")
        ax.text(*(C + np.array([-6, -4])), "C", fontsize=8, color=INK2, ha="right"); ax.text(*(D + np.array([5, 4])), "D", fontsize=8, color=INK2)
        ax.text(-62, -2, "ψ 舵机（髋）固定在 −90°：摇臂 O2→C 竖直向下\nφ 舵机经平行四边形 O2-C-E-D 上下抬 E\n颈件（金色）替换原来的杆，铲头（红）不变", fontsize=8, color=INK2, va="top")
        ax.text(P[0] - 2, P[1] - 10, "颈部铰 P\n（扭簧 + 限位）", fontsize=8.5, color=INK, ha="center", va="top")
        ax.add_patch(Wedge(P, HEAD_L + 6, 180 - STOP, 180 + STOP, fc=GREEN, alpha=0.13, ec=GREEN, lw=0.8, zorder=2))
        ax.text(P[0] - HEAD_L - 12, P[1] + 20, f"被动摆角\n±{STOP:.0f}° 限位\n（工作 ±{best.th0:.0f}°）", fontsize=8, color=GREEN, ha="right")
    return C, D, E, K, P
# 相位残影
for ph, a in ((0.25, 0.25), (0.75, 0.25)):
    phi = 90 + DPHI * math.sin(2 * math.pi * ph); th = best.th0 * math.sin(2 * math.pi * ph + math.radians(best.psi))
    draw_leg(ax, phi, th, alpha=a, lw=1.6)
C, D, E, K, P = draw_leg(ax, 90.0, 0.0, label=True); P_HINGE = P.copy()
# 用户的红色扇区（原概念）
ax.add_patch(Wedge(E, 70, 180, 275, fc=RED, alpha=0.06, ec=RED, lw=1.0, ls=":", zorder=2))
ax.text(E[0] + 45, E[1] - 78, "你画的红扇区：\n整根桨绕 E 被动摆\n→ 推力 ≤ 纯沉浮（图 B）", fontsize=8, color=RED, ha="center", va="top")
# 尺寸
dim(ax, E, K, f"下折 {DROP:.0f}", off=8); dim(ax, K, P, f"{L_S:.0f}", off=-8)
dim(ax, (P[0], P[1] - 36), (P[0] - HEAD_L, P[1] - 36), f"铲头 {HEAD_L:.0f}×{HEAD_W:.0f}（不变）", off=0)
Et, Eb = fk(90 + DPHI)[1], fk(90 - DPHI)[1]
ax.annotate("", xy=(Et[0] + 40, Et[1]), xytext=(Eb[0] + 40, Eb[1]), arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.2))
ax.text(Et[0] + 45, (Et[1] + Eb[1]) / 2, f"E 沉浮 ±{H0:.0f}\nφ = 90°±{DPHI:.0f}°\n{F} Hz", fontsize=8.5, color=BLUE, va="center")
ax.text(P[0] - HEAD_L / 2, P[1] - 50, f"铰轴心水线下 {WL - P[1]:.0f} mm；冲程顶端头离水面 ≈ {WL - (P[1] + H0 + HEAD_L * math.sin(math.radians(best.th0))):.0f} mm", fontsize=8, color=INK2, ha="center")
ax.annotate("", xy=(-30, -125), xytext=(20, -125), arrowprops=dict(arrowstyle="->", color=INK2)); ax.text(-5, -121, "来流 U = 0.25 m/s（机体向前 →）", fontsize=8, color=INK2, ha="center")
ax.set_xlim(-70, 235); ax.set_ylim(-140, 50); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("A  设计侧视：ψ 固定向下，φ 单舵机沉浮；颈件下折 + 铲头颈部扭簧铰（头朝后拖）", fontsize=10.5, loc="left", color=INK)

# ============ B 被动铰在 E
ax = axB
for L, col in ((30, GOLD), (50, PUR), (82, RED)):
    for f, ls in ((1.5, "-"), (2.0, "--")):
        s = pend[(pend.L_stem * 1e3).round() .eq(L) & (pend.f == f) & np.isclose(pend.c, 0.3e-3)].sort_values("k")
        ax.plot(s.k * 1e3, s["T"] * 1e3, ls, color=col, lw=1.8, marker="o", ms=3.5, label=f"杆 {L} mm, {f} Hz")
for f, ls in ((1.5, "-"), (2.0, "--")):
    ax.axhline(PURE[f], color=INK2, lw=1, ls=ls); ax.text(51, PURE[f] + 0.4, f"刚性纯沉浮 {f} Hz", fontsize=7.5, color=INK2, ha="right")
ax.axhline(0, color=INK, lw=0.6)
ax.set_xlabel("扭簧刚度 k [mN·m/rad]"); ax.set_ylabel("推力 / 桨 [mN]")
ax.set_title("B  你的概念：整桨绕 E 被动摆（h0 20 mm，c 0.3）\n俯仰滞后 70–90°，任何 k 都不超过刚性纯沉浮", fontsize=10, loc="left", color=INK)
ax.legend(fontsize=7.5, ncol=2, frameon=False, loc="upper left")
ax.text(0.98, 0.04, "软簧 → 负推力（桨随流摆，成阻力）", transform=ax.transAxes, fontsize=8, color=RED, ha="right")

# ============ C 主动（ψ 驱动）铰在 E
ax = axC
for L, col in ((50, PUR), (82, RED)):
    s = act[np.isclose(act.L, L * 1e-3)]
    ax.scatter(s.alpha_max, s["T"] * 1e3, s=14, color=col, alpha=0.55, label=f"杆 {L} mm（h0 15/20，1.5–2 Hz）")
ax.scatter(actB.alpha_max, actB["T"] * 1e3, s=22, color=GREEN, marker="^", label="杆 50，h0 30–40（现连杆达不到）")
ax.axvspan(0, 30, color=GREEN, alpha=0.07); ax.axvline(30, color=GREEN, lw=1); ax.text(29, ax.get_ylim()[1] * 0.02 + 150, "α ≤ 30°\n不失速区", fontsize=8, color=GREEN, ha="right", va="bottom")
ok = act[act.alpha_max <= 36]; ax.annotate(f"α ≤ 36° 内最大 {ok['T'].max()*1e3:.0f} mN", xy=(36, ok['T'].max() * 1e3), xytext=(45, 60), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="->", color=INK2))
ax.text(0.98, 0.96, "α > 40° 的点是 UVLM 无失速假象：\n真实是拍水阻力型，效率 ~0.25", transform=ax.transAxes, fontsize=8, color=INK2, ha="right", va="top")
ax.set_xlabel("头部最大运动学攻角 α_max [°]"); ax.set_ylabel("推力 / 桨 [mN]"); ax.set_ylim(-5, 215)
ax.set_title("C  改进版：铰在 E，ψ 舵机经硬簧主动俯仰\n远铰 → 抬头即抬头部，相位锁死，不失速只有 10–17 mN", fontsize=10, loc="left", color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="center right")

# ============ D 颈部铰被动：k 扫描
ax = axD
for h0, col, ls in ((0.02, MUTED, "--"), (0.03, RED, "-")):
    for xp, mk in ((-5.0, "s"), (0.0, "o"), (5.0, "^")):
        s = neck[(neck.f == F) & (neck.x_p == xp) & np.isclose(neck.h0, h0)].sort_values("k")
        if len(s): ax.plot(s.k * 1e3, s["T"] * 1e3, ls, color=col, marker=mk, ms=4, lw=1.4, label=f"h0 {h0*1e3:.0f}, 铰 {xp:+.0f} mm")
ax.plot(best.k * 1e3, best["T"] * 1e3, "o", ms=13, mfc="none", mec=GREEN, mew=2)
ax.set_xlabel("扭簧刚度 k [mN·m/rad]"); ax.set_ylabel("推力 / 桨 [mN]", color=INK)
ax2 = ax.twinx(); ax2.set_facecolor("none")
for xp, mk in ((0.0, "o"),):
    s = neck[(neck.f == F) & (neck.x_p == xp) & np.isclose(neck.h0, 0.03)].sort_values("k")
    ax2.plot(s.k * 1e3, s.alpha_max, ":", color=BLUE, marker=mk, ms=3, lw=1.2, label="α_max（h0 30, 铰 0）")
    ax2.plot(s.k * 1e3, s.eta * 100, ":", color=GOLD, marker=mk, ms=3, lw=1.2, label="η ×100")
ax2.set_ylim(0, 100); ax2.set_ylabel("α_max [°]  /  η [%]", color=INK2); ax2.tick_params(colors=INK2, labelsize=8.5)
for s_ in ("top",): ax2.spines[s_].set_visible(False)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=7.2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4)
ax.set_title(f"D  推荐：颈部扭簧铰（头朝后拖）UVLM 流固耦合，{F} Hz，c 0.1\n设计点 k {best.k*1e3:.0f}：θ0 ±{best.th0:.0f}° 相位 {best.psi:.0f}° α {best.alpha_max:.0f}° T {best['T']*1e3:.0f} mN η {best.eta:.2f}", fontsize=10, loc="left", color=INK)

# ============ E 时程
ax = axE
t = np.linspace(0, 2 / F, 400); w = 2 * math.pi * F
phi = 90 + DPHI * np.sin(w * t); z = H0 * np.sin(w * t)
th = best.th0 * np.sin(w * t + math.radians(best.psi))
gam = np.degrees(np.arctan2(H0 * 1e-3 * w * np.cos(w * t), 0.25)); al = gam - th
ax.plot(t * F, phi - 90, color=STEEL, lw=2, label=f"φ − 90° [°]（φ 舵机，峰值 {DPHI*w:.0f} °/s）")
ax.plot(t * F, z, color=BLUE, lw=2, label="E 沉浮 z [mm]")
ax.plot(t * F, th, color=RED, lw=2, label=f"被动俯仰 θ [°]（领先 {best.psi:.0f}°）")
ax.plot(t * F, al, color=GREEN, lw=1.5, ls="--", label=f"头部攻角 α [°]（峰 {best.alpha_max:.0f}°）")
ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel("t / T"); ax.set_ylim(-50, 50)
ax.legend(fontsize=7.5, frameon=False, loc="upper right", ncol=2)
ax.set_title(f"E  一个设计周期（{F} Hz）：ψ 不动，只有 φ 在动，俯仰全靠弹簧", fontsize=10, loc="left", color=INK)

# ============ F 对比
ax = axF
rows = [("现在：铲桨阻力划水\n(V3trap 最优, 4 号)", 87, 21, MUTED),
        ("你的概念：整桨绕 E\n被动扭簧（最好 k）", float(pend[(pend.f == 1.5)]["T"].max() * 1e3), float(pend[(pend.f == 1.5)].loc[pend[(pend.f == 1.5)]["T"].idxmax(), "P"] * 1e3), RED),
        ("改进：铰在 E，ψ 主动\n经硬簧（α ≤ 36°）", float(ok["T"].max() * 1e3), float(ok.loc[ok["T"].idxmax(), "P"] * 1e3), PUR),
        (f"推荐：颈部扭簧铰\n头朝后 h0 30 {F} Hz", float(best["T"] * 1e3), float(best["P"] * 1e3), GREEN),
        (f"同上 2 Hz", float(best2["T"] * 1e3), float(best2["P"] * 1e3), GREEN),
        ("参考：13 号月牙脚\n(18.5 cm², h0 30)", 69, 21, GOLD)]
y = np.arange(len(rows))[::-1]
for yi, (lab, T, P, col) in zip(y, rows):
    ax.barh(yi, T, color=col, alpha=0.85, height=0.62); ax.text(T + 1.5, yi, f"{T:.0f} mN  ·  {P:.0f} mW", va="center", fontsize=8.5, color=INK)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8); ax.set_xlabel("单桨周期平均推力 [mN]（1.5 Hz，铲头 12.9 cm²，UVLM）"); ax.set_xlim(0, 115)
ax.axvspan(21 / 4, 36 / 4, color=BLUE, alpha=0.12); ax.text(36 / 4 + 1, y[0] + 0.45, "整机巡航需求 / 4 桨（腿整流后）", fontsize=7.5, color=BLUE, va="center")
ax.set_title("F  三种做法对比（同一铲头、同一舵机）", fontsize=10, loc="left", color=INK)
ax.grid(axis="y", visible=False)

fig.text(0.035, 0.955, "14 · 现有铲形桨 + 平行四边形末端 E 铰：能不能做成升力型？", fontsize=16, weight="bold", color=INK)
fig.text(0.035, 0.925, f"结论：整根桨绕 E 被动摆 → 俯仰滞后、推力 ≤ 刚性纯沉浮；把弹簧铰移到铲头颈部、头朝后拖、E 用颈件下折 {DROP:.0f} mm → 被动俯仰成立，"
         f"{F} Hz 单桨 {best['T']*1e3:.0f} mN（η {best.eta:.2f}），2 Hz {best2['T']*1e3:.0f} mN。铲头小（12.9 cm²）+ 沉浮 ±30 是上限；"
         "照片桨头更大（≈64×50）→ 推力约 ×2。连杆尺寸为照片估计，请实测后重代。", fontsize=9.5, color=INK2, wrap=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"fig_paddle_E_design.{ext}"), dpi=160 if ext == "png" else None, facecolor=SURF)
json.dump(dict(L_RC=L_RC, L_PAR=L_PAR, PSI_FIX=PSI_FIX, H0=H0, F=F, DPHI=DPHI, DROP=DROP, L_S=L_S, X_P=X_P, STOP=STOP, K=K_PICK, C=C_PICK,
               E_center=list(map(float, fk(90)[1])), P_center=list(map(float, P_HINGE)), best=best.to_dict(), best2=best2.to_dict(),
               pend_best_T=rows[1][1], active_best_T=rows[2][1]), open(os.path.join(HERE, "paddle_E_design.json"), "w"), indent=1, default=float)
print("saved; best:", best.to_dict()); print("best2:", best2.to_dict())
