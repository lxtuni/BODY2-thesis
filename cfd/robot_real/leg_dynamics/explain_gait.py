#!/usr/bin/env python3
"""
步态参数图解：把 leg_dynamics.py 的步态参数（T、DUTY、SWEEP、PSI_START、PHI_EXT、PHI_RET、SWITCH_FRAC、FAN_DEG、FAN_CLOSE、CLOSED_W）
画在一张图上，并打印派生量（划水/回收时长、蜷缩/伸展时长、髋角速度峰值、桨尖速度与深度、桨面积、舵机指令范围）。
    T=0.81 DUTY=0.45 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5 python3 explain_gait.py [--name V3] [--out results_verify]
不带环境变量 = 现步态默认值。输出 <out>/fig_gait_params_<name>.png 与 gait_params_<name>.json
"""
import os, sys, json, argparse, glob, importlib.util
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--name", default="V3"); ap.add_argument("--out", default=None); a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = a.out or os.path.join(HERE, "results_verify"); os.makedirs(OUT, exist_ok=True)
os.environ.setdefault("OUT", OUT)
sp = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(sp); sp.loader.exec_module(ld)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc, FancyArrowPatch
INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"
PW, RC = "#2a78d6", "#eb6834"                      # 划水相 / 回收相

T, DUTY, SWEEP, PSI0 = ld.T_PER, ld.DUTY, ld.SWEEP, np.degrees(ld.PSI0)
PE, PR, SF = ld.PHI_EXT, ld.PHI_RET, min(max(ld.SWITCH_FRAC, 0.05), 0.5)
FAN, FC, CW = ld.FAN_DEG, ld.FAN_CLOSE, ld.CLOSED_W
N = 2000; tau = np.linspace(0, 1, N, endpoint=False); t = tau * T
psi, phi, power = ld.gait(t); ff = ld.fan_fraction(t, power); psi_d, phi_d = np.degrees(psi), np.degrees(phi)
om = np.gradient(np.unwrap(psi), t); PP = [ld.pose(p, q) for p, q in zip(psi, phi)]
E = np.array([p["E"] for p in PP]); U_ = np.array([p["u"] for p in PP]); tip = E + (ld.STEM + ld.HEAD) * U_
vtip = np.linalg.norm(np.gradient(tip, t, axis=0), axis=1); WL = -0.01
tp, tr = DUTY * T, (1 - DUTY) * T
A_c = ld.STEM * ld.STEM_W + ld.HEAD_W * ld.HEAD - (1 - np.pi / 4) * ld.HEAD_W * 0.025 / 2; A_o = ld.STEM * ld.STEM_W + 0.5 * ld.FAN_R ** 2 * np.deg2rad(FAN)
th1 = np.degrees([ld.crank1_from_rocker(p + np.pi) for p in psi]); th2 = np.degrees(psi + phi)
D = dict(name=a.name, T=T, f=1 / T, DUTY=DUTY, SWEEP=SWEEP, PSI_START=PSI0, PSI_END=PSI0 - SWEEP, PSI_CENTER=PSI0 - SWEEP / 2, PHI_EXT=PE, PHI_RET=PR, SWITCH_FRAC=SF,
         FAN_DEG=FAN, FAN_CLOSE=FC, CLOSED_W_mm=CW * 1e3, t_power=tp, t_recovery=tr, t_tuck=SF * tr, t_extend=SF * tr, t_hold_tucked=(1 - 2 * SF) * tr,
         w_power_max_deg_s=float(np.degrees(np.abs(om[power])).max()), w_rec_max_deg_s=float(np.degrees(np.abs(om[~power])).max()),
         w_power_mean_deg_s=SWEEP / tp, vtip_max=float(vtip.max()), vtip_power_mean=float(vtip[power].mean()),
         tip_depth_power_mm=[float((WL - tip[power, 1]).min() * 1e3), float((WL - tip[power, 1]).max() * 1e3)],
         tip_z_range_mm=[float(tip[:, 1].min() * 1e3), float(tip[:, 1].max() * 1e3)], E_lift_mm=float((E[:, 1].max() - E[:, 1].min()) * 1e3),
         tip_x_range_mm=[float(tip[:, 0].min() * 1e3), float(tip[:, 0].max() * 1e3)],
         A_closed_cm2=A_c * 1e4, A_open_cm2=A_o * 1e4, fan_tip_width_mm=2 * ld.FAN_R * np.sin(np.deg2rad(FAN) / 2) * 1e3,
         servo2_crank1_deg=[float(np.nanmin(th1)), float(np.nanmax(th1))], servo1_crank2_deg=[float(th2.min()), float(th2.max())])
json.dump(D, open(os.path.join(OUT, f"gait_params_{a.name}.json"), "w"), indent=1, ensure_ascii=False)
print(f"[{a.name}] T {T} s（{1/T:.2f} Hz）  划水 {tp:.3f} s / 回收 {tr:.3f} s（DUTY {DUTY}）  ψ {PSI0:.0f}→{PSI0-SWEEP:.0f}°（中心 {PSI0-SWEEP/2:.0f}°，SWEEP {SWEEP:.0f}°）")
print(f"     φ 伸展 {PE:.0f}° / 蜷缩 {PR:.0f}°  蜷缩用 {SF*tr:.3f} s，蜷着 {(1-2*SF)*tr:.3f} s，伸展用 {SF*tr:.3f} s（SWITCH_FRAC {SF}）")
print(f"     髋角速度峰值：划水 {D['w_power_max_deg_s']:.0f}°/s（平均 {SWEEP/tp:.0f}°/s）  回收 {D['w_rec_max_deg_s']:.0f}°/s   桨尖速度峰 {vtip.max():.2f} m/s，划水相均 {vtip[power].mean():.2f} m/s")
print(f"     桨尖水下深度（划水相）{D['tip_depth_power_mm'][0]:.0f}–{D['tip_depth_power_mm'][1]:.0f} mm  回收抬起 E 上升 {D['E_lift_mm']:.0f} mm  桨尖 z {D['tip_z_range_mm'][0]:.0f}…{D['tip_z_range_mm'][1]:.0f} mm（水线 −10）")
print(f"     桨面积 收拢 {A_c*1e4:.1f} cm² / 展开 {A_o*1e4:.1f} cm²（扇 {FAN:.0f}°，尖端宽 {D['fan_tip_width_mm']:.0f} mm）  舵机：servo_2(曲柄1) {D['servo2_crank1_deg'][0]:.0f}…{D['servo2_crank1_deg'][1]:.0f}°  servo_1(曲柄2) {D['servo1_crank2_deg'][0]:.0f}…{D['servo1_crank2_deg'][1]:.0f}°")

# ------------------------------------------------------------------ 图
fig = plt.figure(figsize=(15.5, 9.4), facecolor=SURF)
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.3], height_ratios=[1, 1], hspace=0.45, wspace=0.42)
ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[0, 1])
ax4 = fig.add_subplot(gs[1, 1]); ax5 = fig.add_subplot(gs[:, 2])
def style(ax, xlab="周期相位 τ = t / T"):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8); ax.set_xlabel(xlab)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
def shade(ax):
    ax.axvspan(0, DUTY, color=PW, alpha=0.07, lw=0); ax.axvspan(DUTY, 1, color=RC, alpha=0.07, lw=0)
    ax.axvline(DUTY, color=INK2, lw=0.9, ls="--")
# ① ψ(τ)
style(ax1); shade(ax1)
ax1.plot(tau, psi_d, color=INK, lw=2)
ax1.axhline(PSI0, color=GRAY, lw=0.8, ls=":"); ax1.axhline(PSI0 - SWEEP, color=GRAY, lw=0.8, ls=":")
ax1.annotate("", xy=(0.02, PSI0 - SWEEP), xytext=(0.02, PSI0), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2))
ax1.text(0.035, PSI0 - SWEEP / 2, f"SWEEP = {SWEEP:.0f}°", color=RED, fontsize=9.5, va="center")
ax1.text(0.0, PSI0 + 1.5, f"PSI_START = {PSI0:.0f}°（划水起点，桨在前）", fontsize=8.5, color=INK2, va="bottom", zorder=6)
ax1.text(DUTY + 0.02, PSI0 - SWEEP - 2, f"{PSI0-SWEEP:.0f}° = PSI_START − SWEEP（桨在后）", fontsize=8.5, color=INK2, va="top")
ax1.text(0.14, PSI0 - SWEEP + 1, f"划水相\nDUTY·T = {tp:.2f} s", ha="center", va="bottom", color=PW, fontsize=9.5)
ax1.text(0.82, PSI0 - SWEEP + 6, f"回收相\n(1−DUTY)·T = {tr:.2f} s", ha="center", va="bottom", color=RC, fontsize=9.5)
ax1.text(DUTY - 0.01, PSI0 + 9.5, f"DUTY = {DUTY}", color=INK2, fontsize=9, ha="right", va="top")
ax1.set_ylabel("摇臂角 ψ [°]（O2→C，从 +x 逆时针；90° = 摇臂竖直）"); ax1.set_title(f"① 髋摆动 ψ(τ)：T = {T} s（{1/T:.2f} Hz），余弦律，峰值 {D['w_power_max_deg_s']:.0f}°/s", loc="left", fontsize=10)
ax1.set_xlim(0, 1); ax1.set_ylim(PSI0 - SWEEP - 12, PSI0 + 14)
# ② φ(τ)
style(ax2); shade(ax2)
ax2.plot(tau, phi_d, color=INK, lw=2)
ax2.axhline(PE, color=GRAY, lw=0.8, ls=":"); ax2.axhline(PR, color=GRAY, lw=0.8, ls=":")
ax2.text(0.01, PE + 2, f"PHI_EXT = {PE:.0f}°（伸展：桨垂在最下）", fontsize=9, color=INK2, va="bottom")
ax2.text(0.99, PR - 2.5, f"PHI_RET = {PR:.0f}°（蜷缩：桨抬起 {D['E_lift_mm']:.0f} mm）", fontsize=9, color=INK2, va="top", ha="right")
t1, t2 = DUTY + SF * (1 - DUTY), 1 - SF * (1 - DUTY)
ax2.annotate("", xy=(t1, PR + 12), xytext=(DUTY, PR + 12), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2))
ax2.text((DUTY + t1) / 2, PR + 15, f"蜷缩 SWITCH_FRAC×回收 = {SF*tr:.2f} s", ha="center", color=RED, fontsize=9)
ax2.annotate("", xy=(1, PR + 26), xytext=(t2, PR + 26), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2))
ax2.text((1 + t2) / 2, PR + 29, f"伸展 {SF*tr:.2f} s", ha="center", color=RED, fontsize=9)
if (1 - 2 * SF) > 0.02: ax2.text((t1 + t2) / 2, PR + 4, f"保持蜷缩 {(1-2*SF)*tr:.2f} s", ha="center", color=INK2, fontsize=9)
ax2.set_ylabel("平行四边形张角 φ = θ2 − ψ [°]"); ax2.set_title(f"② 伸展/蜷缩 φ(τ)：SWITCH_FRAC = {SF}（{'慢蜷缩，收完立刻放' if SF >= 0.5 else '前 %.0f %% 收、后 %.0f %% 放' % (100*SF, 100*SF)}）", loc="left", fontsize=10)
ax2.set_xlim(0, 1); ax2.set_ylim(PR - 12, PE + 12)
# ③ 脚蹼展开度 + 两态合成
style(ax3); shade(ax3)
ax3.plot(tau, ff, color=INK, lw=2)
tc = DUTY if FC < 0 else FC
ax3.text(0.02, 1.05, "划水开始：线把扇拉开（过渡 8 % T）", fontsize=8.5, color=INK2)
ax3.text(tc + 0.02, 0.62, f"FAN_CLOSE = {'−1' if FC < 0 else FC}\n{'→ 划水结束（τ = DUTY）' if FC < 0 else ''}\n松线自动收拢", fontsize=8.5, color=INK2)
ax3.text(DUTY / 2, 0.22, f"CFD 取 open 算例\n扇 {FAN:.0f}°，{A_o*1e4:.1f} cm²", ha="center", color=PW, fontsize=8.5)
ax3.text(DUTY + (1 - DUTY) / 2 + 0.03, 0.22, f"CFD 取 closed 算例\n{CW*1e3:.0f} mm 宽，{A_c*1e4:.1f} cm²", ha="center", color=RC, fontsize=8.5)
ax3.set_ylabel("脚蹼展开度 0…1"); ax3.set_title("③ 脚蹼：FAN_DEG / FAN_CLOSE / CLOSED_W，两态合成规则", loc="left", fontsize=10); ax3.set_xlim(0, 1); ax3.set_ylim(-0.05, 1.2)
# ④ 桨的正视图：收拢 vs 展开
style(ax4, xlab="y [mm]（横向）"); ax4.set_aspect("equal")
s = np.linspace(0, 1, 30); sh = [(CW / 2 * np.sin(np.pi / 2 * k), ld.STEM + 0.025 * (1 - np.cos(np.pi / 2 * k))) for k in s]
outl = [(-ld.STEM_W / 2, 0), (ld.STEM_W / 2, 0), (ld.STEM_W / 2, ld.STEM)] + [(w, z) for w, z in sh if w > ld.STEM_W / 2] + [(CW / 2, ld.STEM + ld.HEAD), (-CW / 2, ld.STEM + ld.HEAD)] + [(-w, z) for w, z in sh[::-1] if w > ld.STEM_W / 2] + [(-ld.STEM_W / 2, ld.STEM)]
ang = np.deg2rad(FAN); th = np.linspace(ang / 2, -ang / 2, 40); fan = [(ld.FAN_R * np.sin(q), ld.STEM + ld.FAN_R * np.cos(q)) for q in th]
outl_o = [(-ld.STEM_W / 2, 0), (ld.STEM_W / 2, 0), (ld.STEM_W / 2, ld.STEM)] + fan + [(-ld.STEM_W / 2, ld.STEM)]
off = 0.055
ax4.add_patch(Polygon([(x * 1e3 - off * 1e3, -z * 1e3) for x, z in outl], closed=True, fc=RC, ec=INK, lw=1, alpha=0.55))
ax4.add_patch(Polygon([(x * 1e3 + off * 1e3, -z * 1e3) for x, z in outl_o], closed=True, fc=PW, ec=INK, lw=1, alpha=0.55))
ax4.text(-off * 1e3, 8, f"收拢（closed）\nCLOSED_W = {CW*1e3:.0f} mm\n{A_c*1e4:.1f} cm²", ha="center", fontsize=8.5, color=INK2, va="bottom")
ax4.text(off * 1e3, 8, f"展开（open）\nFAN_DEG = {FAN:.0f}°，R = {ld.FAN_R*1e3:.0f} mm\n尖端宽 {D['fan_tip_width_mm']:.0f} mm\n{A_o*1e4:.1f} cm²（×{A_o/A_c:.2f}）", ha="center", fontsize=8.5, color=INK2, va="bottom")
ax4.annotate("", xy=(-off * 1e3 - 22, -ld.STEM * 1e3), xytext=(-off * 1e3 - 22, 0), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=1)); ax4.text(-off * 1e3 - 24, -ld.STEM * 1e3 / 2, f"细杆\n{ld.STEM*1e3:.0f}", fontsize=8, color=INK2, ha="right", va="center")
ax4.annotate("", xy=(-off * 1e3 - 22, -(ld.STEM + ld.HEAD) * 1e3), xytext=(-off * 1e3 - 22, -ld.STEM * 1e3), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=1)); ax4.text(-off * 1e3 - 24, -(ld.STEM + ld.HEAD / 2) * 1e3, f"桨头\n{ld.HEAD*1e3:.0f}", fontsize=8, color=INK2, ha="right", va="center")
ax4.set_xlim(-118, 118); ax4.set_ylim(-100, 48); ax4.set_ylabel("沿桨轴 [mm]（从 E 向桨尖）"); ax4.set_title("④ 桨的正视图（沿桨轴看；脚蹼只在 ±y 张开）", loc="left", fontsize=10)
# ⑤ 侧视图：机构位姿
style(ax5, xlab="x [mm]（+x = 前进方向）"); ax5.set_aspect("equal"); ax5.set_ylabel("z [mm]（机体坐标，水线 z = −10）")
O2 = np.array(ld.O2); ax5.axhline(WL * 1e3, color=BLUE, lw=1.2); ax5.text(-70, WL * 1e3 + 1.5, "水线 WL（z = −10 mm）", color=BLUE, fontsize=9)
ax5.plot(tip[power, 0] * 1e3, tip[power, 1] * 1e3, color=PW, lw=1.2, ls="--"); ax5.plot(tip[~power, 0] * 1e3, tip[~power, 1] * 1e3, color=RC, lw=1.2, ls="--")
def draw_pose(tau0, col, lab, lw=1.6, alpha=1.0, linkage=True, ldx=0, ha="center", at_E=False):
    i = int(tau0 * N) % N; p = PP[i]
    C, Dd, Ee = p["C"], p["D"], p["E"]; B = p["B"]; u = p["u"]
    if linkage:
        ax5.plot([B[0] * 1e3, C[0] * 1e3], [B[1] * 1e3, C[1] * 1e3], color=INK2, lw=lw, alpha=alpha)                 # 摇臂 B–O2–C
        for P1, P2 in ((O2, Dd), (Dd, Ee), (C, Ee)): ax5.plot([P1[0] * 1e3, P2[0] * 1e3], [P1[1] * 1e3, P2[1] * 1e3], color=INK2, lw=lw * 0.8, alpha=alpha)  # 平行四边形
        for P_, nm, dx, dy in ((C, "C", 2, 2), (Dd, "D", -6, 2), (Ee, "E", -7, -1)): ax5.plot(*(P_ * 1e3), "o", color=INK2, ms=3, alpha=alpha); ax5.text(P_[0] * 1e3 + dx, P_[1] * 1e3 + dy, nm, fontsize=8, color=INK2, alpha=alpha)
    st = Ee + ld.STEM * u; tp_ = Ee + (ld.STEM + ld.HEAD) * u
    ax5.plot([Ee[0] * 1e3, st[0] * 1e3], [Ee[1] * 1e3, st[1] * 1e3], color=col, lw=lw, alpha=alpha)            # 细杆
    ax5.plot([st[0] * 1e3, tp_[0] * 1e3], [st[1] * 1e3, tp_[1] * 1e3], color=col, lw=lw * 3.2, alpha=alpha, solid_capstyle="butt")  # 桨头（侧视为板的棱）
    if at_E: ax5.text(Ee[0] * 1e3 + ldx, Ee[1] * 1e3 + 3, lab, color=col, fontsize=8.5, ha=ha, va="bottom", alpha=min(1, alpha + 0.3))
    else: ax5.text(tp_[0] * 1e3 + ldx, tp_[1] * 1e3 - 4, lab, color=col, fontsize=8.5, ha=ha, va="top", alpha=min(1, alpha + 0.3))
draw_pose(0.0, PW, "τ = 0 / 1（ψ = %.0f°，伸展）" % PSI0)
draw_pose(DUTY / 2, PW, "τ = %.2f" % (DUTY / 2), alpha=0.5, linkage=False, ldx=-12, ha="right")
draw_pose(DUTY - 1e-3, PW, "τ = DUTY（ψ = %.0f°）" % (PSI0 - SWEEP), linkage=False)
draw_pose(DUTY + SF * (1 - DUTY), RC, "τ = %.2f 蜷缩到位（φ = %.0f°，E 抬高 %.0f mm，桨头顶出水）" % (DUTY + SF * (1 - DUTY), PR, D["E_lift_mm"]), ldx=10, ha="left", at_E=True)
ax5.plot(*(O2 * 1e3), "o", color=INK, ms=5); ax5.text(O2[0] * 1e3 + 2, O2[1] * 1e3 + 2, "O2（髋轴，下舵机）", fontsize=9, color=INK)
O1 = np.array(ld.O1); ax5.plot(*(O1 * 1e3), "o", color=GRAY, ms=4); ax5.text(O1[0] * 1e3 + 2, O1[1] * 1e3 + 2, "O1（上舵机→曲柄1→摇臂）", fontsize=8, color=INK2)
# ψ 角标注（τ=0 位姿）
ax5.plot([O2[0] * 1e3, O2[0] * 1e3 + 30], [O2[1] * 1e3, O2[1] * 1e3], color=GRAY, lw=0.8, ls=":")
ax5.add_patch(Arc(O2 * 1e3, 40, 40, theta1=0, theta2=PSI0, color=RED, lw=1.2)); ax5.text(O2[0] * 1e3 + 17, O2[1] * 1e3 + 15, f"ψ = {PSI0:.0f}°", color=RED, fontsize=9)
ax5.set_title(f"⑤ 侧视图：{a.name} 的机构位姿（蓝 = 划水相，橙 = 回收相；虚线 = 桨尖轨迹）", loc="left", fontsize=10)
ax5.set_xlim(-75, 130); ax5.set_ylim(-108, 60)
fig.suptitle(f"BODY2 单腿步态参数图解 —— {a.name}：T={T} s  DUTY={DUTY}  SWEEP={SWEEP:.0f}°  PSI_START={PSI0:.0f}°  PHI_EXT={PE:.0f}°  PHI_RET={PR:.0f}°  SWITCH_FRAC={SF}  FAN_DEG={FAN:.0f}°  FAN_CLOSE={FC:.0f}  CLOSED_W={CW*1e3:.0f} mm", fontsize=11.5, x=0.01, ha="left")
fn = os.path.join(OUT, f"fig_gait_params_{a.name}.png"); plt.savefig(fn, dpi=120, facecolor=SURF, bbox_inches="tight"); print("→", fn)
