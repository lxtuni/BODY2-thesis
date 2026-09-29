#!/usr/bin/env python3
"""
结构无关的“目标轨迹”（target trajectory）：单腿在水中、以最高周期平均前进推力为目标的桨头运动 + 开合时机
    python3 optimal_trajectory.py [--out .] [--mp4]
约束（不是结构，是包络）：舵机速率上限 → 桨头中心峰值速度 VMAX；桨头 = 现有 120° 扇（42 mm）；可达范围 ≈ 现有（划水线深度 ZP、行程 LP）；
             水线 WL = −10 mm（绝对坐标，与 leg_dynamics.py 动画同一坐标系：O2 = (0.0634, 0.004)）。
规则（来自 03_最优行程理论与开合时机.md）：划水相直线向后、桨面垂直运动、梯形速度律；减速一开始就收扇；停顿；回收相收拢、抬高、羽化、同速率返回；
             回收末段减速一开始就开扇。
输出：optimal_trajectory.csv（每 2 ms 一行：绝对与相对 O2 的坐标、速度、桨轴角、开度、准定常示意力、相位、事件）
      optimal_trajectory_events.json、fig_optimal_trajectory_storyboard.png、anim_opt/f_%03d.png（--mp4 时再合成 optimal_trajectory_anim.mp4）
"""
import os, json, argparse, glob, subprocess
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=HERE); ap.add_argument("--mp4", action="store_true"); ap.add_argument("--no-anim", action="store_true", help="只写 CSV/JSON，不画分镜和动画帧"); ap.add_argument("--strip", action="store_true", help="另存 1×4 精简分镜 fig_target_traj_strip.png（报告用）")
ap.add_argument("--vmax", type=float, default=0.35, help="桨头中心峰值速度 [m/s]（舵机速率上限 ≈ 215 °/s × r ≈ 95 mm）")
ap.add_argument("--tramp", type=float, default=0.045, help="加/减速时间 [s]")
ap.add_argument("--lp", type=float, default=0.084, help="划水直线行程 [m]（≈ 2 个桨弦）")
ap.add_argument("--zp", type=float, default=-0.085, help="划水线：桨头中心绝对 z [m]（桨竖直时桨根 E 在水线附近、桨尖 −106 mm；现动画桨尖 ≈ −89 mm）")
ap.add_argument("--lift", type=float, default=0.030, help="回收相抬升 [m]")
ap.add_argument("--pause", type=float, default=0.030, help="划水结束后停顿 [s]（Kwak 2019 的 relaxation phase，待 CFD C3 验证）")
ap.add_argument("--tfan", type=float, default=0.060, help="脚蹼开/合过渡时间 [s]")
ap.add_argument("--feather", type=int, default=1, help="回收相收拢桨羽化（长轴顺流）1/0")
a = ap.parse_args()
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False

# ---------------- 几何 / 常数 ----------------
O2 = np.array([0.0634, 0.004]); WL = -0.010
HEAD, FAN_DEG, CLOSED_W, THK, STEM = 0.042, 120.0, 0.012, 0.0056, 0.050
RHO, CN_OPEN, CN_CLOSED, CD_EDGE = 998.8, 1.3, 1.2, 1.0
A_OPEN = 0.5 * HEAD ** 2 * np.radians(FAN_DEG)                      # 扇面积 18.5 cm²
A_CLOSED = CLOSED_W * HEAD                                           # 收拢正面 5.0 cm²
A_EDGE = THK * HEAD                                                  # 羽化后迎流 2.4 cm²
W_OPEN = 2 * HEAD * np.sin(np.radians(FAN_DEG) / 2)                  # 张开时扇尾宽 73 mm
XC = O2[0]                                                           # 行程居中在髋轴正下方
S = np.array([XC + a.lp / 2, a.zp]); E = np.array([XC - a.lp / 2, a.zp])   # 划水起点（前）/ 终点（后）
R1 = E + np.array([0.020, a.lift]); R2 = S + np.array([-0.020, a.lift])    # 回收：斜上 → 水平 → 斜下

def trapezoid_s(L, v, tr, t):
    """梯形速度律：0→v 用 tr，匀速，v→0 用 tr；返回 (弧长 s(t), 速度 |v|(t), 总时长)"""
    tau = L / v + tr
    t = np.asarray(t, float); s = np.where(t < tr, 0.5 * v / tr * t ** 2, np.where(t < tau - tr, 0.5 * v * tr + v * (t - tr), L - 0.5 * v / tr * np.clip(tau - t, 0, None) ** 2))
    vv = np.where(t < tr, v / tr * t, np.where(t < tau - tr, v, v / tr * np.clip(tau - t, 0, None)))
    return np.clip(s, 0, L), vv, tau

def polyline(pts, ds=0.0005, rc=0.008):
    """折线加密并把拐角倒圆（弧长 ±rc 内做移动平均），返回点列与累计弧长"""
    P = []
    for p, q in zip(pts[:-1], pts[1:]):
        n = max(2, int(np.linalg.norm(q - p) / ds)); P.append(np.linspace(p, q, n, endpoint=False))
    P = np.vstack(P + [pts[-1][None]]); k = max(1, int(rc / ds))
    Q = P.copy()
    for i in range(k, len(P) - k): Q[i] = P[i - k:i + k + 1].mean(0)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]); return Q, s

# ---------------- 时序 ----------------
tau_p = a.lp / a.vmax + a.tramp                                       # 划水相
Rpts, Rs = polyline(np.array([E, R1, R2, S])); LR = Rs[-1]; tau_r = LR / a.vmax + a.tramp   # 回收相
T = tau_p + a.pause + tau_r; DUTY = tau_p / T
t_close_cmd = tau_p - a.tramp; t_stop = tau_p; t_rec = tau_p + a.pause; t_open_cmd = T - a.tramp
dt = 0.002; t = np.arange(0, T, dt); tau = t / T

x = np.zeros_like(t); z = np.zeros_like(t); vmag = np.zeros_like(t); phase = np.empty(len(t), dtype=object); theta = np.zeros_like(t)
# 划水相：直线 S→E
m = t < tau_p; s_, v_, _ = trapezoid_s(a.lp, a.vmax, a.tramp, t[m]); x[m] = S[0] - s_; z[m] = S[1]; vmag[m] = v_; phase[m] = "power"; theta[m] = 0.0
# 停顿
m = (t >= tau_p) & (t < t_rec); x[m] = E[0]; z[m] = E[1]; vmag[m] = 0; phase[m] = "pause"; theta[m] = 0.0
# 回收相：沿倒圆折线 E→R1→R2→S
m = t >= t_rec; s_, v_, _ = trapezoid_s(LR, a.vmax, a.tramp, t[m] - t_rec); x[m] = np.interp(s_, Rs, Rpts[:, 0]); z[m] = np.interp(s_, Rs, Rpts[:, 1]); vmag[m] = v_; phase[m] = "recovery"
# 桨轴角 θ 始终 0（桨尖朝下、桨头竖直）；羽化 = 收拢后绕桨轴转 90°（像风向标），侧视图仍是竖直条，只是 x 向迎流从 12 mm 变成 5.6 mm
theta[:] = 0.0
# 速度矢量
vx = np.gradient(x, t); vz = np.gradient(z, t)
# 开度 f：收扇指令 = 划水减速开始；开扇指令 = 回收减速开始；过渡 tfan
fan = np.ones_like(t)
fan = np.where(t >= t_close_cmd, np.clip(1 - (t - t_close_cmd) / a.tfan, 0, 1), fan)
fan = np.where(t >= t_open_cmd, np.clip((t - t_open_cmd) / a.tfan, 0, 1), fan)
fan = np.where(t < (t_open_cmd + a.tfan - T), np.clip((t + T - t_open_cmd) / a.tfan, 0, 1), fan)   # 跨周期：开扇在下一周期开头完成
# 准定常示意力（U = 0）：法向力 ½ρC_N A (v·n)|v·n| + 边缘阻力；n = 桨面法向（垂直桨轴）
u = np.stack([-np.sin(np.radians(theta)), -np.cos(np.radians(theta))], 1)      # 桨轴单位向量（从根指向桨尖）
n = np.stack([-u[:, 1], u[:, 0]], 1)                                          # 桨面法向
vn = vx * n[:, 0] + vz * n[:, 1]
A_face = fan * A_OPEN + (1 - fan) * (A_EDGE if a.feather else A_CLOSED); CN = fan * CN_OPEN + (1 - fan) * CN_CLOSED   # 羽化：收拢后只剩边缘迎流
Fn = -0.5 * RHO * CN * A_face * vn * np.abs(vn)                                # 力与相对来流同向 = 与桨速反向
A_z = THK * (fan * W_OPEN + (1 - fan) * CLOSED_W)                             # z 向迎流：桨头底边（很小）
Fx = Fn * n[:, 0]; Fz = Fn * n[:, 1] - 0.5 * RHO * CD_EDGE * A_z * vz * np.abs(vz)
# 事件表
ev = {"T": T, "DUTY": DUTY, "f_Hz": 1 / T, "tau_p": tau_p, "tau_r": tau_r, "pause": a.pause, "L_power_mm": a.lp * 1e3, "L_recovery_mm": LR * 1e3, "vmax": a.vmax, "tramp": a.tramp, "tfan": a.tfan,
      "events_s": {"open_cmd": t_open_cmd, "power_start": 0.0, "fully_open": (t_open_cmd + a.tfan) % T, "close_cmd": t_close_cmd, "stop": t_stop, "fully_closed": t_close_cmd + a.tfan, "recovery_start": t_rec},
      "frame": "绝对坐标同 leg_dynamics.py 动画：O2 = (0.0634, 0.004) m，WL = −0.010 m；相对坐标 = 绝对 − O2（固件坐标：+x 前，+z 上）",
      "path_abs_m": {"S_power_start": S.tolist(), "E_power_end": E.tolist(), "R1": R1.tolist(), "R2": R2.tolist()},
      "quasi_steady_mean_mN": {"Fx_cycle": float(Fx.mean() * 1e3), "Fx_power": float(Fx[phase == "power"].mean() * 1e3), "Fx_recovery": float(Fx[phase == "recovery"].mean() * 1e3)}}
ev["events_tau"] = {k: v / T for k, v in ev["events_s"].items()}
json.dump(ev, open(os.path.join(a.out, "optimal_trajectory_events.json"), "w"), indent=1, ensure_ascii=False)
with open(os.path.join(a.out, "optimal_trajectory.csv"), "w") as fh:
    fh.write("t_s,tau,x_abs_m,z_abs_m,x_rel_O2_mm,z_rel_O2_mm,vx_m_s,vz_m_s,v_m_s,theta_deg,fan_open,A_face_cm2,Fx_qs_mN,Fz_qs_mN,phase\n")
    for i in range(len(t)):
        fh.write(f"{t[i]:.3f},{tau[i]:.4f},{x[i]:.5f},{z[i]:.5f},{(x[i]-O2[0])*1e3:.2f},{(z[i]-O2[1])*1e3:.2f},{vx[i]:.4f},{vz[i]:.4f},{vmag[i]:.4f},{theta[i]:.1f},{fan[i]:.3f},{A_face[i]*1e4:.2f},{Fx[i]*1e3:.1f},{Fz[i]*1e3:.1f},{phase[i]}\n")
print(f"T = {T:.3f} s ({1/T:.2f} Hz)  DUTY = {DUTY:.2f}  τ_p = {tau_p:.3f}  pause = {a.pause}  τ_r = {tau_r:.3f}  L_p = {a.lp*1e3:.0f} mm  L_r = {LR*1e3:.0f} mm")
for k, v in ev["events_s"].items(): print(f"  {k:15s} t = {v:6.3f} s  τ = {v/T:.3f}")
print(f"准定常示意：周期平均 Fx = {Fx.mean()*1e3:.1f} mN（划水相 {Fx[phase=='power'].mean()*1e3:.1f}，回收相 {Fx[phase=='recovery'].mean()*1e3:.1f}）")

# ---------------- 画帧 ----------------
RED, BLUE, GREY, ORANGE, PURPLE = "#c0392b", "#2980b9", "#7f8c8d", "#e67e22", "#8e44ad"
def draw(ax_main, ax_tl, i, small=False):
    ax = ax_main; pw = phase[i] == "power"; col = RED if pw else (GREY if phase[i] == "pause" else BLUE)
    # 轨迹环：划水红、回收蓝、停顿灰
    for ph, c in (("power", RED), ("recovery", BLUE)):
        mm = phase == ph; ax.plot(x[mm], z[mm], "-", color=c, lw=1.0, alpha=0.35)
    # 事件标记
    for key, lab, mk, off, ha in (("open_cmd", "开扇指令", "^", (7, 4), "left"), ("close_cmd", "收扇指令", "v", (0, -15), "center"), ("stop", "停顿·全闭", "s", (-7, -15), "right"), ("fully_open", "全开", "o", (7, -12), "left")):
        j = int(round(ev["events_s"][key] / dt)) % len(t); ax.plot(x[j], z[j], mk, color="k", ms=6, mfc="white", zorder=6)
        if not small: ax.annotate(lab, (x[j], z[j]), xytext=off, textcoords="offset points", fontsize=8, color="#333", ha=ha)
    # 髋轴与可达范围示意
    ax.plot(O2[0], O2[1], "k+", ms=9, mew=1.5); ax.add_patch(plt.Circle(O2, 0.105, fill=False, ls=":", lw=0.8, color="#999"))
    if not small: ax.text(O2[0] + 0.004, O2[1] + 0.004, "髋轴 O2\n（可达 r ≲ 105 mm 示意）", fontsize=8, color="#666")
    # 桨：桨头 42 mm（线宽示意开度）+ 细杆 50 mm 示意
    P = np.array([x[i], z[i]]); ui = u[i]; ni = n[i]; f = fan[i]
    tp = P + 0.5 * HEAD * ui; sh = P - 0.5 * HEAD * ui
    ax.plot(sh[0], sh[1], "o", color="#a63603", ms=4, alpha=0.6)                      # 桨根铰点
    wfront = f * W_OPEN + (1 - f) * (THK if a.feather else CLOSED_W); w = 0.0008 + 0.0042 * wfront / W_OPEN
    poly = np.array([sh - ni * w, tp - ni * w, tp + ni * w, sh + ni * w]); ax.add_patch(Polygon(poly, closed=True, color="#a63603", alpha=0.9, zorder=5))
    wtip = f * W_OPEN + (1 - f) * CLOSED_W
    st = f"脚蹼张开 宽 {wtip*1e3:.0f} mm ({f*100:.0f} %)" if f > 0.5 else (f"脚蹼收拢 {wtip*1e3:.0f} mm{'，羽化 迎流 ' + format(wfront*1e3, '.0f') + ' mm' if a.feather else ''}" if f < 0.02 else f"脚蹼{'张开中' if t[i] > T/2 or t[i] < 0.05 else '收拢中'} 宽 {wtip*1e3:.0f} mm ({f*100:.0f} %)")
    ax.text(tp[0] + 0.004, tp[1] - 0.004, st, color=col, fontsize=8 if small else 9)
    # 速度与力
    vs = 0.06; ax.arrow(P[0], P[1], vx[i] * vs, vz[i] * vs, width=0.0006, color="#555", length_includes_head=True, alpha=0.7)
    sc = 0.25; ax.arrow(P[0], P[1], Fx[i] * sc, Fz[i] * sc, width=0.0012, color=col, length_includes_head=True, zorder=7)
    lab = {"power": "划水相 (power)", "pause": "停顿 (pause)", "recovery": "回收相 (recovery)"}[phase[i]]
    ax.text(-0.016, 0.046, f"t/T = {tau[i]:.2f}   {lab}   v = {vmag[i]:.2f} m/s\nF_x = {Fx[i]*1e3:+.1f} mN   F_z = {Fz[i]*1e3:+.1f} mN（准定常示意）", fontsize=8 if small else 10, va="top", bbox=dict(fc="white", ec=col))
    ax.axhline(WL, color="tab:blue", ls="--", lw=1); ax.text(0.135, WL + 0.002, "WL", fontsize=8, color="tab:blue", ha="right")
    ax.set_xlim(-0.02, 0.14); ax.set_ylim(-0.12, 0.05); ax.set_aspect("equal"); ax.grid(alpha=0.25)
    ax.set_xlabel("x [m] → 前", fontsize=9); ax.set_ylabel("z [m]", fontsize=9); ax.tick_params(labelsize=8)
    if not small: ax.set_title(f"目标轨迹（结构无关，最高推力）· T = {T:.2f} s, DUTY = {DUTY:.2f}, v_max = {a.vmax} m/s", fontsize=10.5)
    # 时间线
    if ax_tl is not None:
        ax_tl.clear()
        ax_tl.axvspan(0, DUTY, color=RED, alpha=0.08, lw=0); ax_tl.axvspan(DUTY, t_rec / T, color=GREY, alpha=0.18, lw=0); ax_tl.axvspan(t_rec / T, 1, color=BLUE, alpha=0.08, lw=0)
        ax_tl.plot(tau, vmag / a.vmax * 0.9 + 2.1, color="#555", lw=1.5); ax_tl.plot(tau, fan * 0.9 + 1.05, color=ORANGE, lw=1.5)
        zz = (z - a.zp) / a.lift; ax_tl.plot(tau, zz * 0.9 + 0.0, color=PURPLE, lw=1.5)
        for yy, lb in ((2.55, "v / v_max"), (1.5, "开度"), (0.45, "抬升")): ax_tl.text(-0.01, yy, lb, ha="right", va="center", fontsize=7.5, color="#444")
        for key, lb in (("open_cmd", "开扇"), ("close_cmd", "收扇"), ("stop", "停"), ("recovery_start", "回收")):
            tv = ev["events_tau"][key]; ax_tl.axvline(tv, color="k", lw=0.7, ls=":"); ax_tl.text(tv, 3.15, lb, ha="center", fontsize=7.5, color="#222")
        ax_tl.axvline(tau[i], color=col, lw=2); ax_tl.set_xlim(-0.005, 1.0); ax_tl.set_ylim(-0.1, 3.5); ax_tl.set_yticks([]); ax_tl.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax_tl.tick_params(labelsize=7.5)
        ax_tl.set_xlabel("τ = t / T", fontsize=8, labelpad=1)
        for s_ in ("top", "right", "left"): ax_tl.spines[s_].set_visible(False)

# 分镜：12 帧
if a.no_anim: raise SystemExit(0)
if a.strip:
    fig, axs = plt.subplots(1, 4, figsize=(22, 6.2))
    for ax, tp_ in zip(axs, [0.08, 0.24, t_rec + 0.10, t_open_cmd + 0.01]):
        i = int(round((tp_ % T) / dt)) % len(t); draw(ax, None, i, small=True); ax.set_title(f"t = {t[i]:.3f} s   τ = {tau[i]:.2f}", fontsize=11)
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "fig_target_traj_strip.png"), dpi=80); plt.close(); print("→ fig_target_traj_strip.png")
picks = [0.0, 0.08, 0.16, t_close_cmd, t_stop + 0.5 * a.pause, t_rec + 0.02, t_rec + 0.10, t_rec + 0.18, t_rec + 0.26, t_open_cmd - 0.03, t_open_cmd + 0.01, T - 0.012]
fig, axs = plt.subplots(3, 4, figsize=(22, 17))
for ax, tp_ in zip(axs.ravel(), picks):
    i = int(round((tp_ % T) / dt)) % len(t); draw(ax, None, i, small=True); ax.set_title(f"t = {t[i]:.3f} s   τ = {tau[i]:.2f}", fontsize=10)
fig.suptitle(f"△ 开扇指令  ▽ 收扇指令  □ 停顿·全闭  ○ 全开\n结构无关的目标轨迹分镜：直线划水（梯形律，v_max {a.vmax} m/s）→ 减速即收扇 → 停顿 {a.pause*1e3:.0f} ms → 收拢抬升 {a.lift*1e3:.0f} mm 羽化返回 → 减速即开扇   |   T = {T:.2f} s, f = {1/T:.2f} Hz, DUTY = {DUTY:.2f}", fontsize=13)
plt.tight_layout(rect=(0, 0, 1, 0.955)); plt.savefig(os.path.join(a.out, "fig_optimal_trajectory_storyboard.png"), dpi=70); plt.close()
print("→ fig_optimal_trajectory_storyboard.png")

# 动画帧（16 fps、一个周期放慢到 5 s = 80 帧）
fr_dir = os.path.join(a.out, "anim_opt"); os.makedirs(fr_dir, exist_ok=True)
for old in glob.glob(os.path.join(fr_dir, "f_*.png")): os.remove(old)
nfr = 80
for j in range(nfr):
    i = int(round(j / nfr * len(t))) % len(t)
    fig = plt.figure(figsize=(7, 8.2)); gs = fig.add_gridspec(2, 1, height_ratios=[6.8, 1.4], hspace=0.28, left=0.1, right=0.98, top=0.95, bottom=0.07)
    ax = fig.add_subplot(gs[0]); ax_tl = fig.add_subplot(gs[1]); draw(ax, ax_tl, i)
    plt.savefig(os.path.join(fr_dir, f"f_{j:03d}.png"), dpi=80); plt.close()
print(f"动画帧 {nfr} 张 → {fr_dir}/")
if a.mp4:
    out = os.path.join(a.out, "optimal_trajectory_anim.mp4")
    subprocess.run(["ffmpeg", "-y", "-v", "error", "-framerate", "16", "-i", os.path.join(fr_dir, "f_%03d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", out], check=True); print("→", out)
