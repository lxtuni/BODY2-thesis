#!/usr/bin/env python3
"""
单腿动态 CFD —— 三向力 (x/y/z)、力矩、阻力预算、四腿对角步态叠加、推力–阻力平衡 估算
    python3 analyze_leg_xyz.py <results_dir> [out_dir]
读 results/force_{open,closed}.dat、moment_*.dat、gait_*.json、linkage.json（analyze_leg.py 之后跑），
写 xyz_data.json 与 fig_xyz.png / fig_moments.png / fig_drag_budget.png / fig_gait4.png / fig_balance.png。
坐标：x 前进 (+x = 推力)，y 左，z 上；力矩参考点 = O2（髋轴，单腿算例）；整机叠加时换算到船体中点、水线高度。
"""
import os, sys, json, importlib.util
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

def _cjk():
    import glob as g
    c = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]
    c += g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + g.glob("/mnt/c/Windows/Fonts/simhei.ttf")
    for f in c:
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return
        except Exception: pass
_cjk(); plt.rcParams["axes.unicode_minus"] = False
BLUE, RED, GRAY, INK, INK2, SURF, ORANGE, GREEN, PURPLE = "#2a78d6", "#e34948", "#a8a7a1", "#0b0b0b", "#52514e", "#fcfcfb", "#eb6834", "#3a9d5d", "#7b5cd6"
GRID = "#e6e5e1"

R = os.path.abspath(sys.argv[1]); OUT = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else R
os.makedirs(OUT, exist_ok=True); HERE = os.path.dirname(os.path.abspath(__file__))
RHO = 998.8

# ---------------- 读数据 ----------------
def load(tag):
    fn = os.path.join(R, f"force_{tag}.dat")
    if not os.path.exists(fn): return None
    d = np.loadtxt(fn, comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
    m = np.loadtxt(os.path.join(R, f"moment_{tag}.dat"), comments="#"); m = m[np.argsort(m[:, 0])]; _, i = np.unique(m[:, 0], return_index=True); m = m[i]
    g = json.load(open(os.path.join(R, f"gait_{tag}.json")))
    return dict(t=d[:, 0], F=d[:, 1:4], P=d[:, 4:7], V=d[:, 7:10], M=m[:, 1:4], MP=m[:, 4:7], MV=m[:, 7:10], gait=g)
cases = {k: v for k in ("open", "closed") if (v := load(k)) is not None}
assert "open" in cases and "closed" in cases, "需要 open 与 closed 两个算例"
g = cases["open"]["gait"]; T = g["T"]; DUTY = g["DUTY"]; U0 = g.get("U", 0.0)
FC = g.get("FAN_CLOSE", -1); tc = DUTY if (FC is None or FC < 0) else FC
ncyc = int(np.floor(min(c["t"][-1] for c in cases.values()) / T + 1e-6)); t0, t1 = (ncyc - 1) * T, ncyc * T
N = 2500; tt = np.linspace(t0, t1, N, endpoint=False); tau = (tt - t0) / T; dt = T / N
def smooth(y, w=0.02):
    k = max(1, int(round(w / dt))); ker = np.ones(k) / k
    yp = np.concatenate([y[-k:], y, y[:k]]); return np.convolve(yp, ker, "same")[k:-k]
def interp_cols(c, key):
    return np.stack([np.interp(tt, c["t"], c[key][:, i]) for i in range(3)], 1)
res = {}
for tag, c in cases.items():
    res[tag] = {k: interp_cols(c, k) for k in ("F", "P", "V", "M")}
fan_open = tau < tc
comp = {k: np.where(fan_open[:, None], res["open"][k], res["closed"][k]) for k in res["open"]}
res["composite"] = comp
tr = g.get("SWITCH_FRAC", 0.3) * (1 - DUTY)
phases = [("划水相", 0.0, DUTY), ("收拢切换", DUTY, DUTY + tr), ("蜷缩前扫", DUTY + tr, 1 - tr), ("伸展切换", 1 - tr, 1.0)]
def phase_mean_contrib(y):     # 各相位冲量 / T  [mN]
    return [float(np.trapezoid(y[(tau >= a) & (tau < b)], tt[(tau >= a) & (tau < b)]) / T * 1e3) for _, a, b in phases]

# ---------------- 运动学：E(t)、ψ(t) → 功率 ----------------
ld = None
try:
    for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
        if k in g: os.environ[k] = str(g[k])
    if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
    spec = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(spec); spec.loader.exec_module(ld)
except Exception as e:
    print("运动学加载失败:", e)
O2 = np.array([0.0634, 0.0040]) if ld is None else np.array(ld.O2)
power = None
if ld is not None:
    psi, phi, _ = ld.gait(tt); P = [ld.pose(p, q) for p, q in zip(psi, phi)]
    E = np.array([p["E"] for p in P])                                 # (x, z) 销孔 E
    vE = np.gradient(E, tt, axis=0); omega_y = -np.gradient(np.unwrap(psi), tt)   # 绕 +y 的角速度（x→z 旋转为负）
    rEx, rEz = E[:, 0] - O2[0], E[:, 1] - O2[1]
    power = {}
    for tag in ("open", "closed", "composite"):
        F, M = res[tag]["F"], res[tag]["M"]
        M_E_y = M[:, 1] - (rEz * F[:, 0] - rEx * F[:, 2])              # 关于 E 的力矩 = M_O2 − (r_E−r_O2) × F
        p_fluid = F[:, 0] * vE[:, 0] + F[:, 2] * vE[:, 1] + M_E_y * omega_y   # 流体对桨做功功率
        power[tag] = -p_fluid                                           # 舵机输入给水的功率 [W]

# ---------------- 三向力、力矩统计 ----------------
def stats(y):
    ys = smooth(y)
    return dict(mean=float(y.mean() * 1e3), rms=float(np.sqrt(((y - y.mean()) ** 2).mean()) * 1e3), max_s=float(ys.max() * 1e3), min_s=float(ys.min() * 1e3),
                max_raw=float(y.max() * 1e3), min_raw=float(y.min() * 1e3), contrib=phase_mean_contrib(y))
data = dict(T=T, DUTY=DUTY, U=U0, tc=tc, ncycles=ncyc, phases=[p[0] for p in phases], phase_ranges=[[a, b] for _, a, b in phases], cases={})
for tag, r in res.items():
    d = {}
    for i, ax in enumerate("xyz"):
        d[f"F{ax}"] = stats(r["F"][:, i]); d[f"F{ax}"]["press_mean"] = float(r["P"][:, i].mean() * 1e3); d[f"F{ax}"]["visc_mean"] = float(r["V"][:, i].mean() * 1e3)
        d[f"F{ax}"]["visc_absmean"] = float(np.abs(r["V"][:, i]).mean() * 1e3); d[f"F{ax}"]["press_absmean"] = float(np.abs(r["P"][:, i]).mean() * 1e3)
        d[f"M{ax}"] = stats(r["M"][:, i])                               # [mN·m]
    Fx = r["F"][:, 0]
    d["impulse_pos"] = float(np.trapezoid(np.maximum(Fx, 0), tt) * 1e3); d["impulse_neg"] = float(np.trapezoid(np.minimum(Fx, 0), tt) * 1e3)   # mN·s
    d["impulse_pos_phase"] = phase_mean_contrib(np.maximum(Fx, 0)); d["impulse_neg_phase"] = phase_mean_contrib(np.minimum(Fx, 0))
    Fmag = np.linalg.norm(r["F"], axis=1); d["Fmag_max_s"] = float(smooth(Fmag).max() * 1e3)
    # 合力方向：划水相内 F 与 +x 的夹角（平滑后、|F| 大于峰值 20% 的时刻）
    Fs = np.stack([smooth(r["F"][:, i]) for i in range(3)], 1); k = (tau < DUTY) & (np.linalg.norm(Fs, axis=1) > 0.2 * np.linalg.norm(Fs[tau < DUTY], axis=1).max())
    d["power_angle_deg"] = float(np.degrees(np.arctan2(Fs[k, 2], Fs[k, 0])).mean()) if k.any() else None
    if power is not None:
        pw = power[tag]; d["power_mean_mW"] = float(pw.mean() * 1e3); d["power_peak_mW"] = float(smooth(pw).max() * 1e3); d["power_min_mW"] = float(smooth(pw).min() * 1e3)
        d["power_phase_mW"] = phase_mean_contrib(pw); d["power_power_phase_mean_mW"] = float(pw[tau < DUTY].mean() * 1e3)
        d["thrust_per_power_mN_per_mW"] = float(Fx.mean() / pw.mean()) if pw.mean() > 0 else None
    data["cases"][tag] = d

# ---------------- 连杆（条带法）----------------
LK = json.load(open(os.path.join(R, "linkage.json"))) if os.path.exists(os.path.join(R, "linkage.json")) else None
lk_mean = LK["total"]["mean_Fx"] if LK else 0.0; lk_peak = LK["total"]["peak"] if LK else 0.0
lk_meanFz = LK["total"].get("mean_Fz", 0.0) if LK else 0.0
data["linkage"] = dict(mean_Fx=lk_mean, mean_Fz=lk_meanFz, peak=lk_peak, links={k: dict(mean_Fx=v["mean_Fx"], peak=v["peak"]) for k, v in LK["links"].items()}) if LK else None

# ---------------- 阻力预算（单腿，U = 0）----------------
c = data["cases"]["composite"]
budget = dict(thrust_impulse=c["impulse_pos"], drag_impulse=c["impulse_neg"], net=c["impulse_pos"] + c["impulse_neg"],
              drag_share=float(-c["impulse_neg"] / c["impulse_pos"]),
              thrust_phase=c["impulse_pos_phase"], drag_phase=c["impulse_neg_phase"],
              linkage_mean=lk_mean, viscous_mean=c["Fx"]["visc_mean"], pressure_mean=c["Fx"]["press_mean"],
              net_with_linkage=c["Fx"]["mean"] + lk_mean,
              closed_vs_open_switch=dict(open=data["cases"]["open"]["Fx"]["contrib"][1], closed=data["cases"]["closed"]["Fx"]["contrib"][1]),
              recovery_mean_force=float(sum(c["Fx"]["contrib"][1:]) / (1 - DUTY)), power_mean_force=float(c["Fx"]["contrib"][0] / DUTY))
data["budget"] = budget

# ---------------- 四腿对角步态叠加（忽略腿间干扰与船体）----------------
LEGS = {"FL": dict(phase=0.0, dx=0.0, y=+0.0688), "FR": dict(phase=0.5, dx=0.0, y=-0.0688), "RL": dict(phase=0.5, dx=-0.1435, y=+0.070), "RR": dict(phase=0.0, dx=-0.1435, y=-0.070)}
GAITS = {"diag": {"FL": 0.0, "FR": 0.5, "RL": 0.5, "RR": 0.0}, "sync": {"FL": 0.0, "FR": 0.0, "RL": 0.0, "RR": 0.0}, "bound": {"FL": 0.0, "FR": 0.0, "RL": 0.5, "RR": 0.5}, "walk": {"FL": 0.0, "RR": 0.25, "FR": 0.5, "RL": 0.75}}
CofR = np.array([0.005, 0.0, -0.010])                                   # 船体中点、水线（与整机算例一致）
def shift(y, ph):                                                       # 按相位平移一个周期的信号（周期性）
    return np.roll(y, -int(round(ph * N)), axis=0)
gait4 = {}
for gname, ph in GAITS.items():
    Ft = np.zeros((N, 3)); Mt = np.zeros((N, 3))
    for L, l in LEGS.items():
        sgn = 1.0 if l["y"] > 0 else -1.0                               # 右腿 = 左腿镜像：Fy、Mx、Mz 反号
        F = shift(comp["F"], ph[L]) * np.array([1, sgn, 1]); M = shift(comp["M"], ph[L]) * np.array([sgn, 1, sgn])
        rO2 = np.array([O2[0] + l["dx"], l["y"], O2[1]]) - CofR
        Ft += F; Mt += M + np.cross(np.broadcast_to(rO2, F.shape), F)
    gait4[gname] = dict(F=Ft, M=Mt)
data["gait4"] = {}
for gname, r in gait4.items():
    d = {}
    for i, ax in enumerate("xyz"):
        d[f"F{ax}"] = dict(mean=float(r["F"][:, i].mean() * 1e3), max_s=float(smooth(r["F"][:, i]).max() * 1e3), min_s=float(smooth(r["F"][:, i]).min() * 1e3), rms=float(np.std(r["F"][:, i]) * 1e3))
        d[f"M{ax}"] = dict(mean=float(r["M"][:, i].mean() * 1e3), max_s=float(smooth(r["M"][:, i]).max() * 1e3), min_s=float(smooth(r["M"][:, i]).min() * 1e3), rms=float(np.std(r["M"][:, i]) * 1e3))
    data["gait4"][gname] = d
data["gait4_legs"] = LEGS; data["gait4_CofR"] = CofR.tolist()

# ---------------- 推力–阻力平衡（估算）----------------
# 阻力：静止浮态报告 U = 0.2 m/s 全机（四腿收拢在水下）77.6 mN，其中船体本身按迎流面积比例（23.3 / 83.9 cm²，同 C_D）+ 层流摩擦估 ≈ 25 mN
R02_all, R02_hull, n_exp = 77.6, 25.0, 1.95
Ugrid = np.linspace(0, 0.35, 141)
R_all = R02_all * (Ugrid / 0.2) ** n_exp; R_hull = R02_hull * (Ugrid / 0.2) ** n_exp
# 推力随航速：两相"滑移"缩放  T(U) = T_pow (1 − U/v_p)² − D_rec (1 + U/v_r)² + L(U)
T_pow = c["Fx"]["contrib"][0]; D_rec = -float(sum(c["Fx"]["contrib"][1:]))
v_p = v_r = None
if ld is not None:
    head = E + (ld.STEM + ld.HEAD / 2) * np.array([p_["u"] for p_ in P])          # 桨头中心
    vh = np.linalg.norm(np.gradient(head, tt, axis=0), axis=1); Fx = comp["F"][:, 0]
    kp = (tau < DUTY) & (Fx > 0); kr = (tau >= DUTY) & (Fx < 0)
    v_p = float((vh[kp] * Fx[kp]).sum() / Fx[kp].sum()); v_r = float((vh[kr] * -Fx[kr]).sum() / (-Fx[kr]).sum())
v_p = v_p or 0.30; v_r = v_r or 0.25
v_tip_pow_max = float(np.linalg.norm(np.gradient(E + (ld.STEM + ld.HEAD) * np.array([p_["u"] for p_ in P]), tt, axis=0), axis=1)[tau < DUTY].max()) if ld is not None else None
lk_U = None
try:
    import importlib.util as _iu
    sp = _iu.spec_from_file_location("lkd", os.path.join(HERE, "linkage_drag.py")); lkd = _iu.module_from_spec(sp); sp.loader.exec_module(lkd)
    tl = np.linspace(0, T, 400, endpoint=False)
    lk_U = np.array([float(np.mean(lkd.leg_forces(ld, tl, 0.0, 0.0, -0.010, float(u))[0])) * 1e3 for u in Ugrid])   # mN，含 U 的相对速度
except Exception as e:
    print("连杆随 U 估算失败:", e); lk_U = np.full_like(Ugrid, lk_mean)
def T_leg(U, fp=1.0):
    return T_pow * np.clip(1 - U / (v_p * fp), 0, None) ** 2 - D_rec * (1 + U / (v_r * fp)) ** 2 + np.interp(U, Ugrid, lk_U)
T4 = 4 * T_leg(Ugrid); T4_lo = 4 * T_leg(Ugrid, 0.8); T4_hi = 4 * T_leg(Ugrid, 1.25); T4_const = np.full_like(Ugrid, 4 * c["Fx"]["mean"])
def cross(Tc, Rc):
    dlt = Tc - Rc; k = np.where((dlt[:-1] > 0) & (dlt[1:] <= 0))[0]
    if len(k) == 0: return None
    i = k[0]; return float(Ugrid[i] + (Ugrid[i + 1] - Ugrid[i]) * dlt[i] / (dlt[i] - dlt[i + 1]))
data["balance"] = dict(R02_all=R02_all, R02_hull=R02_hull, n_exp=n_exp, T0_4legs=float(4 * c["Fx"]["mean"]), T_pow=T_pow, D_rec=D_rec, v_p=v_p, v_r=v_r, v_tip_pow_max=v_tip_pow_max,
                       Ueq_center=cross(T4, R_hull), Ueq_lo=cross(T4_lo, R_all), Ueq_hi=cross(T4_hi, R_hull), Ueq_const_hull=cross(T4_const, R_hull), Ueq_const_all=cross(T4_const, R_all),
                       table=[dict(U=float(u), R_hull=float(R02_hull * (u / 0.2) ** n_exp), R_all=float(R02_all * (u / 0.2) ** n_exp), T4=float(4 * T_leg(u)), T4_lo=float(4 * T_leg(u, 0.8)), T4_hi=float(4 * T_leg(u, 1.25)), L1=float(np.interp(u, Ugrid, lk_U))) for u in (0.0, 0.05, 0.1, 0.15, 0.2)],
                       P_hull_02_mW=R02_hull * 0.2, P_all_02_mW=R02_all * 0.2, power_4legs_mW=4 * c.get("power_mean_mW", 0.0))
json.dump(data, open(os.path.join(OUT, "xyz_data.json"), "w"), indent=1, ensure_ascii=False)

# ================= 图 =================
cols = {"open": ORANGE, "closed": BLUE, "composite": INK}; names = {"open": "全程张开 open", "closed": "全程收拢 closed", "composite": "两态合成 composite"}
def phase_marks(ax, labels=True):
    for name, a, b in phases:
        ax.axvline(a, color=INK2, lw=0.6, alpha=0.35)
        if labels: ax.text((a + b) / 2, 0.98, name, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8.5, color=INK2)
    ax.axvspan(0, tc, color=ORANGE, alpha=0.05); ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_facecolor(SURF)

# 1) 三向力
fig, axs = plt.subplots(3, 1, figsize=(12, 10.5), facecolor=SURF, sharex=True)
lab = {"x": "F_x [mN]  (+ = 推力)", "y": "F_y [mN]  (+ = 左)", "z": "F_z [mN]  (+ = 上)"}
for i, ax in enumerate(axs):
    for tag in ("open", "closed", "composite"):
        r = res[tag]; y = r["F"][:, i]
        ax.plot(tau, y * 1e3, color=cols[tag], lw=0.5, alpha=0.18)
        ax.plot(tau, smooth(y) * 1e3, color=cols[tag], lw=2.4 if tag == "composite" else 1.4, label=f"{names[tag]}（均值 {data['cases'][tag]['F'+'xyz'[i]]['mean']:+.1f}）")
    phase_marks(ax, labels=(i == 0)); ax.set_ylabel(lab["xyz"[i]], color=INK); ax.legend(frameon=False, fontsize=8.5, loc="upper left" if i == 2 else "lower right", bbox_to_anchor=(0.0, 0.92) if i == 2 else None)
    if i == 1: ax.set_ylim(-3, 3)
axs[0].set_title(f"单腿三向力时间历程（第 {ncyc} 周期，20 ms 平滑，淡线为原始）  U = {U0} m/s", loc="left", fontsize=11)
axs[1].text(0.01, 0.9, f"单腿算例桨面在 y = 0 对称面内运动，F_y 只有 overset 数值噪声（原始 |F_y| < {max(abs(data['cases']['open']['Fy']['max_raw']),abs(data['cases']['open']['Fy']['min_raw'])):.1f} mN，平滑后 < 0.1 mN）；整机的侧向力与偏航来自左右腿相位差和船体不对称绕流（§ 四腿叠加）", transform=axs[1].transAxes, fontsize=8.5, color=INK2, va="top")
axs[2].set_xlabel("t / T", color=INK)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_xyz.png"), dpi=110, facecolor=SURF); plt.close()

# 2) 力矩（关于 O2）+ 功率
fig, axs = plt.subplots(1, 2, figsize=(14, 5), facecolor=SURF)
ax = axs[0]
for i, (lb, col) in enumerate((("M_x（横滚）", GREEN), ("M_y（俯仰，绕髋轴）", PURPLE), ("M_z（偏航）", RED))):
    y = comp["M"][:, i]; ax.plot(tau, y * 1e3, color=col, lw=0.5, alpha=0.18); ax.plot(tau, smooth(y) * 1e3, color=col, lw=2, label=f"{lb}  均值 {y.mean()*1e3:+.2f}，峰 {np.abs(smooth(y)).max()*1e3:.1f}")
phase_marks(ax); ax.set_ylabel("关于髋轴 O2 的力矩 [mN·m]", color=INK); ax.set_xlabel("t / T", color=INK); ax.legend(frameon=False, fontsize=8.5, loc="lower right")
ax.set_title("两态合成：桨对髋轴的力矩（舵机负载）", loc="left", fontsize=11)
ax = axs[1]
if power is not None:
    for tag in ("open", "closed", "composite"):
        y = power[tag] * 1e3; ax.plot(tau, y, color=cols[tag], lw=0.5, alpha=0.18); ax.plot(tau, smooth(y / 1e3) * 1e3, color=cols[tag], lw=2.4 if tag == "composite" else 1.4, label=f"{names[tag]}  均值 {y.mean():.1f} mW")
    phase_marks(ax); ax.set_ylabel("桨传给水的机械功率 [mW]", color=INK); ax.set_xlabel("t / T", color=INK); ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 0.9))
    ax.set_title("水动力功率 P = −(F·v_E + M_E·ω)（正 = 舵机做功）", loc="left", fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_moments.png"), dpi=110, facecolor=SURF); plt.close()

# 3) 阻力预算：分相位的推力冲量 / 阻力冲量 + 连杆 + 净值
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.25, 1]))
ax = axs[0]; ax.set_facecolor(SURF)
x = np.arange(len(phases)); w = 0.38
tp, dp = budget["thrust_phase"], budget["drag_phase"]
ax.bar(x - w / 2, tp, w, color=GREEN, label="推力冲量（F_x > 0 部分）/ T")
ax.bar(x + w / 2, dp, w, color=RED, label="阻力冲量（F_x < 0 部分）/ T")
for xi, v in zip(x - w / 2, tp): ax.text(xi, v + 0.5, f"{v:+.1f}", ha="center", fontsize=8.5, color=INK)
for xi, v in zip(x + w / 2, dp): ax.text(xi, v - 0.5, f"{v:+.1f}", ha="center", va="top", fontsize=8.5, color=INK)
ax.set_xticks(x); ax.set_xticklabels([f"{p[0]}\nτ {p[1]:.2f}–{p[2]:.2f}" for p in phases], fontsize=9); ax.axhline(0, color=INK, lw=0.8)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color=GRID, lw=0.8); ax.set_ylabel("折算到周期平均 [mN]", color=INK); ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.set_title(f"两态合成：正负冲量分相位   推力 {budget['thrust_impulse']/T:+.1f}  阻力 {budget['drag_impulse']/T:+.1f}  → 净 {budget['net']/T:+.1f} mN（阻力吃掉 {100*budget['drag_share']:.0f}%）", loc="left", fontsize=10.5)
ax = axs[1]; ax.set_facecolor(SURF)
items = [("划水相净推力", c["Fx"]["contrib"][0], GREEN), ("收拢切换", c["Fx"]["contrib"][1], RED), ("蜷缩前扫", c["Fx"]["contrib"][2], RED), ("伸展切换", c["Fx"]["contrib"][3], GREEN if c["Fx"]["contrib"][3] > 0 else RED),
         ("连杆（条带法）", lk_mean, RED)]
run = 0; xs = []
for i, (nm, v, col) in enumerate(items):
    ax.bar(i, v, 0.6, bottom=run, color=col, edgecolor=SURF); ax.text(i, run + v + (0.6 if v >= 0 else -0.6), f"{v:+.1f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8.5, color=INK)
    run += v; xs.append(i)
ax.bar(len(items), run, 0.6, color=INK); ax.text(len(items), run + 0.6, f"{run:+.1f}", ha="center", fontsize=9, color=INK, fontweight="bold")
ax.set_ylim(min(0, run) - 2, max(items[0][1], run) * 1.15)
ax.set_xticks(xs + [len(items)]); ax.set_xticklabels([it[0] for j, it in enumerate(items) if j in xs] + ["单腿净推力\n(含连杆)"], fontsize=8.5, rotation=15)
ax.axhline(0, color=INK, lw=0.8)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color=GRID, lw=0.8); ax.set_ylabel("周期平均 F_x [mN]", color=INK)
ax.set_title(f"阻力预算瀑布图（U = {U0} m/s）：划水相产生，回收相与连杆消耗", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_drag_budget.png"), dpi=110, facecolor=SURF); plt.close()

# 4) 四腿步态叠加
fig, axs = plt.subplots(2, 3, figsize=(15, 7.8), facecolor=SURF)
r = gait4["diag"]; rs = gait4["sync"]; rw = gait4["walk"]
def three(ax, i, key, ylab, unit_scale=1e3):
    ax.plot(tau, smooth(rs[key][:, i]) * unit_scale, color=GRAY, lw=1.2, label=f"同步 sync（峰 {np.abs(smooth(rs[key][:, i])).max()*unit_scale:.0f}）")
    ax.plot(tau, smooth(rw[key][:, i]) * unit_scale, color=PURPLE, lw=1.1, ls="--", label=f"四拍 walk（峰 {np.abs(smooth(rw[key][:, i])).max()*unit_scale:.1f}）")
    ax.plot(tau, smooth(r[key][:, i]) * unit_scale, color=INK, lw=2.3, label=f"对角 diag（均值 {r[key][:, i].mean()*unit_scale:+.1f}，峰 {np.abs(smooth(r[key][:, i])).max()*unit_scale:.1f}）")
    phase_marks(ax, labels=False); ax.set_ylabel(ylab, color=INK); ax.legend(frameon=False, fontsize=8.3, loc="lower right")
three(axs[0][0], 0, "F", "整机 F_x [mN]（+ = 推力）")
# 中上：对角步态两对腿各自的 F_x
ax = axs[0][1]
pairA = shift(comp["F"][:, 0], 0.0) * 2; pairB = shift(comp["F"][:, 0], 0.5) * 2
ax.plot(tau, smooth(pairA) * 1e3, color=ORANGE, lw=1.8, label="FL + RR（相位 0）")
ax.plot(tau, smooth(pairB) * 1e3, color=BLUE, lw=1.8, label="FR + RL（相位 0.5）")
ax.plot(tau, smooth(r["F"][:, 0]) * 1e3, color=INK, lw=2.3, label="合计 diag")
phase_marks(ax, labels=False); ax.set_ylabel("对角步态：两对腿的 F_x [mN]", color=INK); ax.legend(frameon=False, fontsize=8.3, loc="lower right")
three(axs[0][2], 2, "F", "整机 F_z [mN]（+ = 上）")
three(axs[1][0], 0, "M", "横滚 M_x [mN·m]")
three(axs[1][1], 1, "M", "俯仰 M_y [mN·m]（− = 抬头）")
three(axs[1][2], 2, "M", "偏航 M_z [mN·m]")
for ax in axs[1]: ax.set_xlabel("t / T（以 FL 相位计）", color=INK)
axs[1][0].set_ylim(-25, 25); axs[1][2].set_ylim(-25, 25)
axs[0][0].set_title("四腿叠加（单腿结果按相位、位置平移，右腿镜像；忽略腿间干扰与船体；F_y 左右抵消 ≡ 0）", loc="left", fontsize=10.5)
axs[1][0].set_title(f"力矩参考点 = 船体中点、水线 ({CofR[0]:.3f}, 0, {CofR[2]:.3f})；对角步态横滚/偏航几乎为零，四拍步态则明显", loc="left", fontsize=9.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_gait4.png"), dpi=110, facecolor=SURF); plt.close()

# 5) 推力–阻力平衡
fig, ax = plt.subplots(figsize=(9, 5.4), facecolor=SURF); ax.set_facecolor(SURF)
ax.fill_between(Ugrid, T4_lo, T4_hi, color=GREEN, alpha=0.15, lw=0)
ax.plot(Ugrid, T4, color=GREEN, lw=2.4, label=f"四腿推力 T(U)：两相滑移缩放（v_p = {v_p:.2f}，v_r = {v_r:.2f} m/s；带 = ±25%）")
ax.plot(Ugrid, T4_const, color=GRAY, lw=1.2, ls="--", label=f"四腿推力上界：系泊值 {4*c['Fx']['mean']:.0f} mN 不衰减")
ax.plot(Ugrid, R_hull, color=RED, lw=2.4, label=f"船体阻力 R_hull(U) ≈ {R02_hull:.0f} mN·(U/0.2)^{n_exp}（迎流面积比例估算）")
ax.plot(Ugrid, R_all, color=RED, lw=1.2, ls="--", label=f"全机静止浮态阻力 {R02_all} mN·(U/0.2)^{n_exp}（四腿收拢拖行，偏保守）")
b_ = data["balance"]
for ue, col, nm, yy in ((b_["Ueq_center"], GREEN, "估计", 0.62), (b_["Ueq_lo"], RED, "下界", 0.50), (b_["Ueq_hi"], GREEN, "上界", 0.74)):
    if ue: ax.axvline(ue, color=col, lw=0.8, ls=":"); ax.text(ue + 0.003, yy, f"U_eq {nm} ≈ {ue:.2f} m/s", color=col, fontsize=9, transform=ax.get_xaxis_transform())
ax.set_xlabel("航速 U [m/s]", color=INK); ax.set_ylabel("力 [mN]", color=INK); ax.set_xlim(0, 0.25); ax.set_ylim(0, 120)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.grid(color=GRID, lw=0.8); ax.legend(frameon=False, fontsize=8.3, loc="upper right")
ax.set_title("推力–阻力平衡（估算；整机 CFD 会给出真值）", loc="left", fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_balance.png"), dpi=110, facecolor=SURF); plt.close()

# ---------------- 摘要 ----------------
print(f"第 {ncyc} 周期，T = {T} s，U = {U0} m/s")
for tag in ("open", "closed", "composite"):
    d = data["cases"][tag]
    print(f"[{tag:<9}] Fx {d['Fx']['mean']:+6.1f}  Fy {d['Fy']['mean']:+5.2f}  Fz {d['Fz']['mean']:+6.1f} mN | Fx 平滑峰/谷 {d['Fx']['max_s']:.0f}/{d['Fx']['min_s']:.0f} | Fz 峰/谷 {d['Fz']['max_s']:.0f}/{d['Fz']['min_s']:.0f} | My 均 {d['My']['mean']:+.2f} 峰 {max(abs(d['My']['max_s']),abs(d['My']['min_s'])):.1f} mN·m | 功率 {d.get('power_mean_mW', float('nan')):.1f} mW | 推力冲量 {d['impulse_pos']/T:+.1f} 阻力冲量 {d['impulse_neg']/T:+.1f}")
print("阻力预算:", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in budget.items() if not isinstance(v, (list, dict))})
for gname, d in data["gait4"].items():
    print(f"[{gname:<5}] 4 腿 Fx {d['Fx']['mean']:+.1f} (峰 {d['Fx']['max_s']:.0f}/{d['Fx']['min_s']:.0f})  Fy 峰 {max(abs(d['Fy']['max_s']),abs(d['Fy']['min_s'])):.1f}  Fz {d['Fz']['mean']:+.1f} (峰 {d['Fz']['max_s']:.0f}/{d['Fz']['min_s']:.0f})  Mx 峰 {max(abs(d['Mx']['max_s']),abs(d['Mx']['min_s'])):.2f}  My {d['My']['mean']:+.2f} (峰 {max(abs(d['My']['max_s']),abs(d['My']['min_s'])):.2f})  Mz 峰 {max(abs(d['Mz']['max_s']),abs(d['Mz']['min_s'])):.2f}")
print("平衡:", {k: v for k, v in data["balance"].items() if k != "table"}); [print("   ", r) for r in data["balance"]["table"]]
