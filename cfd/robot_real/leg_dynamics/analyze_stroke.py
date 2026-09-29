#!/usr/bin/env python3
"""
行程方案分析：为什么全行程不如半行程 —— 行程分辨的推力分析 + CFD 标定的行程模型 + 候选方案比较
    python3 analyze_stroke.py <results_dir> [out_dir]
输入：results/force_{open,closed}.dat、gait_open.json（单腿 CFD，两态合成）；leg_dynamics.py（运动学 + 叶素）
输出：stroke_data.json；fig_stroke_resolved.png / fig_stroke_calib.png / fig_stroke_map.png / fig_stroke_best.png / fig_stroke_tradeoff.png
模型：叶素准定常力 + 平板附加质量，乘以从 CFD 标定的"行程校正函数" ε(s)（s = 桨头中心自划水开始的行程），
      回收相损失按 CFD 收拢态相位冲量随蜷缩角、回收速度、来流 U 缩放。用于比较行程幅度 / 中心角 / 蜷缩深度，不替代 CFD 验证。
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
BLUE, RED, GRAY, INK, INK2, SURF, ORANGE, GREEN, PURPLE, GRID = "#2a78d6", "#e34948", "#a8a7a1", "#0b0b0b", "#52514e", "#fcfcfb", "#eb6834", "#3a9d5d", "#7b5cd6", "#e6e5e1"

R = os.path.abspath(sys.argv[1]); OUT = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else R
os.makedirs(OUT, exist_ok=True); HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- CFD 数据（两态合成，第 ncyc 周期） ----------------
def load(tag):
    d = np.loadtxt(os.path.join(R, f"force_{tag}.dat"), comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); return d[i]
o, cl = load("open"), load("closed"); g = json.load(open(os.path.join(R, "gait_open.json")))
T = g["T"]; DUTY = g["DUTY"]; SWEEP = g["SWEEP"]; PSI0 = g["psi0_deg"]; PHI_EXT, PHI_RET = g["PHI_EXT"], g["PHI_RET"]
ncyc = int(np.floor(min(o[-1, 0], cl[-1, 0]) / T + 1e-6)); t0 = (ncyc - 1) * T
N = 2500; tt = np.linspace(t0, t0 + T, N, endpoint=False); tau = (tt - t0) / T; dt = T / N
def col(d, j): return np.interp(tt, d[:, 0], d[:, j])
Fx = np.where(tau < DUTY, col(o, 1), col(cl, 1)); Fz = np.where(tau < DUTY, col(o, 3), col(cl, 3))
Fx_o, Fz_o, Fx_c, Fz_c = col(o, 1), col(o, 3), col(cl, 1), col(cl, 3)
def smooth(y, w=0.02, dt_=dt):
    k = max(1, int(round(w / dt_))); ker = np.ones(k) / k; yp = np.concatenate([y[-k:], y, y[:k]]); return np.convolve(yp, ker, "same")[k:-k]

# ---------------- 运动学 ----------------
for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
    if k in g: os.environ[k] = str(g[k])
os.environ["PSI_START"] = str(PSI0)
spec = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(spec); spec.loader.exec_module(ld)
RHO = ld.RHO; O2 = np.array(ld.O2); S_HEAD = ld.STEM + ld.HEAD / 2; S_TIP = ld.STEM + ld.HEAD
psi, phi, _ = ld.gait(tt); P = [ld.pose(p, q) for p, q in zip(psi, phi)]
E = np.array([p["E"] for p in P]); u = np.array([p["u"] for p in P]); n = np.array([p["n"] for p in P])
tip = E + S_TIP * u; head = E + S_HEAD * u
v_tip = np.linalg.norm(np.gradient(tip, tt, axis=0), axis=1); v_head = np.linalg.norm(np.gradient(head, tt, axis=0), axis=1)
def travel(X, mask):
    s = np.zeros(len(X)); d = np.linalg.norm(np.diff(X, axis=0), axis=1); s[1:] = np.cumsum(d * mask[1:]); return s
pw = tau < DUTY
s_tip = travel(tip, pw); s_head = travel(head, pw)
# 桨面法向（指向划水时的迎水面 = 与桨运动方向同向的那一面）与 +x 的夹角
n_signed = n.copy()
vh = np.gradient(head, tt, axis=0)
flip = np.einsum("ij,ij->i", n_signed, vh) < 0; n_signed[flip] *= -1               # 让 n 指向运动方向
theta_n = np.degrees(np.arctan2(-n_signed[:, 1], -n_signed[:, 0]))                  # 水作用在桨上的合力方向（= −运动方向法向）与 +x 的夹角（正 = 偏上）
Fn = Fx * n_signed[:, 0] + Fz * n_signed[:, 1]                                       # CFD 合力在法向上的分量
Fmag = np.hypot(Fx, Fz)

# ---------------- 行程分辨表 ----------------
def imp(mask, y=Fx): return float(np.trapezoid(y[mask], tt[mask]) / T * 1e3)
quarters = []
for a in np.arange(0, DUTY, DUTY / 4):
    k = (tau >= a) & (tau < a + DUTY / 4); i0, i1 = np.argmax(k), len(k) - np.argmax(k[::-1]) - 1
    quarters.append(dict(tau=[float(a), float(a + DUTY / 4)], psi=[float(np.degrees(psi[i0])), float(np.degrees(psi[i1]))], s_tip=[float(s_tip[i0] * 1e3), float(s_tip[i1] * 1e3)],
                         v_tip=[float(v_tip[i0]), float(v_tip[i1])], Fx=imp(k), Fz=imp(k, Fz), Fx_open=imp(k, Fx_o), Fx_closed=imp(k, Fx_c)))
Fxs = smooth(Fx); ip = int(np.argmax(Fxs[pw])); zc = next((i for i in range(ip, int(DUTY * N)) if Fxs[i] <= 0), None)
cum = np.cumsum(Fx * dt) / T * 1e3                                                   # 累计冲量（折算周期平均）
tot_pw = imp(pw)
def at_travel(s_mm):
    i = int(np.argmax(s_tip * 1e3 >= s_mm)); return float(cum[i])
resolved = dict(peak_tau=float(tau[ip]), peak_Fx=float(Fxs[ip] * 1e3), peak_psi=float(np.degrees(psi[ip])), peak_s_tip=float(s_tip[ip] * 1e3),
                zero_tau=float(tau[zc]) if zc else None, zero_s_tip=float(s_tip[zc] * 1e3) if zc else None, total_travel=float(s_tip[int(DUTY * N) - 1] * 1e3),
                power_total=tot_pw, first_half=imp(pw & (tau < DUTY / 2)), second_half=imp(pw & (tau >= DUTY / 2)),
                cum_50mm=at_travel(50), cum_25mm=at_travel(25), cum_75mm=at_travel(75), Fz_first_half=imp(pw & (tau < DUTY / 2), Fz), Fz_second_half=imp(pw & (tau >= DUTY / 2), Fz),
                equal_speed=[dict(tau=float(t_), psi=float(np.degrees(psi[int(t_ * N)])), v_tip=float(v_tip[int(t_ * N)]), Fx=float(Fxs[int(t_ * N)] * 1e3), Fz=float(smooth(Fz)[int(t_ * N)] * 1e3), theta_n=float(theta_n[int(t_ * N)])) for t_ in (0.10, 0.40)],
                quarters=quarters, theta_n_start=float(theta_n[0]), theta_n_peak=float(theta_n[ip]), theta_n_end=float(theta_n[int(DUTY * N) - 1]))

# ---------------- 叶素模型 + 行程校正标定 ----------------
# 法向力模型：N(t) = κ_d(s)·D(t) + κ_a·AM(t)
#   D  = Σ ½ρ c_n |w_n| w_n · w·ds   （叶素准定常压差力，w_n = 水相对叶素的法向速度，含 U）
#   AM = Σ −ρ π w²/4 · ds · a_n       （平板附加质量，减速时为负）
#   κ_d(s)：随桨头中心自划水开始的行程 s 分箱的校正系数（起动涡增强 → 尾迹亏损），κ_a：附加质量系数（三维效应）
#   两者用 open 算例划水相的 CFD 法向力最小二乘（非负）标定；用 closed 算例独立检验。
BINS = np.array([0, 10, 25, 45, 70, 200]) * 1e-3
def strip_terms(psi_t, phi_t, t, U, open_):
    """返回每个时刻的 D、AM（沿 n_signed 的标量）、n_signed、以及供力矩/功率用的分叶素数据"""
    P_ = [ld.pose(p, q) for p, q in zip(psi_t, phi_t)]
    E_ = np.array([p["E"] for p in P_]); u_ = np.array([p["u"] for p in P_]); n_ = np.array([p["n"] for p in P_])
    hd = E_ + S_HEAD * u_; vh_ = np.gradient(hd, t, axis=0)
    ns = n_.copy(); fl = np.einsum("ij,ij->i", ns, vh_) < 0; ns[fl] *= -1
    s_ = np.zeros(len(t)); s_[1:] = np.cumsum(np.linalg.norm(np.diff(hd, axis=0), axis=1))
    D = np.zeros(len(t)); AM = np.zeros(len(t)); parts = []
    for s0, ds, w, cn in ld.blade_strips(open_):
        X = E_ + (s0 + ds / 2) * u_; V = np.gradient(X, t, axis=0); A_ = np.gradient(V, t, axis=0)
        wrel = np.stack([-U - V[:, 0], -V[:, 1]], 1); wn = np.einsum("ij,ij->i", wrel, ns)
        d = 0.5 * RHO * cn * np.abs(wn) * wn * (w * ds); am = -RHO * np.pi * w ** 2 / 4 * ds * np.einsum("ij,ij->i", A_, ns)
        D += d; AM += am; parts.append((X, V, d, am))
    return dict(D=D, AM=AM, n=ns, s=s_, parts=parts, E=E_, u=u_)
def assemble(terms, kd, ka):
    """由 κ_d(s)、κ_a 合成力、力矩、功率"""
    e = kd(terms["s"]); Fx_ = np.zeros(len(e)); Fz_ = np.zeros(len(e)); M = np.zeros(len(e)); Pw = np.zeros(len(e))
    for X, V, d, am in terms["parts"]:
        dN = e * d + ka * am; dF = dN[:, None] * terms["n"]
        Fx_ += dF[:, 0]; Fz_ += dF[:, 1]; r = X - O2; M += r[:, 0] * dF[:, 1] - r[:, 1] * dF[:, 0]; Pw += -(dF[:, 0] * V[:, 0] + dF[:, 1] * V[:, 1])
    return Fx_, Fz_, M, Pw
tp = tt[pw]; term_o = strip_terms(psi[pw], phi[pw], tp, 0.0, True); term_c = strip_terms(psi[pw], phi[pw], tp, 0.0, False)
Fn_cfd = (Fx * term_o["n"][0, 0] * 0 + 0)  # placeholder
ns_o = term_o["n"]; Fn_o = Fx_o[pw] * ns_o[:, 0] + Fz_o[pw] * ns_o[:, 1]; Fn_c = Fx_c[pw] * ns_o[:, 0] + Fz_c[pw] * ns_o[:, 1]
bidx = np.clip(np.searchsorted(BINS, term_o["s"], side="right") - 1, 0, len(BINS) - 2)
Amat = np.zeros((pw.sum(), len(BINS) - 1 + 1))
for j in range(len(BINS) - 1): Amat[:, j] = np.where(bidx == j, term_o["D"], 0.0)
Amat[:, -1] = term_o["AM"]
from scipy.optimize import nnls
coef, _ = nnls(Amat, Fn_o)
kd_bin, ka = coef[:-1], float(coef[-1])
s_mid = 0.5 * (BINS[:-1] + np.minimum(BINS[1:], term_o["s"].max()))
def kd_fn(s): return np.interp(s, s_mid, kd_bin, left=kd_bin[0], right=kd_bin[-1])
Fxm1, Fzm1, Mm1, Pm1 = assemble(term_o, kd_fn, ka); Fxm0, Fzm0, _, _ = assemble(term_o, lambda s: np.ones_like(s), 1.0)
Fxc1, Fzc1, _, _ = assemble(term_c, kd_fn, ka)
def impT(y, t): return float(np.trapezoid(y, t) / T * 1e3)
calib = dict(bins_mm=(BINS[:-1] * 1e3).tolist(), kappa_d=kd_bin.tolist(), kappa_a=ka, s_mid_mm=(s_mid * 1e3).tolist(),
             open_power_cfd=impT(Fx_o[pw], tp), open_power_qs=impT(Fxm0, tp), open_power_model=impT(Fxm1, tp),
             open_peak_cfd=float(smooth(Fx_o)[pw].max() * 1e3), open_peak_model=float(smooth(Fxm1, 0.02, tp[1] - tp[0]).max() * 1e3),
             closed_power_cfd=impT(Fx_c[pw], tp), closed_power_model=impT(Fxc1, tp), closed_peak_cfd=float(smooth(Fx_c)[pw].max() * 1e3), closed_peak_model=float(smooth(Fxc1, 0.02, tp[1] - tp[0]).max() * 1e3),
             r2_open=float(1 - np.sum((Fn_o - (Amat @ coef)) ** 2) / np.sum((Fn_o - Fn_o.mean()) ** 2)))
Nm0 = np.abs(term_o["D"]) + 0 * term_o["AM"]; N_cfd_abs = np.abs(Fn_o); Nm1 = Amat @ coef
# 分解：压差项 κ_d·D 与附加质量项 κ_a·AM 在 x 向的贡献（open 划水相）
Fx_drag = np.zeros(len(tp)); Fx_am = np.zeros(len(tp)); e_o = kd_fn(term_o["s"])
for X, V, d, am in term_o["parts"]:
    Fx_drag += e_o * d * term_o["n"][:, 0]; Fx_am += ka * am * term_o["n"][:, 0]
h1 = tp - tp[0] < DUTY * T / 2; h2 = ~h1
decomp = dict(drag_first=impT(Fx_drag[h1], tp[h1]) * 2 * 0 + float(np.trapezoid(Fx_drag[h1], tp[h1]) / T * 1e3), drag_second=float(np.trapezoid(Fx_drag[h2], tp[h2]) / T * 1e3),
              am_first=float(np.trapezoid(Fx_am[h1], tp[h1]) / T * 1e3), am_second=float(np.trapezoid(Fx_am[h2], tp[h2]) / T * 1e3),
              drag_total=impT(Fx_drag, tp), am_total=impT(Fx_am, tp), cfd_first=float(np.trapezoid(Fx_o[pw][h1], tp[h1]) / T * 1e3), cfd_second=float(np.trapezoid(Fx_o[pw][h2], tp[h2]) / T * 1e3))
calib["decomp"] = decomp

# ---------------- 回收相损失的 CFD 基准 ----------------
tr = g.get("SWITCH_FRAC", 0.3) * (1 - DUTY)
ph = [(DUTY, DUTY + tr), (DUTY + tr, 1 - tr), (1 - tr, 1.0)]
J_sw_in, J_back, J_sw_out = [float(np.trapezoid(Fx[(tau >= a) & (tau < b)], tt[(tau >= a) & (tau < b)]) * 1e3) for a, b in ph]   # mN·s
J_rec0 = J_sw_in + J_back + J_sw_out; J_SW0 = J_sw_in + J_sw_out
T_SW0 = tr * T; DPHI0 = PHI_EXT - PHI_RET; V_SW0 = 0.31; V_BACK0 = 0.15; A_CLOSED = 13.3e-4
WMAX = np.deg2rad(SWEEP) / 2 * np.pi / (DUTY * T)                                     # 现步态髋角峰值角速度（余弦剖面）
try:
    x = json.load(open(os.path.join(R, "xyz_data.json"))); P_rec0 = float(sum(x["cases"]["composite"]["power_phase_mW"][1:])); P_pow0 = float(x["cases"]["composite"]["power_phase_mW"][0])
except Exception:
    P_rec0, P_pow0 = 2.2, 4.5

# ---------------- 候选行程模型 ----------------
def candidate(A_deg, psi_c_deg, dphi_deg=DPHI0, U=0.0, mode="rate", recovery="tuck", t_p=None, n_t=400, keep=False, rf=1.0):
    """行程幅度 A、中心角 psi_c、蜷缩深度 dphi。
       mode='rate'：峰值髋角速度与现步态相同（t_p ∝ A）；mode='period'：划水相时长固定 = 现值（半程 = 慢一半）。
       recovery：'tuck' 舵机最快速率蜷缩（现）、'slowtuck' 蜷缩摊到整个回收相、'notuck' 收拢桨伸展着回扫。"""
    A = np.deg2rad(A_deg)
    W_ = WMAX * rf                                                                    # 峰值髋角速度 = 现值 × rf
    if t_p is None: t_p = A * np.pi / (2 * W_) if mode == "rate" else DUTY * T
    t = np.linspace(0, t_p, n_t, endpoint=False)
    psi_t = np.deg2rad(psi_c_deg) + A / 2 * np.cos(np.pi * t / t_p); phi_t = np.full_like(psi_t, np.deg2rad(PHI_EXT))
    tm = strip_terms(psi_t, phi_t, t, U, True); Fxm, Fzm, Mm, Pw = assemble(tm, kd_fn, ka)
    J_pow = float(np.trapezoid(Fxm, t) * 1e3); Jz = float(np.trapezoid(Fzm, t) * 1e3); P_pow = float(np.trapezoid(Pw, t) * 1e3)
    # 回收相时长：回扫同速率；切换按舵机速率（Δφ 越小越快）
    t_back = A * np.pi / (2 * W_) if mode == "rate" else (1 - DUTY) * T
    t_sw = T_SW0 * (dphi_deg / DPHI0)
    if recovery == "tuck":
        t_rec = max(t_back, 2 * t_sw + 0.05); J_sw = J_SW0 * (dphi_deg / DPHI0) * ((V_SW0 + U) / V_SW0) ** 2
    elif recovery == "slowtuck":
        t_rec = max(t_back, 2 * t_sw + 0.05); t_sw_eff = t_rec / 2
        v_sw = V_SW0 * (dphi_deg / DPHI0) * (T_SW0 / t_sw_eff)                     # 甩出速度 ∝ Δφ / t_sw
        J_sw = J_SW0 * (dphi_deg / DPHI0) ** 2 * (T_SW0 / t_sw_eff) * ((v_sw + U) / v_sw) ** 2 if v_sw > 0 else 0.0
    else:
        t_rec = t_back; J_sw = 0.0
    v_back = V_BACK0 * (A_deg / SWEEP) * (DUTY * T / t_back)                          # 回扫时桨头速度 ∝ 幅度/时间
    if recovery == "notuck":
        J_bk = -0.5 * RHO * 1.2 * A_CLOSED * 0.5 * (v_back * 1.3 + U) ** 2 * t_back * 1e3   # 收拢桨伸展着回扫（余弦剖面 v² 均值 = ½ 峰值²）
    else:
        J_bk = J_back * (A_deg / SWEEP) * (v_back / V_BACK0) * ((v_back + U) / v_back) ** 2 if v_back > 0 else 0.0
    J_rec = J_sw + J_bk; T_ = t_p + t_rec
    P_rec = P_rec0 * T * (abs(J_rec) / abs(J_rec0))                                   # mW·s，按损失冲量比例
    Fmean = (J_pow + J_rec) / T_; Pmean = (P_pow + P_rec) / T_
    out = dict(A=A_deg, psi_c=psi_c_deg, dphi=dphi_deg, U=U, mode=mode, recovery=recovery, rf=rf, t_p=t_p, t_rec=t_rec, T=T_, f=1 / T_, duty=t_p / T_,
               J_pow=J_pow, J_rec=J_rec, J_sw=J_sw, J_bk=J_bk, Fmean=Fmean, Fz_mean=Jz / T_, Pmean=Pmean, eff=Fmean / Pmean if Pmean > 0 else None,
               peak_Fx=float(smooth(Fxm, 0.02, t[1] - t[0]).max() * 1e3), M_peak=float(np.abs(Mm).max() * 1e3),
               v_tip_max=float(np.linalg.norm(np.gradient(tm["E"] + S_TIP * tm["u"], t, axis=0), axis=1).max()), psi_range=[psi_c_deg + A_deg / 2, psi_c_deg - A_deg / 2])
    if keep: out["t"] = t.tolist(); out["Fx"] = (Fxm * 1e3).tolist()
    return out
base = candidate(SWEEP, PSI0 - SWEEP / 2, keep=True)
base_cfd = dict(Fmean=float(Fx.mean() * 1e3), Pmean=P_pow0 + P_rec0, J_pow=tot_pw * T, J_rec=J_rec0)

# ---------------- 扫描 ----------------
PSI_MIN, PSI_MAX = 70.0, 150.0                                                        # 模型扫描允许的髋角范围（机构可达；船体间隙需另查）
PSI_MIN0, PSI_MAX0 = PSI0 - SWEEP, PSI0                                               # 现步态范围（整机间隙已核对）
As = np.arange(15, 91, 5); psics = np.arange(80, 136, 3)
def scan(U, rec, rf=1.0, mode="rate", dphi=DPHI0):
    M_ = np.full((len(As), len(psics)), np.nan)
    for i, A in enumerate(As):
        for j, pc in enumerate(psics):
            if pc + A / 2 > PSI_MAX + 1e-6 or pc - A / 2 < PSI_MIN - 1e-6: continue
            M_[i, j] = candidate(A, pc, dphi, U, mode, rec, n_t=160, rf=rf)["Fmean"]
    return M_
maps = {}; bests = {}
for U in (0.0, 0.05):
    for rf in (1.0, 1.5):
        for rec in ("tuck", "slowtuck"):
            M_ = scan(U, rec, rf); maps[(U, rf, rec)] = M_
            i, j = np.unravel_index(np.nanargmax(M_), M_.shape); bests[f"U{U}_rf{rf}_{rec}"] = candidate(As[i], psics[j], DPHI0, U, "rate", rec, keep=True, rf=rf)
            # 现步态范围内的最优
            Mc = M_.copy()
            for i2, A in enumerate(As):
                for j2, pc in enumerate(psics):
                    if pc + A / 2 > PSI_MAX0 + 1e-6 or pc - A / 2 < PSI_MIN0 - 1e-6: Mc[i2, j2] = np.nan
            i, j = np.unravel_index(np.nanargmax(Mc), Mc.shape); bests[f"U{U}_rf{rf}_{rec}_inrange"] = candidate(As[i], psics[j], DPHI0, U, "rate", rec, keep=True, rf=rf)
def at_all_U(c_):
    return {str(U): {kk: vv for kk, vv in candidate(c_["A"], c_["psi_c"], c_["dphi"], U, c_["mode"], c_["recovery"], n_t=200, rf=c_["rf"]).items() if kk in ("Fmean", "Pmean", "eff", "J_pow", "J_rec", "peak_Fx", "M_peak")} for U in (0.0, 0.05, 0.1)}
schemes = {
    "current":           candidate(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "tuck", keep=True),
    "first_half_rate":   candidate(SWEEP / 2, PSI0 - SWEEP / 4, DPHI0, 0.0, "rate", "tuck", keep=True),
    "first_half_period": candidate(SWEEP / 2, PSI0 - SWEEP / 4, DPHI0, 0.0, "period", "tuck", keep=True),
    "mid_half_rate":     candidate(SWEEP / 2, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "tuck", keep=True),
    "half_rate2":        candidate(SWEEP / 2, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "tuck", keep=True, rf=2.0),
    "full_slowtuck":     candidate(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", keep=True),
    "full_rate1.5_slow": candidate(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", keep=True, rf=1.5),
    "best_U0_rf1_slow_inrange":  bests["U0.0_rf1.0_slowtuck_inrange"],
    "best_U0_rf1_slow":          bests["U0.0_rf1.0_slowtuck"],
    "best_U0_rf1.5_slow_inrange": bests["U0.0_rf1.5_slowtuck_inrange"],
    "best_U005_rf1_slow_inrange": bests["U0.05_rf1.0_slowtuck_inrange"],
    "best_U005_rf1.5_slow":       bests["U0.05_rf1.5_slowtuck"],
}
for k_ in schemes: schemes[k_]["at_U"] = at_all_U(schemes[k_])
sweep_A = {}
for U in (0.0, 0.05):
    for rec in ("tuck", "slowtuck"):
        pc = bests[f"U{U}_rf1.0_{rec}_inrange"]["psi_c"]
        sweep_A[f"U{U}_{rec}"] = dict(psi_c=float(pc), rows=[dict(A=float(A), **{kk: candidate(A, pc, DPHI0, U, "rate", rec, n_t=160)[kk] for kk in ("Fmean", "Pmean", "f", "J_pow", "J_rec", "M_peak", "T")}) for A in np.arange(10, 91, 5) if pc + A / 2 <= PSI_MAX + 1e-6 and pc - A / 2 >= PSI_MIN - 1e-6])
rate_sens = {f"A{A}": [dict(rf=float(rf), **{kk: candidate(A, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", n_t=160, rf=rf)[kk] for kk in ("Fmean", "Pmean", "M_peak", "f", "v_tip_max")}) for rf in (0.75, 1.0, 1.25, 1.5, 2.0)] for A in (30, 60)}
dphi_sens = {f"A{A}_{rec}": [dict(dphi=float(d), **{kk: candidate(A, PSI0 - SWEEP / 2, d, 0.0, "rate", rec, n_t=160)[kk] for kk in ("Fmean", "J_rec", "T")}) for d in (65, 50, 40, 30, 20)] for A in (30, 60) for rec in ("tuck", "slowtuck")}
# 收拢时机扫描（纯数据：open τ<tc，closed τ≥tc）
tc_scan = [dict(tc=float(tc_), Fmean=float(np.trapezoid(Fx_o[tau < tc_], tt[tau < tc_]) / T * 1e3 + np.trapezoid(Fx_c[tau >= tc_], tt[tau >= tc_]) / T * 1e3)) for tc_ in np.arange(0.2, 0.501, 0.025)]
tc_best = max(tc_scan, key=lambda q: q["Fmean"])

data = dict(T=T, DUTY=DUTY, SWEEP=SWEEP, PSI0=PSI0, PHI_EXT=PHI_EXT, PHI_RET=PHI_RET, ncycles=ncyc, wmax_deg_s=float(np.degrees(WMAX)),
            resolved=resolved, calib=calib, recovery=dict(J_sw_in=J_sw_in, J_back=J_back, J_sw_out=J_sw_out, J_rec=J_rec0, P_rec_mW=P_rec0, P_pow_mW=P_pow0, t_sw=T_SW0, dphi0=DPHI0),
            base_model=base, base_cfd=base_cfd, schemes=schemes, sweep_A=sweep_A, dphi_sens=dphi_sens, rate_sens=rate_sens, tc_scan=tc_scan, tc_best=tc_best, psi_limits=[PSI_MIN, PSI_MAX], psi_limits_current=[PSI_MIN0, PSI_MAX0],
            maps={f"U{U}_rf{rf}_{r}": dict(A=As.tolist(), psi_c=psics.tolist(), F=np.where(np.isnan(M_), None, M_).tolist()) for (U, rf, r), M_ in maps.items()})
def _np(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, np.ndarray): return o.tolist()
    raise TypeError(str(type(o)))
json.dump(data, open(os.path.join(OUT, "stroke_data.json"), "w"), indent=1, ensure_ascii=False, default=_np)

# ================= 图 =================
def style(ax):
    ax.set_facecolor(SURF)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(color=GRID, lw=0.8)
# 1) 行程分辨
fig, axs = plt.subplots(2, 2, figsize=(14, 9), facecolor=SURF)
ax = axs[0][0]; st = s_tip[pw] * 1e3
ax.plot(st, Fxs[pw] * 1e3, color=INK, lw=2.2, label="F_x（推力）"); ax.plot(st, smooth(Fz)[pw] * 1e3, color=BLUE, lw=1.6, label="F_z（竖向）")
ax.plot(st, smooth(Fmag)[pw] * 1e3, color=GRAY, lw=1.2, ls="--", label="|F|")
ax.axvline(resolved["zero_s_tip"] or 0, color=RED, lw=0.9, ls=":"); ax.text((resolved["zero_s_tip"] or 0) - 1, 0.12, f"F_x 过零 {resolved['zero_s_tip']:.0f} mm ", color=RED, fontsize=9, transform=ax.get_xaxis_transform(), ha="right")
ax2 = ax.twinx(); ax2.plot(st, v_tip[pw], color=ORANGE, lw=1.4, alpha=0.9); ax2.set_ylabel("桨尖速度 [m/s]", color=ORANGE); ax2.tick_params(axis="y", colors=ORANGE); ax2.set_ylim(0, 0.6); ax2.spines["top"].set_visible(False)
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel("桨尖行程 s [mm]（自划水开始）", color=INK); ax.set_ylabel("力 [mN]（20 ms 平滑）", color=INK); ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.set_title("划水相：力随行程的分布（两态合成 = open 算例；橙线为桨尖速度）", loc="left", fontsize=10.5)
ax = axs[0][1]
ax.plot(st, cum[pw], color=INK, lw=2.2); ax.axhline(tot_pw, color=GRAY, lw=0.9, ls="--")
for s_mm in (25, 50, 75):
    i = int(np.argmax(s_tip * 1e3 >= s_mm)); ax.plot([s_mm], [cum[i]], "o", color=RED, ms=6); ax.text(s_mm + 1.5, cum[i] - 1.8, f"{s_mm} mm: {cum[i]:.1f} ({100*cum[i]/tot_pw:.0f} %)", fontsize=9, color=INK)
style(ax); ax.set_xlabel("桨尖行程 s [mm]", color=INK); ax.set_ylabel("累计推力冲量 / T [mN]", color=INK); ax.set_title(f"累计推力：前 50 mm 产生 {100*resolved['cum_50mm']/tot_pw:.0f} %，全程 {tot_pw:.1f} mN", loc="left", fontsize=10.5)
ax = axs[1][0]
ax.plot(st, theta_n[pw], color=PURPLE, lw=2); ax.axhline(0, color=INK, lw=0.8)
ax.axhline(0, color=INK, lw=0.8); ax.set_ylim(-30, 60)
style(ax); ax.set_xlabel("桨尖行程 s [mm]", color=INK); ax.set_ylabel("桨受力方向与 +x 的夹角 [°]（正 = 偏上）", color=INK)
ax.set_title(f"力的方向（桨面法向）：起点 {resolved['theta_n_start']:.0f}° 偏上，峰值处 {resolved['theta_n_peak']:.0f}°，终点 {resolved['theta_n_end']:.0f}°", loc="left", fontsize=10.5)
ax = axs[1][1]
qs_ = quarters; xq = np.arange(4)
ax.bar(xq - 0.2, [q["Fx"] for q in qs_], 0.4, color=INK, label="F_x 贡献"); ax.bar(xq + 0.2, [q["Fz"] for q in qs_], 0.4, color=BLUE, label="F_z 贡献")
for i, q in enumerate(qs_): ax.text(i - 0.2, q["Fx"] + 0.3, f"{q['Fx']:+.1f}", ha="center", fontsize=9, color=INK); ax.text(i + 0.2, q["Fz"] + 0.3, f"{q['Fz']:+.1f}", ha="center", fontsize=9, color=BLUE)
ax.set_xticks(xq); ax.set_xticklabels([f"τ {q['tau'][0]:.3f}–{q['tau'][1]:.3f}\nψ {q['psi'][0]:.0f}→{q['psi'][1]:.0f}°\n{q['s_tip'][0]:.0f}–{q['s_tip'][1]:.0f} mm" for q in qs_], fontsize=8.5)
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_ylabel("对周期平均的贡献 [mN]", color=INK); ax.legend(frameon=False, fontsize=9); ax.set_title("按四分之一行程分段的贡献", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_stroke_resolved.png"), dpi=110, facecolor=SURF); plt.close()

# 2) 标定
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), facecolor=SURF)
ax = axs[0]; taup = tau[pw]; dtp = tp[1] - tp[0]
ax.plot(taup, -Fn_o * 1e3, color=INK, lw=0.6, alpha=0.35); ax.plot(taup, -smooth(Fn_o, 0.02, dtp) * 1e3, color=INK, lw=2, label="CFD 法向力（逆桨运动方向为正）")
ax.plot(taup, -(term_o["D"] + term_o["AM"]) * 1e3, color=GRAY, lw=1.4, ls="--", label="叶素准定常 + 附加质量（未标定）")
ax.plot(taup, -Nm1 * 1e3, color=GREEN, lw=1.8, label=f"标定模型 κ_d(s)·D + κ_a·AM（R² = {calib['r2_open']:.2f}）")
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel("t / T", color=INK); ax.set_ylabel("法向力 [mN]", color=INK); ax.legend(frameon=False, fontsize=8.5); ax.set_title("open 算例划水相：模型标定", loc="left", fontsize=10.5)
ax = axs[1]
ax.step(np.append(BINS[:-1], term_o["s"].max()) * 1e3, np.append(kd_bin, kd_bin[-1]), where="post", color=GREEN, lw=2.2)
ax.axhline(1, color=GRAY, lw=0.9, ls="--"); ax.text(0.98, 0.9, f"κ_a = {ka:.2f}", transform=ax.transAxes, ha="right", fontsize=10, color=INK)
style(ax); ax.set_xlabel("桨头中心行程 s [mm]", color=INK); ax.set_ylabel("κ_d(s)：压差项校正系数", color=INK); ax.set_title("行程校正系数 κ_d(s)：沿行程增大（力滞后于速度）", loc="left", fontsize=10.5)
ax = axs[2]
ax.plot(taup, smooth(Fx_c, 0.02)[pw] * 1e3, color=BLUE, lw=2, label=f"closed CFD（{calib['closed_power_cfd']:+.1f} mN）"); ax.plot(taup, smooth(Fxc1, 0.02, dtp) * 1e3, color=BLUE, lw=1.6, ls="--", label=f"同一系数预测（{calib['closed_power_model']:+.1f} mN）")
ax.plot(taup, smooth(Fx_o, 0.02)[pw] * 1e3, color=ORANGE, lw=2, label=f"open CFD（{calib['open_power_cfd']:+.1f}）"); ax.plot(taup, smooth(Fxm1, 0.02, dtp) * 1e3, color=ORANGE, lw=1.6, ls="--", label=f"open 模型（{calib['open_power_model']:+.1f}）")
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel("t / T", color=INK); ax.set_ylabel("F_x [mN]", color=INK); ax.legend(frameon=False, fontsize=8.5); ax.set_title("独立检验：用 open 标定的系数预测 closed 桨", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_stroke_calib.png"), dpi=110, facecolor=SURF); plt.close()

# 3) 方案地图：U = 0 / 0.05 × 速率 1.0 / 1.5 × 慢蜷缩（现蜷缩方式的图在 json 里）
fig, axs = plt.subplots(2, 2, figsize=(13, 9.5), facecolor=SURF)
keys = [(0.0, 1.0), (0.0, 1.5), (0.05, 1.0), (0.05, 1.5)]
vmin = np.nanmin([np.nanmin(maps[(U, rf, "slowtuck")]) for U, rf in keys]); vmax = np.nanmax([np.nanmax(maps[(U, rf, "slowtuck")]) for U, rf in keys])
for ax, (U, rf) in zip(axs.ravel(), keys):
    M_ = maps[(U, rf, "slowtuck")]
    im = ax.imshow(M_, origin="lower", aspect="auto", extent=[psics[0] - 1.5, psics[-1] + 1.5, As[0] - 2.5, As[-1] + 2.5], cmap="viridis", vmin=vmin, vmax=vmax)
    # 现步态范围框：ψ ∈ [76,136] → A ≤ 2·min(ψc−76, 136−ψc)
    pcs = np.linspace(psics[0], psics[-1], 200); ax.plot(pcs, 2 * np.minimum(pcs - PSI_MIN0, PSI_MAX0 - pcs), color="white", lw=1.2, ls="--")
    b = bests[f"U{U}_rf{rf}_slowtuck"]; bi = bests[f"U{U}_rf{rf}_slowtuck_inrange"]
    ax.plot(b["psi_c"], b["A"], "*", color=RED, ms=15, mec="white"); ax.plot(bi["psi_c"], bi["A"], "o", color=RED, ms=9, mec="white")
    ax.plot(PSI0 - SWEEP / 2, SWEEP, "s", color=ORANGE, ms=9, mec="white")
    ax.text(0.02, 0.97, f"★ 全范围最优 A={b['A']:.0f}° ψc={b['psi_c']:.0f}°: {b['Fmean']:+.1f} mN, {b['f']:.2f} Hz\n● 现范围内最优 A={bi['A']:.0f}° ψc={bi['psi_c']:.0f}°: {bi['Fmean']:+.1f} mN\n■ 现步态: {candidate(SWEEP, PSI0 - SWEEP / 2, DPHI0, U, 'rate', 'slowtuck', n_t=160, rf=rf)['Fmean']:+.1f} mN（慢蜷缩）", transform=ax.transAxes, va="top", fontsize=8.3, color=INK)
    ax.set_xlabel("行程中心髋角 ψc [°]（90° = 桨竖直）", color=INK); ax.set_ylabel("行程幅度 A [°]", color=INK)
    ax.set_title(f"U = {U} m/s · 峰值髋角速度 = 现值 × {rf}（{np.degrees(WMAX)*rf:.0f}°/s）", loc="left", fontsize=10)
cb = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.8); cb.set_label("周期平均净推力 [mN]（模型，慢蜷缩）", color=INK)
fig.suptitle(f"行程方案地图（频率 ∝ 速率/A；白虚线以下 = 现步态髋角范围 {PSI_MIN0:.0f}–{PSI_MAX0:.0f}°，整机间隙已核对；空白 = 超出 {PSI_MIN:.0f}–{PSI_MAX:.0f}°）", fontsize=10.5, x=0.42)
plt.savefig(os.path.join(OUT, "fig_stroke_map.png"), dpi=110, facecolor=SURF, bbox_inches="tight"); plt.close()

# 4) 幅度曲线 + 方案时间历程 + 收拢时机
fig, axs = plt.subplots(1, 3, figsize=(18, 5), facecolor=SURF)
ax = axs[0]
recname = {"tuck": "蜷缩（舵机最快，现）", "slowtuck": "慢蜷缩（摊到回收相）"}
for rec, col_ in (("tuck", INK), ("slowtuck", GREEN)):
    for U, ls_ in ((0.0, "-"), (0.05, "--")):
        sa = sweep_A[f"U{U}_{rec}"]; ax.plot([q["A"] for q in sa["rows"]], [q["Fmean"] for q in sa["rows"]], ls_, marker="o", color=col_, lw=1.8, ms=4, label=f"{recname[rec]}  U = {U}（ψc {sa['psi_c']:.0f}°）")
ax.axvline(SWEEP, color=GRAY, lw=0.9, ls=":"); ax.text(SWEEP - 0.5, 0.97, "现幅度 60° ", transform=ax.get_xaxis_transform(), fontsize=9, color=INK2, ha="right")
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel("行程幅度 A [°]（同一峰值髋角速度 → 频率 ∝ 1/A）", color=INK); ax.set_ylabel("周期平均净推力 [mN]", color=INK); ax.legend(frameon=False, fontsize=8)
ax.set_title("行程幅度扫描：实线 U = 0，虚线 U = 0.05 m/s", loc="left", fontsize=10.5)
ax = axs[1]
for key, col_, lab in (("current", ORANGE, "现步态 60°，同速率"), ("first_half_rate", GRAY, "前半程 30°（136→106°），同速率"), ("first_half_period", PURPLE, "前半程 30°，同时长（慢一半）"), ("mid_half_rate", BLUE, "中段半程 30°（121→91°），同速率"), ("best_U0_rf1_slow_inrange", GREEN, "现范围内最优 55°（132→76°），慢蜷缩")):
    c_ = schemes[key]; t_ = np.array(c_["t"]); ax.plot(t_, np.array(c_["Fx"]), color=col_, lw=2.4 if key.startswith("best") else 1.6, label=f"{lab}：{c_['at_U']['0.0']['Fmean']:+.1f} mN")
ax.plot(tau[pw] * T, Fxs[pw] * 1e3, color=INK, lw=1.0, ls="--", label="现步态 CFD")
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel("划水相时间 [s]", color=INK); ax.set_ylabel("F_x [mN]（模型，U = 0）", color=INK); ax.legend(frameon=False, fontsize=8.5); ax.set_title("各方案划水相推力时间历程（周期平均值见图例）", loc="left", fontsize=10.5)
ax = axs[2]
ax.plot([q["tc"] for q in tc_scan], [q["Fmean"] for q in tc_scan], "o-", color=INK, lw=2, ms=5)
ax.plot([tc_best["tc"]], [tc_best["Fmean"]], "*", color=RED, ms=15); ax.text(tc_best["tc"], tc_best["Fmean"] + 0.4, f"最优 τ_c = {tc_best['tc']:.3f}：{tc_best['Fmean']:+.1f} mN", ha="center", fontsize=9, color=INK)
ax.axvline(DUTY, color=GRAY, lw=0.9, ls=":"); ax.text(DUTY - 0.005, 0.05, "现：划水结束收拢 ", transform=ax.get_xaxis_transform(), ha="right", fontsize=9, color=INK2)
style(ax); ax.set_xlabel("扇收拢时刻 τ_c（此前取 open 算例，此后取 closed）", color=INK); ax.set_ylabel("单腿周期平均推力 [mN]", color=INK); ax.set_title("收拢时机扫描（纯 CFD 数据重拼）", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_stroke_best.png"), dpi=110, facecolor=SURF); plt.close()

# 5) 推力 – 功率 / 力矩权衡（U = 0，慢蜷缩，速率 1.0 与 1.5）
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), facecolor=SURF)
for ax, ykey, ylab in ((axs[0], "Pmean", "单腿平均机械功率 [mW]"), (axs[1], "M_peak", "髋轴力矩峰值 |M_y| [mN·m]（舵机负载）")):
    for rf, col_, mk_ in ((1.0, INK, "o"), (1.5, GREEN, "^")):
        pts = []
        for A in As:
            for pc in psics:
                if pc + A / 2 > PSI_MAX + 1e-6 or pc - A / 2 < PSI_MIN - 1e-6: continue
                c_ = candidate(A, pc, DPHI0, 0.0, "rate", "slowtuck", n_t=100, rf=rf); pts.append((c_[ykey], c_["Fmean"], A))
        pts = np.array(pts); ax.scatter(pts[:, 0], pts[:, 1], marker=mk_, color=col_, s=20, alpha=0.55, label=f"速率 ×{rf}")
    c0 = schemes["current"]; ax.plot(c0[ykey], c0["Fmean"], "s", color=ORANGE, ms=11, mec="white"); ax.text(c0[ykey], c0["Fmean"] - 1.5, " 现步态", fontsize=9, color=INK)
    b = schemes["best_U0_rf1_slow_inrange"]; ax.plot(b[ykey], b["Fmean"], "*", color=RED, ms=16, mec="white")
    style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel(ylab, color=INK); ax.set_ylabel("周期平均净推力 [mN]（U = 0）", color=INK); ax.legend(frameon=False, fontsize=9)
axs[0].set_title("推力–功率权衡（所有候选行程，慢蜷缩）", loc="left", fontsize=10.5); axs[1].set_title("推力–舵机力矩权衡", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_stroke_tradeoff.png"), dpi=110, facecolor=SURF); plt.close()

# 6) 压差项 / 附加质量项分解
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURF)
ax = axs[0]; taup = tau[pw]
ax.plot(taup, smooth(Fx_o, 0.02)[pw] * 1e3, color=INK, lw=2.2, label="CFD F_x（open）")
ax.plot(taup, smooth(Fx_drag, 0.02, tp[1] - tp[0]) * 1e3, color=GREEN, lw=1.8, label="模型：压差项 κ_d·D（∝ 速度²）")
ax.plot(taup, smooth(Fx_am, 0.02, tp[1] - tp[0]) * 1e3, color=PURPLE, lw=1.8, label="模型：附加质量项 κ_a·AM（∝ 加速度）")
ax.plot(taup, smooth(Fx_drag + Fx_am, 0.02, tp[1] - tp[0]) * 1e3, color=GRAY, lw=1.2, ls="--", label="模型合计")
ax.axvline(DUTY / 2, color=GRAY, lw=0.9, ls=":")
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_xlabel("t / T", color=INK); ax.set_ylabel("F_x [mN]", color=INK); ax.legend(frameon=False, fontsize=8.5); ax.set_title("划水相推力的两个来源：加速段附加质量为正，减速段为负", loc="left", fontsize=10.5)
ax = axs[1]; d_ = decomp
x_ = np.arange(2); w_ = 0.26
ax.bar(x_ - w_, [d_["drag_first"], d_["drag_second"]], w_, color=GREEN, label="压差项"); ax.bar(x_, [d_["am_first"], d_["am_second"]], w_, color=PURPLE, label="附加质量项"); ax.bar(x_ + w_, [d_["cfd_first"], d_["cfd_second"]], w_, color=INK, label="CFD 合计")
for i, (a_, b_, c_) in enumerate(zip([d_["drag_first"], d_["drag_second"]], [d_["am_first"], d_["am_second"]], [d_["cfd_first"], d_["cfd_second"]])):
    for xx, v in ((i - w_, a_), (i, b_), (i + w_, c_)): ax.text(xx, v + (0.4 if v >= 0 else -0.4), f"{v:+.1f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=9, color=INK)
ax.set_xticks(x_); ax.set_xticklabels(["前半行程 τ 0–0.25（加速）", "后半行程 τ 0.25–0.5（减速）"])
style(ax); ax.axhline(0, color=INK, lw=0.8); ax.set_ylabel("对周期平均推力的贡献 [mN]", color=INK); ax.legend(frameon=False, fontsize=9); ax.set_title("前后半行程分解：'前半 65 %' 里有附加质量的预支", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_stroke_decomp.png"), dpi=110, facecolor=SURF); plt.close()

# ---------------- 摘要 ----------------
r = resolved
print(f"划水相 {tot_pw:+.1f} mN：前半 {r['first_half']:+.1f} / 后半 {r['second_half']:+.1f}；前 50 mm {r['cum_50mm']:+.1f}（{100*r['cum_50mm']/tot_pw:.0f}%）；峰 τ={r['peak_tau']:.2f} s={r['peak_s_tip']:.0f} mm；过零 τ={r['zero_tau']:.2f} s={r['zero_s_tip']:.0f}/{r['total_travel']:.0f} mm")
print("等速对比:", [(q['tau'], round(q['v_tip'], 2), round(q['Fx']), round(q['theta_n'])) for q in r["equal_speed"]])
print("标定: κ_d(s) =", [f"{b:.0f}mm:{k:.2f}" for b, k in zip(calib['bins_mm'], calib['kappa_d'])], f"κ_a = {ka:.2f}  R² = {calib['r2_open']:.3f}")
print("  open 划水相 CFD %.1f / 未标定 %.1f / 模型 %.1f（峰 %.0f / %.0f）；closed 检验 CFD %.1f / 模型 %.1f（峰 %.0f / %.0f）" % (calib["open_power_cfd"], calib["open_power_qs"], calib["open_power_model"], calib["open_peak_cfd"], calib["open_peak_model"], calib["closed_power_cfd"], calib["closed_power_model"], calib["closed_peak_cfd"], calib["closed_peak_model"]))
print(f"回收损失冲量 [mN·s]: 切入 {J_sw_in:+.2f} 回扫 {J_back:+.2f} 切出 {J_sw_out:+.2f} 合计 {J_rec0:+.2f}；峰值髋角速度 {np.degrees(WMAX):.0f}°/s；回收功率 {P_rec0:.1f} mW")
print(f"基准（现步态）模型 {base['Fmean']:+.1f} mN vs CFD {base_cfd['Fmean']:+.1f}；功率 模型 {base['Pmean']:.1f} vs CFD {base_cfd['Pmean']:.1f} mW")
for k_, c_ in schemes.items():
    a_ = c_["at_U"]
    print(f"[{k_:<20}] A {c_['A']:4.0f}° ψc {c_['psi_c']:5.1f}° (ψ {c_['psi_range'][0]:.0f}→{c_['psi_range'][1]:.0f}) {c_['mode']:<6} {c_['recovery']:<8} T {c_['T']:.2f} s ({c_['f']:.2f} Hz) duty {c_['duty']:.2f} | F̄ U0 {a_['0.0']['Fmean']:+.1f} U0.05 {a_['0.05']['Fmean']:+.1f} U0.1 {a_['0.1']['Fmean']:+.1f} mN | P(U0.1) {a_['0.1']['Pmean']:.1f} mW | J_pow {c_['J_pow']:+.1f} J_rec {c_['J_rec']:+.1f} mN·s | M峰 {c_['M_peak']:.1f} | v_tip {c_['v_tip_max']:.2f}")
for k_, rows in dphi_sens.items(): print(k_, [(q['dphi'], round(q['Fmean'], 1)) for q in rows])
for k_, rows in rate_sens.items(): print(k_, [(q['rf'], round(q['Fmean'], 1), round(q['Pmean'], 1), round(q['M_peak'], 1)) for q in rows])
print("tc scan best:", tc_best, " current:", tc_scan[-1])
print("分解:", {k: round(v, 1) for k, v in decomp.items()})
