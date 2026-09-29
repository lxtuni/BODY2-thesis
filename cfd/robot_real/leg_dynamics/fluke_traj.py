#!/usr/bin/env python3
"""
尾鳍（升力型）单腿 CFD 的运动表生成器 → make_fluke_cfd.py --motion <csv>
===========================================================================
输出 CSV（一个周期，等间隔）：t_s, hx_m, hz_m, fin_deg, shank_deg, phase
  hx, hz   铰链（尾鳍销轴）相对“原点”（腿 = 髋轴；翼型基准 = 俯仰轴平均位置）的位置 [m]，
           Flare 同一套约定：hx 朝船尾为正，hz 向上为正
  fin_deg  尾鳍弦线绝对角：0 = 水平朝船尾，正 = 尾缘上翘（Flare kin_actual 的 ch）
  shank_deg 脚蹼杆绝对角（同一约定），只用于记录 / 画图

三种来源：
  --bench  flare-live/bench_LF_active.csv   Flare 单翼基准的“实测”关节角，按相位平均（L0：与 Flare 完全同一运动）
  --alpha  20 --U 0.223 --freq 2            Flare run_float.py 的 fan 步态 + 攻角律（逐行移植，含尾鳍限速 35 rad/s 与 10 ms 低通）
  --heave-pitch H0,THETA0,PSI --freq F      正弦沉浮 + 正弦俯仰（翼型基准 V1）；THETA0 = auto 时按 --alpha-max 反解

例：
  python3 fluke_traj.py --bench ~/flare/bench_LF_active.csv --freq 2 --out motion_L0.csv
  python3 fluke_traj.py --alpha 20 --U 0     --freq 2 --out motion_L2_U000.csv
  python3 fluke_traj.py --heave-pitch 0.075,auto,90 --U 0.4 --St 0.25 --alpha-max 15 --out motion_V1.csv
"""
import argparse, json, sys
import numpy as np

ap = argparse.ArgumentParser()
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--bench", help="Flare bench_LF_active.csv")
g.add_argument("--alpha", type=float, help="目标攻角 [deg]（Flare alpha 模式）")
g.add_argument("--heave-pitch", help="H0[m],THETA0[deg|auto],PSI[deg]")
ap.add_argument("--U", type=float, default=0.0, help="来流 / 航速 [m/s]")
ap.add_argument("--freq", type=float, default=2.0)
ap.add_argument("--St", type=float, default=None, help="heave-pitch：给 St 则 freq = St·U/(2·H0)")
ap.add_argument("--alpha-max", type=float, default=15.0, help="heave-pitch + auto：反解 θ0 使 max|α| = 该值")
ap.add_argument("--amp", type=float, default=1.0, help="脚蹼杆摆幅缩放（Flare AMP）")
ap.add_argument("--u-min", type=float, default=0.05, help="攻角律用的来流下限（Flare U_MIN）；U=0 时尾鳍不会每半周期翻 90°")
ap.add_argument("--n", type=int, default=400, help="每周期点数")
ap.add_argument("--nh", type=int, default=8, help="攻角律尾鳍角的平滑谐波数")
ap.add_argument("--nh-hinge", type=int, default=12, help="攻角律：脚蹼杆角（铰链轨迹）的平滑谐波数；0 = 不平滑（旧版，控制点处铰链速度有折点）")
ap.add_argument("--out", required=True)
a = ap.parse_args()

# ---------------------------------------------------------------- Flare run_float.py 的步态（逐行移植，数值与 Flare 一致）
CP = [(-39.7, -14.5), (-42.0, -24.0), (-40.3, -43.4), (-36.4, -52.4), (-36.4, -52.4), (-40.3, -43.4), (-42.0, -24.0)]
PH = [0.0, 0.167, 0.333, 0.500, 0.530, 0.687, 0.843]
FIN_LIM = (-75.0, 60.0)
L1, L2, TH0, SH0 = 40.0, 46.66, -90.0, -25.0
FIN_RATE, FIN_TAU = 35.0, 0.01          # rad/s, s（Flare 尾鳍指令限速与低通）
RHO, AREA = 1000., 0.034 * 0.060        # 只用来在 ±α 两支里挑推力更大的一支，量级无关

def traj(p):
    p = p % 1.0; pts = CP + [CP[0]]; ph = PH + [1.0]
    for i in range(len(ph) - 1):
        if ph[i] <= p <= ph[i + 1]:
            w = (p - ph[i]) / max(ph[i + 1] - ph[i], 1e-9)
            return np.array(pts[i]) * (1 - w) + np.array(pts[i + 1]) * w
    return np.array(CP[0])
def ik(dx, dz):
    r = min(np.hypot(dx, dz), L1 + L2 - 0.5); b = np.arctan2(dz, dx)
    c2 = np.clip((r * r - L1 * L1 - L2 * L2) / (2 * L1 * L2), -1, 1); s2 = np.sqrt(1 - c2 * c2)
    t1 = b - np.arctan2(L2 * s2, L1 + L2 * c2); return np.degrees(t1), np.degrees(t1 + np.arctan2(s2, c2))
def hinge_cp(p):
    x, y = traj(p); return np.array([-x, y])
def shank_cp(p):
    S = hinge_cp(p); return ik(S[0], S[1])[1]
_TH = np.array([shank_cp(q) for q in np.linspace(0, 1, 400, endpoint=False)])
FAN_CENTER = (_TH.max() + _TH.min()) / 2
def shank_fan(p): return (shank_cp(p) - FAN_CENTER) * a.amp
def hinge(p):
    th = np.radians(shank_fan(p)); K = np.array([0.0, -L1])
    return K + L2 * np.array([np.cos(th), np.sin(th)])        # mm，(朝船尾, 向上)，相对髋
def fwd(ch, Wx, Wz):
    Wm = np.hypot(Wx, Wz)
    if Wm < 1e-9: return 0.0
    phi = np.arctan2(Wz, Wx); al = (ch - phi + np.pi) % (2 * np.pi) - np.pi
    q = 0.5 * RHO * Wm * Wm * AREA
    return -(q * 1.8 * np.sin(2 * al) * (Wz / Wm) + q * (0.05 + 1.0 * (1 - np.cos(2 * al))) * (Wx / Wm))
def fin_alpha(p, U, freq, alpha, V=None, t2=None):
    if V is None:
        d = 1e-3; V = (hinge(p + d) - hinge(p - d)) / (2 * d) * freq / 1000.0
    Wx, Wz = U - V[0], -V[1]
    c = max([np.arctan2(Wz, Wx) + np.radians(alpha), np.arctan2(Wz, Wx) - np.radians(alpha)], key=lambda c_: fwd(c_, Wx, Wz))
    if t2 is None: t2 = shank_fan(p)
    q3 = (np.degrees(c) - t2 + 180) % 360 - 180
    return t2 + float(np.clip(q3, FIN_LIM[0], FIN_LIM[1]))     # 绝对弦角 deg

def periodic_smooth(y, nh=12):
    """周期序列只留前 nh 阶谐波（去掉限速折角带来的高频，保证运动表可导）"""
    Y = np.fft.rfft(y); Y[nh + 1:] = 0; return np.fft.irfft(Y, len(y))

f = a.freq; meta = {}
if a.heave_pitch:
    h0, th0, psi = a.heave_pitch.split(","); h0 = float(h0); psi = np.radians(float(psi))
    if a.St: f = a.St * a.U / (2 * h0)
    ph = np.linspace(0, 1, a.n, endpoint=False); w = 2 * np.pi * f
    hz = h0 * np.sin(2 * np.pi * ph); hzd = h0 * w * np.cos(2 * np.pi * ph)
    def amax(t0):   # 弦角 fin = −θ0·sin(ωt+ψ)（ψ=90° 时上行最快处前缘抬头）；来流 −x、铰链上行 → 相对来流从前上方来
        fin = -t0 * np.sin(2 * np.pi * ph + psi)
        inflow = np.degrees(np.arctan2(hzd, a.U))          # 相对来流偏离水平的角（正 = 来自上方）
        return np.max(np.abs(inflow + fin))                 # 前缘抬头 = fin 取负；两者相消即攻角
    if th0 == "auto":
        from scipy.optimize import brentq
        th0 = brentq(lambda t0: amax(t0) - a.alpha_max, 0.0, np.degrees(np.arctan(w * h0 / a.U)))   # 取 θ0 < 最大来流角的那一支（另一支俯仰过头）
    th0 = float(th0)
    fin = -th0 * np.sin(2 * np.pi * ph + psi)
    hx = np.zeros_like(ph); shank = np.zeros_like(ph)
    meta = dict(mode="heave-pitch", h0=h0, theta0_deg=th0, psi_deg=np.degrees(psi), alpha_max_deg=float(amax(th0)), St=2 * h0 * f / a.U if a.U else None)
elif a.bench:
    import csv
    r = list(csv.DictReader(open(a.bench)))
    d = {k: np.array([float(x[k]) for x in r]) for k in ("t", "LF_shank_act", "LF_fin_rel_act")}
    SETTLE = 1.5; m = d["t"] > 6.0                                    # 周期修正收敛后（bench 6 s 以后）
    p = (f * (d["t"][m] - SETTLE)) % 1.0
    sa = d["LF_shank_act"][m]; fa = sa + d["LF_fin_rel_act"][m]
    nb = a.n; ph = np.linspace(0, 1, nb, endpoint=False)
    def pavg(y):                                                      # 相位平均（圆形插值到均匀相位）
        o = np.argsort(p); pp, yy = p[o], y[o]
        return np.interp(ph, np.r_[pp - 1, pp, pp + 1], np.r_[yy, yy, yy])
    t2 = periodic_smooth(pavg(sa)) + SH0; ch = periodic_smooth(pavg(fa)) + SH0
    shank, fin = t2, ch
    hx = L2 * np.cos(np.radians(t2)) / 1000; hz = (-L1 + L2 * np.sin(np.radians(t2))) / 1000
    meta = dict(mode="bench", source=a.bench, samples=int(m.sum()), cycles=float((d["t"][m][-1] - d["t"][m][0]) * f))
else:
    ph = np.linspace(0, 1, a.n, endpoint=False)
    U_law = max(a.U, a.u_min)
    shank = np.array([shank_fan(q) for q in ph])
    if a.nh_hinge > 0:
        # 分段线性控制点 → 铰链速度在 7 个控制点处跳变（附加质量力尖峰）；零相位截断到 nh_hinge 阶，
        # 与 L0（bench 实测，12 阶）同一处理。尾鳍攻角律用平滑后的铰链速度重新计算。
        shank = periodic_smooth(shank, nh=a.nh_hinge)
    thr = np.radians(shank)
    H = np.c_[L2 * np.cos(thr), -L1 + L2 * np.sin(thr)] / 1000; hx, hz = H[:, 0], H[:, 1]
    k = np.fft.rfftfreq(a.n, 1.0 / a.n)                                  # 谱导数（周期、精确）
    Vh = np.c_[np.fft.irfft(1j * 2 * np.pi * k * np.fft.rfft(hx), a.n), np.fft.irfft(1j * 2 * np.pi * k * np.fft.rfft(hz), a.n)] * f
    raw = np.array([fin_alpha(q, U_law, f, a.alpha, V=Vh[i], t2=shank[i]) for i, q in enumerate(ph)])
    # Flare 里尾鳍是舵机跟踪（限速 35 rad/s + 10 ms 低通 + 相位超前补偿），实测攻角 ≈ 目标值。
    # CFD 里运动是精确给定的，不需要舵机模型；只做零相位平滑（周期傅里叶截断），把 ±α 两支切换处的折角抹圆，
    # 否则切换瞬间角速度无穷大、时间步会被压到 0。NH 越小越圆滑、越偏离目标攻角。
    fin = periodic_smooth(np.unwrap(np.radians(raw)) * 180 / np.pi, nh=a.nh)
    meta = dict(mode="alpha", alpha_deg=a.alpha, U=a.U, U_law=U_law, amp=a.amp, fan_center_deg=float(FAN_CENTER),
                smooth_harmonics=a.nh, smooth_harmonics_hinge=a.nh_hinge)

# ---------------------------------------------------------------- 输出 + 体检
T = 1 / f; t = ph * T
np.savetxt(a.out, np.column_stack([t, hx, hz, fin, shank, ph]), delimiter=",", fmt="%.7f", header="t_s,hx_m,hz_m,fin_deg,shank_deg,phase", comments="")
w = 2 * np.pi * f
dh = np.gradient(np.c_[hx, hz], t, axis=0); dfin = np.gradient(np.unwrap(np.radians(fin)), t)
# 实际攻角（相对来流 = (U 朝船尾吹, 0) − 铰链速度；与 Flare kin_actual 同式）
Wx, Wz = a.U - dh[:, 0], -dh[:, 1]
alpha_act = (np.degrees(np.radians(fin) - np.arctan2(Wz, Wx)) + 90) % 180 - 90
meta.update(freq=f, T=T, U=a.U, n=a.n, hinge_travel_z_mm=float(np.ptp(hz) * 1000), hinge_travel_x_mm=float(np.ptp(hx) * 1000),
            fin_range_deg=[float(fin.min()), float(fin.max())], fin_rate_max_deg_s=float(np.degrees(np.abs(dfin).max())),
            hinge_speed_max=float(np.linalg.norm(dh, axis=1).max()), alpha_act_p75=float(np.percentile(np.abs(alpha_act), 75)),
            alpha_act_max=float(np.abs(alpha_act).max()))
json.dump(meta, open(a.out.rsplit(".", 1)[0] + ".json", "w"), indent=1)
print(f"→ {a.out}  ({meta['mode']}, f {f:.3f} Hz, U {a.U} m/s)  铰链竖向行程 {meta['hinge_travel_z_mm']:.1f} mm  尾鳍角 {fin.min():.1f}…{fin.max():.1f}°  "
      f"尾鳍最大角速度 {meta['fin_rate_max_deg_s']:.0f} °/s  实际攻角 75 分位 {meta['alpha_act_p75']:.1f}°")
