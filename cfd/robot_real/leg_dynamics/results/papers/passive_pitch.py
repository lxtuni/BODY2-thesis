#!/usr/bin/env python3
"""被动俯仰（扭簧铰）水翼 + 平行四边形竖直沉浮：流固耦合（基波相量法）
 结构：I θ̈ + c θ̇ + k θ = M_hydro(t)   （θ 绝对俯仰角，M 绕铰轴，UVLM 给出、含附加质量与尾迹）
 流体：UVLM 对 θ 近似线性 → M1 = M_h + A·Θ（复相量）。识别 M_h、A 后解 Θ = M_h/(k − Iω² + icω − A)，再用完整 UVLM 验证。
 沉浮由平行四边形给出：z = h0 sin ωt，附带二倍频纵荡 x = A_s sin(2ωt+90°)，A_s = 39(1−cos(asin(h0/39)))/2。
 注意 uvlm_foil 的 pivot 参数 = 铰轴在展中前缘之后的距离 / 翼尖弦（14.3 mm）。本文件用 x_p [mm] 直接给。
"""
import os, sys, json, math, time
import numpy as np
from uvlm_foil import run_case
from sweep_fluke import planform

HERE = os.path.dirname(os.path.abspath(__file__))
U = 0.25; RHO = 1000.0
GEO, B_SPAN, C_MID = planform("lunate", 4.0, 0.5, 35.0)
C_TIP = GEO["stations"][0][1]          # 14.3 mm：uvlm_foil 的 pivot 参照弦
R_PAR = 0.039                           # 平行四边形长边 [m]
M_FOIL = 0.0052                         # PLA 翼质量 [kg]
X_CG = 0.45 * C_MID                     # 厚截面翼的质心（距展中前缘）

def kin(f, h0, x_p, th0=0.0, psi=90.0):
    A_s = R_PAR * (1 - math.cos(math.asin(h0 / R_PAR))) / 2
    return dict(f=f, h0=h0, theta0=th0, psi=psi, pivot=x_p * 1e-3 / C_TIP, surge_amp=A_s, surge_phase=90.0)

def phasor(sig, t, f):
    """基波复相量：sig(t) ≈ Re[S e^{iωt}]（最后一个周期）"""
    w = 2 * math.pi * f; n = len(t)
    return 2.0 / n * np.sum(sig * np.exp(-1j * w * t))

def theta_phasor(th0_deg, psi_deg):
    return math.radians(th0_deg) * np.exp(1j * math.radians(psi_deg - 90.0))

def from_phasor(TH):
    return math.degrees(abs(TH)), (math.degrees(np.angle(TH)) + 90.0 + 180) % 360 - 180

def run(f, h0, x_p, th0, psi, **kw):
    r = run_case(GEO, kin(f, h0, x_p, th0, psi), U=U, n_cycles=4, **kw)
    t = np.asarray(r["t_series"]); My = np.asarray(r["My_series"])
    Fz = np.asarray(r["Fz_series"]); Fx = np.asarray(r["Fx_series"]); hd = np.asarray(r["hdot_series"])
    w = 2 * math.pi * f
    xs = kin(f, h0, x_p)["surge_amp"]; xdot = xs * 2 * w * np.cos(2 * w * t + math.pi / 2)
    P_heave = float(np.mean(-Fz * hd) + np.mean(-Fx * xdot))        # 沉浮舵机做的功（含纵荡）
    T = float(r["T_mean"]); CT = T / (0.5 * RHO * U**2 * 18.5e-4) - 0.03
    T_v = CT * 0.5 * RHO * U**2 * 18.5e-4
    return dict(M1=phasor(My, t, f), T=T_v, P=P_heave, eta=T_v * U / max(P_heave, 1e-9), P_pitch=float(np.mean(-My * np.asarray(r["thetadot_series"]))), r=r)

def _job(a):
    f, h0, x_p, k, c = a
    try:
        o = solve(f, h0, x_p, k, c); o.pop("r", None); return o
    except Exception as e:
        return dict(f=f, h0=h0, x_p=x_p, k=k, c=c, error=repr(e))

def solve(f, h0, x_p, k, c, I=None, n_iter=2, th_guess=(25.0, 90.0), verbose=False):
    """给定弹簧 k [N·m/rad]、阻尼 c [N·m·s/rad]，求被动俯仰的稳态基波响应"""
    w = 2 * math.pi * f
    if I is None:
        I = M_FOIL * (C_MID**2 / 12 + (X_CG - x_p * 1e-3)**2)
    base = run(f, h0, x_p, 0.0, 90.0); Mh = base["M1"]
    TH1 = theta_phasor(*th_guess); r1 = run(f, h0, x_p, *th_guess)
    A = (r1["M1"] - Mh) / TH1
    hist = []
    for it in range(n_iter):
        TH = Mh / (k - I * w**2 + 1j * c * w - A)
        th0, psi = from_phasor(TH)
        th0 = min(th0, 80.0)
        rv = run(f, h0, x_p, th0, psi)
        resid = abs(rv["M1"] - (k - I * w**2 + 1j * c * w) * theta_phasor(th0, psi)) / max(abs(rv["M1"]), 1e-9)
        hist.append(dict(th0=th0, psi=psi, resid=resid, T=rv["T"], P=rv["P"], eta=rv["eta"], P_pitch=rv["P_pitch"]))
        if verbose: print(f"  it{it}: θ0={th0:.1f}° ψ={psi:.0f}° resid={resid:.2f} T={rv['T']*1e3:.1f} mN P={rv['P']*1e3:.1f} mW η={rv['eta']:.2f}")
        # 割线更新 A
        A = (rv["M1"] - Mh) / theta_phasor(th0, psi) if th0 > 1 else A
    ind = math.degrees(math.atan(h0 * w / U))
    out = hist[-1]; out.update(f=f, h0=h0, x_p=x_p, k=k, c=c, I=I, chi=out["th0"] / ind, ind=ind,
                               Kf=float(-A.real), Cf=float(-A.imag / w), alpha_max=None)
    # 运动学攻角
    tt = np.linspace(0, 1 / f, 720); gam = np.degrees(np.arctan2(h0 * w * np.cos(w * tt), U))
    th = out["th0"] * np.sin(w * tt + math.radians(out["psi"]))
    out["alpha_max"] = float(np.abs(gam - th).max())
    return out

if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "sweep"
    if stage == "ident":
        for f, h0 in ((1.5, 0.03), (2.0, 0.03)):
            for x_p in (6.0, 9.0, 12.0):
                w = 2 * math.pi * f
                base = run(f, h0, x_p, 0.0, 90.0); r1 = run(f, h0, x_p, 25.0, 90.0)
                A = (r1["M1"] - base["M1"]) / theta_phasor(25.0, 90.0)
                print(f"f{f} h0{h0*1e3:.0f} x_p{x_p:.0f}mm: M_h={abs(base['M1'])*1e3:.2f} mN·m ∠{np.degrees(np.angle(base['M1'])):.0f}°  K_f(流体刚度)={-A.real*1e3:.2f} mN·m/rad  C_f(流体阻尼)={-A.imag/w*1e3:.3f} mN·m·s/rad")
    else:
        from multiprocessing import Pool
        ks = [1e-3, 2e-3, 3e-3, 4e-3, 6e-3, 9e-3, 14e-3]
        jobs = [(f, h0, x_p, k, c) for f, h0 in ((1.5, 0.03), (2.0, 0.03)) for x_p in (6.0, 9.0, 12.0)
                for c in (0.05e-3, 0.15e-3, 0.3e-3) for k in ks]
        res = []
        with Pool(2) as pool:
            for i, o in enumerate(pool.imap_unordered(_job, jobs)):
                res.append(o)
                if "error" in o: print(i, "ERR", o["error"], flush=True)
                else: print(f"{i+1}/{len(jobs)} f{o['f']} x_p{o['x_p']:.0f} c{o['c']*1e3:.2f} k{o['k']*1e3:.0f}: θ0={o['th0']:.1f}° ψ={o['psi']:.0f}° χ={o['chi']:.2f} α_max={o['alpha_max']:.0f}° T={o['T']*1e3:.1f} mN P={o['P']*1e3:.1f} mW η={o['eta']:.2f} resid={o['resid']:.2f}", flush=True)
        json.dump(res, open(os.path.join(HERE, "passive_pitch_sweep.json"), "w"), indent=1, default=float)
        print("→ passive_pitch_sweep.json", len(res))
