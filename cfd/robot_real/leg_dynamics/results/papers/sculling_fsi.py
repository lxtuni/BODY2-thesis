#!/usr/bin/env python3
"""“摇橹式”：整根杆绕 E 附近摆（ψ 舵机），尾鳍在杆端颈部铰上被动俯仰，铰中立位随杆转
 铰点沉浮 h = r·sinθ_s（θ_s = 杆摆角，幅值 θs0），二阶纵荡 r(1−cosθ_s)
 铰方程：I δ̈ + c δ̇ + k δ = M_hydro(θ_abs)，θ_abs = θ_s + δ
 相量：Θ_abs = [M_h + (k − Iω² + icω)Θ_s] / (k − Iω² + icω − A)，Θ_s 与沉浮同相（相位 0°）
"""
import os, sys, json, math
import numpy as np
import passive_pitch as pp
from uvlm_foil import run_case

HERE = os.path.dirname(os.path.abspath(__file__))
U = pp.U; RHO = pp.RHO

def kin(f, r, ths0, x_p, th0=0.0, psi=90.0):
    h0 = r * math.sin(math.radians(ths0)); A_s = r * (1 - math.cos(math.radians(ths0))) / 2
    return dict(f=f, h0=h0, theta0=th0, psi=psi, pivot=x_p * 1e-3 / pp.C_TIP, surge_amp=A_s, surge_phase=90.0)

def run(f, r, ths0, x_p, th0, psi):
    k_ = kin(f, r, ths0, x_p, th0, psi)
    rr = run_case(pp.GEO, k_, U=U, n_cycles=4)
    t = np.asarray(rr["t_series"]); My = np.asarray(rr["My_series"]); Fz = np.asarray(rr["Fz_series"]); Fx = np.asarray(rr["Fx_series"]); hd = np.asarray(rr["hdot_series"])
    w = 2 * math.pi * f; xs = k_["surge_amp"]; xdot = xs * 2 * w * np.cos(2 * w * t + math.pi / 2)
    P_heave = float(np.mean(-Fz * hd) + np.mean(-Fx * xdot))
    # 杆转动带来的俯仰功：舵机经弹簧对铰做功 = mean(−My·θ̇_s)（弹簧力矩 = 水动力矩，稳态）
    thd_s = math.radians(ths0) * w * np.cos(w * t)          # θ_s 与沉浮同相
    P_stem = float(np.mean(-My * thd_s))
    T = float(rr["T_mean"]); CT = T / (0.5 * RHO * U**2 * 18.5e-4) - 0.03; T_v = CT * 0.5 * RHO * U**2 * 18.5e-4
    P = P_heave + P_stem
    return dict(M1=pp.phasor(My, t, f), T=T_v, P=P, P_heave=P_heave, P_stem=P_stem, eta=T_v * U / max(P, 1e-9))

def solve(f, r, ths0, x_p, k, c, n_iter=2, verbose=False):
    w = 2 * math.pi * f; h0 = r * math.sin(math.radians(ths0))
    I = pp.M_FOIL * (pp.C_MID**2 / 12 + (pp.X_CG - x_p * 1e-3)**2)
    S = k - I * w**2 + 1j * c * w
    TH_s = pp.theta_phasor(ths0, 0.0)
    base = run(f, r, ths0, x_p, 0.0, 90.0); Mh = base["M1"]
    g = (25.0, 90.0); r1 = run(f, r, ths0, x_p, *g); A = (r1["M1"] - Mh) / pp.theta_phasor(*g)
    hist = []
    for it in range(n_iter):
        TH = (Mh + S * TH_s) / (S - A)
        th0, psi = pp.from_phasor(TH); th0 = min(th0, 80.0)
        rv = run(f, r, ths0, x_p, th0, psi)
        resid = abs(rv["M1"] - S * (pp.theta_phasor(th0, psi) - TH_s)) / max(abs(rv["M1"]), 1e-9)
        d = pp.from_phasor(pp.theta_phasor(th0, psi) - TH_s)
        hist.append(dict(th0=th0, psi=psi, delta0=d[0], delta_psi=d[1], resid=resid, T=rv["T"], P=rv["P"], P_heave=rv["P_heave"], P_stem=rv["P_stem"], eta=rv["eta"]))
        if verbose: print(f"  it{it}: θ_abs={th0:.1f}°∠{psi:.0f}° δ={d[0]:.1f}° resid={resid:.2f} T={rv['T']*1e3:.1f} mN P={rv['P']*1e3:.1f} mW η={rv['eta']:.2f}", flush=True)
        A = (rv["M1"] - Mh) / pp.theta_phasor(th0, psi) if th0 > 1 else A
    out = hist[-1]
    tt = np.linspace(0, 1 / f, 720); gam = np.degrees(np.arctan2(h0 * w * np.cos(w * tt), U))
    th = out["th0"] * np.sin(w * tt + math.radians(out["psi"]))
    out.update(f=f, r=r, ths0=ths0, h0=h0, x_p=x_p, k=k, c=c, alpha_max=float(np.abs(gam - th).max()), Kf=float(-A.real), Cf=float(-A.imag / w), Mh=abs(Mh))
    return out

def _job(a):
    try:
        return solve(*a)
    except Exception as e:
        return dict(args=a, error=repr(e))

if __name__ == "__main__":
    from multiprocessing import Pool
    jobs = [(f, 0.075, ths0, -3.0, k, 0.1e-3) for f in (1.5, 2.0) for ths0 in (20.0, 25.0, 30.0) for k in (2e-3, 3e-3, 5e-3, 8e-3)]
    res = []
    with Pool(2) as pool:
        for i, o in enumerate(pool.imap_unordered(_job, jobs)):
            res.append(o)
            if "error" in o: print(i, "ERR", o["error"], flush=True)
            else: print(f"{i+1}/{len(jobs)} f{o['f']} r{o['r']*1e3:.0f} θs{o['ths0']:.0f} (h0 {o['h0']*1e3:.0f}) k{o['k']*1e3:.0f}: θ_abs={o['th0']:.1f}°∠{o['psi']:.0f}° δ={o['delta0']:.1f}° α={o['alpha_max']:.0f}° T={o['T']*1e3:.1f} mN P={o['P']*1e3:.1f} (杆 {o['P_stem']*1e3:.1f}) η={o['eta']:.2f} res={o['resid']:.2f}", flush=True)
    json.dump(res, open(os.path.join(HERE, "sculling_fsi.json"), "w"), indent=1, default=float)
    print("→ sculling_fsi.json", len(res))
