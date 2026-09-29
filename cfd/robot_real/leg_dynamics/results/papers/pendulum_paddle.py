#!/usr/bin/env python3
"""现有铲形桨（杆 + 头）绕 E 点（平行四边形末端蓝点）被动摆动：流固耦合（基波相量法）
 结构：I θ̈ + c θ̇ + k θ = M_hydro(t)，θ = 桨轴相对水平的角度（绕 E，向上为正）
 沉浮：E 由 ψ+φ 两自由度逆解沿竖直线运动 z_E = h0 sin ωt（无纵荡）
 水动力：UVLM（头 = 铲形平板，杆的升力忽略）；铰轴在头展中前缘之前 L_stem 处（pivot < 0）
 输出：θ0、相位、头部有效沉浮、攻角、推力、功率、效率；扫 L_stem / f / h0 / k / c
"""
import os, sys, json, math, time
import numpy as np
import sweep_fluke as sf
from uvlm_foil import run_case

HERE = os.path.dirname(os.path.abspath(__file__))
U = 0.25; RHO = 1000.0
HEAD_W, HEAD_L = 0.044, 0.042          # 铲头：宽（展向）×长（弦向），CAD 值；照片件更大，按比例缩放
sf.S_FIX = HEAD_W * HEAD_L * 1.4 / 2   # 鱼尾平面形：前窄后宽，面积 ≈ 12.9 cm²
GEO, B_SPAN, C_MID = sf.planform("fishtail", HEAD_W**2 / sf.S_FIX, 0.4)
C_TIP = GEO["stations"][0][1]
S_HEAD = sf.S_FIX
M_PADDLE = 0.012                        # 杆 + 头 PLA 质量 [kg]（照片件估计，CAD 件 ≈ 8 g）

def kin(f, h0, L_stem, th0=0.0, psi=90.0):
    return dict(f=f, h0=h0, theta0=th0, psi=psi, pivot=-L_stem / C_TIP)

def phasor(sig, t, f):
    w = 2 * math.pi * f; n = len(t)
    return 2.0 / n * np.sum(sig * np.exp(-1j * w * t))

def theta_phasor(th0_deg, psi_deg):
    return math.radians(th0_deg) * np.exp(1j * math.radians(psi_deg - 90.0))

def from_phasor(TH):
    return math.degrees(abs(TH)), (math.degrees(np.angle(TH)) + 90.0 + 180) % 360 - 180

def run(f, h0, L_stem, th0, psi, **kw):
    r = run_case(GEO, kin(f, h0, L_stem, th0, psi), U=U, n_cycles=4, **kw)
    t = np.asarray(r["t_series"]); My = np.asarray(r["My_series"])
    Fz = np.asarray(r["Fz_series"]); hd = np.asarray(r["hdot_series"])
    P_heave = float(np.mean(-Fz * hd))
    T = float(r["T_mean"]); CT = T / (0.5 * RHO * U**2 * S_HEAD) - 0.05   # 平板钝缘，阻力偏置取 0.05
    T_v = CT * 0.5 * RHO * U**2 * S_HEAD
    return dict(M1=phasor(My, t, f), T=T_v, T_raw=T, P=P_heave, eta=T_v * U / max(P_heave, 1e-9),
                Fz1=phasor(Fz, t, f), r=r)

def inertia(L_stem):
    r_cg = L_stem + 0.5 * HEAD_L
    return M_PADDLE * (r_cg**2 + HEAD_L**2 / 12)

def solve(f, h0, L_stem, k, c, I=None, n_iter=2, th_guess=(20.0, 90.0), verbose=False):
    w = 2 * math.pi * f
    if I is None: I = inertia(L_stem)
    base = run(f, h0, L_stem, 0.0, 90.0); Mh = base["M1"]
    TH1 = theta_phasor(*th_guess); r1 = run(f, h0, L_stem, *th_guess)
    A = (r1["M1"] - Mh) / TH1
    hist = []
    for it in range(n_iter):
        TH = Mh / (k - I * w**2 + 1j * c * w - A)
        th0, psi = from_phasor(TH); th0 = min(th0, 70.0)
        rv = run(f, h0, L_stem, th0, psi)
        resid = abs(rv["M1"] - (k - I * w**2 + 1j * c * w) * theta_phasor(th0, psi)) / max(abs(rv["M1"]), 1e-9)
        hist.append(dict(th0=th0, psi=psi, resid=resid, T=rv["T"], T_raw=rv["T_raw"], P=rv["P"], eta=rv["eta"]))
        if verbose: print(f"  it{it}: θ0={th0:.1f}° ψ={psi:.0f}° resid={resid:.2f} T={rv['T']*1e3:.1f} mN P={rv['P']*1e3:.1f} mW η={rv['eta']:.2f}", flush=True)
        A = (rv["M1"] - Mh) / theta_phasor(th0, psi) if th0 > 1 else A
    out = hist[-1]
    # 头部运动学：r_ac = L_stem + 0.5 c（低展弦比平板气动中心近弦中）
    r_ac = L_stem + 0.5 * HEAD_L
    TH = theta_phasor(out["th0"], out["psi"])
    H_head = h0 + r_ac * TH                       # 头部沉浮相量（E 沉浮相量 = h0，相位 0）
    tt = np.linspace(0, 1 / f, 720)
    zdot = np.real(1j * w * H_head * np.exp(1j * w * tt))
    th = np.real(TH * np.exp(1j * w * tt))
    alpha = np.degrees(np.arctan2(zdot, U)) - np.degrees(th)
    out.update(f=f, h0=h0, L_stem=L_stem, k=k, c=c, I=I, Kf=float(-A.real), Cf=float(-A.imag / w),
               Mh=abs(Mh), head_h0=float(abs(H_head)), head_phase=float(np.degrees(np.angle(H_head))),
               St_head=2 * f * float(abs(H_head)) / U, alpha_max=float(np.abs(alpha).max()),
               pitch_induced=float(r_ac * abs(TH)))
    return out

def _job(a):
    f, h0, L, k, c = a
    try:
        o = solve(f, h0, L, k, c); return o
    except Exception as e:
        return dict(f=f, h0=h0, L_stem=L, k=k, c=c, error=repr(e))

if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "sweep"
    print(f"head: b={B_SPAN*1e3:.1f} c_mid={C_MID*1e3:.1f} c_tip={C_TIP*1e3:.1f} S={S_HEAD*1e4:.1f} cm²")
    if stage == "ident":
        for L in (0.03, 0.05, 0.082):
            for f, h0 in ((1.5, 0.02), (2.0, 0.02)):
                w = 2 * math.pi * f
                base = run(f, h0, L, 0.0, 90.0); r1 = run(f, h0, L, 20.0, 90.0)
                A = (r1["M1"] - base["M1"]) / theta_phasor(20.0, 90.0)
                I = inertia(L)
                print(f"L{L*1e3:.0f} f{f} h0{h0*1e3:.0f}: M_h={abs(base['M1'])*1e3:.2f} mN·m ∠{np.degrees(np.angle(base['M1'])):.0f}°  "
                      f"K_f={-A.real*1e3:.2f} mN·m/rad  C_f={-A.imag/w*1e3:.3f} mN·m·s/rad  Iω²={I*w*w*1e3:.2f}  T_h0={base['T']*1e3:.1f} mN", flush=True)
    elif stage == "one":
        L, f, h0, k, c = [float(x) for x in sys.argv[2:7]]
        o = solve(f, h0, L, k, c, verbose=True)
        print({kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in o.items()})
    else:
        from multiprocessing import Pool
        Ls = [0.03, 0.05, 0.082]
        ks = [0.0, 3e-3, 6e-3, 12e-3, 25e-3, 50e-3]
        cs = [0.3e-3, 1.0e-3]
        jobs = [(f, h0, L, k, c) for L in Ls for f, h0 in ((1.5, 0.02), (2.0, 0.02)) for c in cs for k in ks]
        res = []
        with Pool(2) as pool:
            for i, o in enumerate(pool.imap_unordered(_job, jobs)):
                res.append(o)
                if "error" in o: print(i, "ERR", o["error"], flush=True)
                else: print(f"{i+1}/{len(jobs)} L{o['L_stem']*1e3:.0f} f{o['f']} c{o['c']*1e3:.1f} k{o['k']*1e3:.0f}: θ0={o['th0']:.1f}° ψ={o['psi']:.0f}° head_h0={o['head_h0']*1e3:.0f}mm St_h={o['St_head']:.2f} α={o['alpha_max']:.0f}° T={o['T']*1e3:.1f} mN P={o['P']*1e3:.1f} mW η={o['eta']:.2f} res={o['resid']:.2f}", flush=True)
        json.dump(res, open(os.path.join(HERE, "pendulum_paddle_sweep.json"), "w"), indent=1, default=float)
        print("→ pendulum_paddle_sweep.json", len(res))
