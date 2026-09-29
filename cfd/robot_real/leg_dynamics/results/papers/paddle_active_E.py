#!/usr/bin/env python3
"""现有铲形桨绕 E 主动俯仰（ψ 舵机经弹簧驱动）+ E 竖直沉浮：UVLM 扫 θ0、相位
 取向：桨头朝前（照片状态），铰轴 E 在头尾缘之后 L_stem 处 → 相对来流平面形 = delta（直前缘、尾缘前掠成颈）
 头部运动 = E 沉浮 + 俯仰诱导沉浮 r·θ（r = L_stem + c/2，头在铰前 → 抬头则头向上）
 pivot_z > 0：颈部下折 d，头在 E 之下（俯仰附带纵荡 d·θ）
"""
import os, sys, json, math, itertools
import numpy as np
import sweep_fluke as sf
from uvlm_foil import run_case

HERE = os.path.dirname(os.path.abspath(__file__))
U = 0.25; RHO = 1000.0
HEAD_W, HEAD_L = 0.044, 0.042
sf.S_FIX = HEAD_W * HEAD_L * 1.4 / 2
GEO, B_SPAN, C_MID = sf.planform("delta", HEAD_W**2 / sf.S_FIX, 0.4)   # 直前缘（宽边迎流），窄颈在后
C_TIP = GEO["stations"][0][1]; S_HEAD = sf.S_FIX
CD_OFF = 0.05

def one(job):
    f, h0, L, th0, psi, dz = job["f"], job["h0"], job["L"], job["th0"], job["psi"], job.get("dz", 0.0)
    kin = dict(f=f, h0=h0, theta0=th0, psi=psi, pivot=(C_MID + L) / C_TIP, pivot_z=dz)
    try:
        r = run_case(GEO, kin, U=U, n_cycles=4)
    except Exception as e:
        return dict(job, error=repr(e))
    t = np.asarray(r["t_series"]); Fz = np.asarray(r["Fz_series"]); My = np.asarray(r["My_series"])
    hd = np.asarray(r["hdot_series"]); thd = np.asarray(r["thetadot_series"])   # rad/s
    P_h = float(np.mean(-Fz * hd)); P_p = float(np.mean(-My * thd))
    T = float(r["T_mean"]) - CD_OFF * 0.5 * RHO * U**2 * S_HEAD
    w = 2 * math.pi * f; rr = L + 0.5 * HEAD_L
    tt = np.linspace(0, 1 / f, 720)
    th = math.radians(th0) * np.sin(w * tt + math.radians(psi))
    zh = h0 * np.sin(w * tt) + rr * th                       # 头在铰前：抬头 → 头上移
    zd = np.gradient(zh, tt)
    alpha = np.degrees(np.arctan2(zd, U)) - np.degrees(th)
    Hh = h0 + rr * math.radians(th0) * np.exp(1j * math.radians(psi))   # 头部沉浮相量（相对 E 沉浮）
    My_pk = float(np.abs(My).max())
    return dict(job, T=T, T_raw=float(r["T_mean"]), P_h=P_h, P_p=P_p, P=P_h + P_p,
                eta=T * U / max(P_h + P_p, 1e-9), alpha_max=float(np.abs(alpha).max()),
                head_h0=float(abs(Hh)), St_head=2 * f * float(abs(Hh)) / U, My_pk=My_pk,
                z_top=float(zh.max()), z_bot=float(zh.min()))

if __name__ == "__main__":
    from multiprocessing import Pool
    stage = sys.argv[1] if len(sys.argv) > 1 else "A"
    print(f"head(delta): b={B_SPAN*1e3:.1f} c_mid={C_MID*1e3:.1f} c_tip={C_TIP*1e3:.1f} S={S_HEAD*1e4:.1f} cm²", flush=True)
    if stage == "A":
        jobs = [dict(f=f, h0=h0, L=L, th0=th0, psi=psi) for L in (0.05, 0.082) for f in (1.5, 2.0) for h0 in (0.015, 0.02)
                for th0 in (10, 15, 20, 25, 30) for psi in (60, 90, 120, 150, 180)]
        out = "paddle_active_E_A.json"
    else:
        jobs = [dict(f=f, h0=h0, L=L, th0=th0, psi=psi, dz=dz) for (L, f, h0, th0, psi) in json.load(open(os.path.join(HERE, "paddle_active_E_best.json")))
                for dz in (0.0, 0.03, 0.05)]
        out = "paddle_active_E_B.json"
    res = []
    with Pool(2) as pool:
        for i, o in enumerate(pool.imap_unordered(one, jobs)):
            res.append(o)
            if "error" in o: print(i, "ERR", o["error"], flush=True)
            else: print(f"{i+1}/{len(jobs)} L{o['L']*1e3:.0f} f{o['f']} h0{o['h0']*1e3:.0f} θ0{o['th0']} ψ{o['psi']} dz{o.get('dz',0)*1e3:.0f}: T={o['T']*1e3:.1f} mN P={o['P']*1e3:.1f} mW (pitch {o['P_p']*1e3:.1f}) η={o['eta']:.2f} α={o['alpha_max']:.0f}° head_h0={o['head_h0']*1e3:.0f} St_h={o['St_head']:.2f} My_pk={o['My_pk']*1e3:.1f} mN·m", flush=True)
    json.dump(res, open(os.path.join(HERE, out), "w"), indent=1, default=float)
    print("→", out, len(res))
