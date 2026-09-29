#!/usr/bin/env python3
"""
尾鳍 / 翼型算例（make_fluke_cfd.py）与阻力型对照（make_leg_cfd.py 的 open+closed）的后处理
============================================================================================
  python3 fluke_forces.py collect ~/run/fluke_L0 L0            # 把 postProcessing 与运动表拷到 results_fluke/L0/
  python3 fluke_forces.py summary                              # 汇总 results_fluke/*，写 summary.json，打印推力–航速对比表

尾鳍：最后 NAVG 个周期平均
  T   = F̄x（+x = 前进方向，即推力）         Fz、峰值
  M_h = 流体对尾鳍绕铰链的力矩（y 分量）= M_O − (H − O) × F
  P   = −⟨F·v_H + M_h·ω⟩（驱动输入功率，杆驱动 + 尾鳍俯仰合计）
  C_T = T / (½ρU²S)，S = gait.json 的 AREA（A2 = 11.5 cm²）；η = T·U / P
阻力型：results_verify/<方案>/force_open|closed.dat 两态合成（划水相取 open：τ ∈ [τ_o, τ_c)），
        报告 V3trap 最优时机（τ_c 0.33, τ_o −0.06）与该航速下重新扫出的最优时机两个值
"""
import argparse, glob, json, os, shutil, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results_fluke"); RHO = 1000.0
FLARE = {"L0": dict(T_mN=9.72, note="Flare Live 单翼基准（整腿：尾鳍 + 脚蹼杆，3 mm 网格）")}
TARGET = {"V1": dict(CT=0.32, eta=0.73, note="Schouveiler 2005 NACA0012, h0/c 0.75, St 0.25, α_max 15°（C_T 0.32 ± 0.01, η 0.73 ± 0.04；引用前请对原文核对）")}

def cat_dat(case, name):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", "forces", "*", name)))
    if not files: return None
    d = np.concatenate([np.atleast_2d(np.loadtxt(f, comments="#")) for f in files])
    d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); return d[i]

def collect(case, name):
    case = os.path.expanduser(case); out = os.path.join(RES, name); os.makedirs(out, exist_ok=True)
    for fn in ("force.dat", "moment.dat"):
        d = cat_dat(case, fn)
        if d is None: print(f"!! {case}: 没有 postProcessing/forces/*/{fn}"); continue
        np.savetxt(os.path.join(out, fn), d, fmt="%.9e"); print(f"{name}/{fn} ← {len(d)} 行，t = {d[0,0]:.4f} … {d[-1,0]:.4f} s")
    for fn in ("gait.json", "kinematics_check.csv"):
        if os.path.exists(os.path.join(case, fn)): shutil.copyfile(os.path.join(case, fn), os.path.join(out, fn))
    lg = os.path.join(case, "log.overPimpleDyMFoam")
    if os.path.exists(lg):
        import re; txt = open(lg, errors="ignore").read()
        m = re.findall(r"ClockTime = ([\d.]+) s", txt); n = len(re.findall(r"^Time = ", txt, re.M))
        json.dump(dict(steps=n, clock_s=float(m[-1]) if m else None, end="End" in txt[-2000:]), open(os.path.join(out, "run.json"), "w"))

def fluke_one(R, navg):
    g = json.load(open(os.path.join(R, "gait.json"))); T = g["T"]; U = g["U"]; S = g["AREA"]; OX, OZ = g["origin"]
    f = np.loadtxt(os.path.join(R, "force.dat")); m = np.loadtxt(os.path.join(R, "moment.dat"))
    k = np.genfromtxt(os.path.join(R, "kinematics_check.csv"), delimiter=",", names=True)
    t = f[:, 0]; tend = min(t[-1], g["cycles"] * T); sel = (t >= tend - navg * T) & (t <= tend)
    if sel.sum() < 20: return None
    ts = t[sel]; F = f[sel, 1:4]; Mo = np.array([np.interp(ts, m[:, 0], m[:, j]) for j in (1, 2, 3)]).T
    Hx = np.interp(ts, k["t"], k["Hx"]); Hz = np.interp(ts, k["t"], k["Hz"])
    vx = np.interp(ts, k["t"], np.gradient(k["Hx"], k["t"])); vz = np.interp(ts, k["t"], np.gradient(k["Hz"], k["t"]))
    om = np.interp(ts, k["t"], np.gradient(np.radians(k["ry_deg"]), k["t"]))           # 绕 +y 的角速度
    rx, rz = Hx - OX, Hz - OZ
    Mh = Mo[:, 1] - (rz * F[:, 0] - rx * F[:, 2])                                        # 绕铰链的 y 力矩
    Pf = F[:, 0] * vx + F[:, 2] * vz + Mh * om                                           # 流体对物体做功
    w = np.gradient(ts); w /= w.sum()                                                    # 非均匀步长加权平均
    avg = lambda y: float(np.sum(y * w))
    Tm = avg(F[:, 0]); P = -avg(Pf); q = 0.5 * RHO * U * U * S
    # 实际攻角（该算例的 U 下）：Flare 约定（x 朝船尾）里相对来流 = (U + v_x,世界, −v_z)，α = 弦角 − 来流角
    fin = np.interp(ts, k["t"], k["fin_deg"])
    al = (fin - np.degrees(np.arctan2(-vz, U + vx)) + 90) % 180 - 90
    return dict(U=U, T_mN=Tm * 1e3, Fz_mN=avg(F[:, 2]) * 1e3, Fx_pp_mN=float(np.ptp(F[:, 0])) * 1e3, Fz_pp_mN=float(np.ptp(F[:, 2])) * 1e3,
                Mh_peak_mNm=float(np.abs(Mh).max()) * 1e3, P_mW=P * 1e3, CT=(Tm / q if U > 0 else None), eta=(Tm * U / P if U > 0 and P > 0 else None),
                T_per_P=(Tm / P if P > 0 else None), cycles_avg=navg, alpha_p75=float(np.percentile(np.abs(al), 75)), S_cm2=S * 1e4)

def drag_one(R, navg, tc_fix=0.33, to_fix=-0.06, own_timing=False):
    """阻力型两态合成（与 fold_timing_sweep.py 同一规则：τ ∈ [τ_o, τ_c) 用 open，其余 closed）
    own_timing=True：按该步态自己的开合时机合成（τ_o = 0，τ_c = FAN_CLOSE，未设则 = DUTY），与 compare_verify.py 的 composite 相同。
    真机现步态必须用它；默认的 0.33 / −0.06 是 V3trap 扫出来的最优时机，套在别的步态上会算错。"""
    def load(tag):
        fn = os.path.join(R, f"force_{tag}.dat")
        if not os.path.exists(fn): return None
        d = np.loadtxt(fn, comments="#"); d = d[np.argsort(d[:, 0])]; return d
    o, c = load("open"), load("closed")
    if o is None or c is None: return None
    g = json.load(open(os.path.join(R, "gait_open.json"))); T = g["T"]; U = g["U"]; cyc = g["cycles"]
    if own_timing:
        FC = g.get("FAN_CLOSE", -1); tc_fix = g["DUTY"] if (FC is None or FC < 0) else FC; to_fix = 0.0
    tq = np.linspace((cyc - navg) * T, cyc * T, 2000, endpoint=False); tau = (tq / T) % 1.0
    Fo = np.interp(tq, o[:, 0], o[:, 1]); Fc = np.interp(tq, c[:, 0], c[:, 1])
    def comp(tc_, to_):
        L = (tc_ - to_) % 1.0; inwin = ((tau - to_) % 1.0) < L
        return float(np.where(inwin, Fo, Fc).mean())
    fix = comp(tc_fix, to_fix)
    best = max(((comp(a_, b_), a_, b_) for a_ in np.arange(0.20, 0.50, 0.01) for b_ in np.arange(-0.20, 0.10, 0.01)), key=lambda x: x[0])
    return dict(U=U, T_fix_mN=fix * 1e3, tau_c_fix=round(float(tc_fix), 3), tau_o_fix=round(float(to_fix), 3), T_best_mN=best[0] * 1e3, tau_c_best=round(best[1], 2), tau_o_best=round(best[2], 2),
                T_open_mN=float(Fo.mean()) * 1e3, T_closed_mN=float(Fc.mean()) * 1e3)

def summary(navg):
    out = {}
    print(f"\n===== 尾鳍 / 翼型（最后 {navg} 个周期） =====")
    print(f"{'算例':<10}{'U m/s':>7}{'T mN':>9}{'Fz mN':>9}{'C_T':>7}{'η':>7}{'P mW':>8}{'T/P N/W':>9}{'|M_h| mN·m':>12}{'α75°':>7}   对照")
    for R in sorted(glob.glob(os.path.join(RES, "*"))):
        n = os.path.basename(R)
        if not os.path.exists(os.path.join(R, "force.dat")): continue
        r = fluke_one(R, navg)
        if r is None: print(f"{n:<10} 数据不够"); continue
        out[n] = r; ref = ""
        if n in FLARE: ref = f"Flare {FLARE[n]['T_mN']:.1f} mN → CFD/Flare = {r['T_mN']/FLARE[n]['T_mN']:.2f}"
        if n in TARGET: ref = f"目标 C_T {TARGET[n]['CT']}、η {TARGET[n]['eta']}"
        fmt = lambda x, p=2: "—" if x is None else f"{x:.{p}f}"
        print(f"{n:<10}{r['U']:>7.3f}{r['T_mN']:>9.1f}{r['Fz_mN']:>9.1f}{fmt(r['CT']):>7}{fmt(r['eta']):>7}{r['P_mW']:>8.1f}{fmt(r['T_per_P']):>9}{r['Mh_peak_mNm']:>12.2f}{fmt(r['alpha_p75'],1):>7}   {ref}")
    RV = os.path.join(HERE, "results_verify"); drag = {}
    for name in ["V3trap"] + sorted(os.path.basename(p) for p in glob.glob(os.path.join(RV, "V3trap_U*"))):
        r = drag_one(os.path.join(RV, name), 1)          # 与 05 号笔记一致：两态合成取最后 1 个周期
        if r: drag[name] = r
    if drag:
        print("\n===== 阻力型 V3trap 两态合成 =====")
        print(f"{'方案':<14}{'U m/s':>7}{'T(τ 0.33/−0.06) mN':>20}{'T 最优 mN':>11}{'τ_c':>6}{'τ_o':>6}")
        for n, r in drag.items():
            print(f"{n:<14}{r['U']:>7.3f}{r['T_fix_mN']:>20.1f}{r['T_best_mN']:>11.1f}{r['tau_c_best']:>6.2f}{r['tau_o_best']:>6.2f}")
    cur = {}
    for name in sorted(os.path.basename(p) for p in glob.glob(os.path.join(RV, "cur_U*"))):
        r = drag_one(os.path.join(RV, name), 1, own_timing=True)
        if r: cur[name] = r
    if cur:
        print("\n===== 阻力型 真机现步态（按自身开合时机合成：τ_o 0 / τ_c = DUTY） =====")
        print(f"{'方案':<14}{'U m/s':>7}{'T mN':>9}{'τ_c':>6}{'T 重扫最优 mN':>15}   （重扫最优只作参考：真机开合由线驱动，时机不可调）")
        for n, r in cur.items():
            print(f"{n:<14}{r['U']:>7.3f}{r['T_fix_mN']:>9.1f}{r['tau_c_fix']:>6.2f}{r['T_best_mN']:>15.1f}")
    # 推力–航速对照（升力型用 L0 同一运动的 L2 系列）
    lift = {r["U"]: r["T_mN"] for n, r in out.items() if n == "L0" or n.startswith("L2")}
    if lift or drag or cur:
        print("\n===== 推力–航速（单腿，mN） =====")
        Us = sorted(set(list(lift) + [r["U"] for r in drag.values()] + [r["U"] for r in cur.values()]))
        dmap = {r["U"]: r["T_best_mN"] for r in drag.values()}
        cmap = {r["U"]: r["T_fix_mN"] for r in cur.values()}
        f1 = lambda m, U: f'{m[U]:.1f}' if U in m else '—'
        print(f"{'U m/s':>7}{'升力型 尾鳍':>12}{'阻力型 现步态':>14}{'阻力型 V3trap':>15}")
        for U in Us:
            print(f"{U:>7.3f}{f1(lift, U):>12}{f1(cmap, U):>14}{f1(dmap, U):>15}")
    os.makedirs(RES, exist_ok=True); json.dump(dict(fluke=out, drag=drag, drag_current=cur), open(os.path.join(RES, "summary.json"), "w"), indent=1)
    print("\n→", os.path.join(RES, "summary.json"))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["collect", "summary"]); ap.add_argument("case", nargs="?"); ap.add_argument("name", nargs="?")
    ap.add_argument("--navg", type=int, default=2)
    a = ap.parse_args()
    if a.cmd == "collect": collect(a.case, a.name)
    else: summary(a.navg)
