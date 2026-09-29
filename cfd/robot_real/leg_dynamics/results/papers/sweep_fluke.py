#!/usr/bin/env python3
"""海豚尾鳍式脚蹼：几何 × 运动学 扫描（UVLM，uvlm_foil.run_case）
目标：U = 0.25 m/s 巡航，S = 18.5 cm² 固定，f ≤ 3 Hz，h0 ≤ 50 mm，α_max ≤ 30°，求最大周期平均推力。
用法：python3 sweep_fluke.py A|B|C   → results/papers/sweep_<stage>.csv
"""
import sys, os, json, time, itertools, math
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
S_FIX = 18.5e-4          # m²
U = 0.25
RHO = 1000.0

# ---------------- 平面形族（面积固定；y 从一个翼尖 0 到另一个翼尖 b；支杆在展中）
def planform(family, AR, taper=1.0, sweep_deg=0.0, n_st=7):
    """返回 stations [(y, c, x_le)]，弦分布关于展中对称。
    family: rect | fishtail(尾缘直、前缘后掠，前窄后宽) | delta(前缘直、尾缘前掠) | ellipse | lunate(前后缘同时后掠)
    taper = 翼尖弦 / 展中弦；sweep_deg 只对 lunate 有意义。"""
    b = math.sqrt(AR * S_FIX)
    ys = np.linspace(0, b, n_st)
    s = np.abs(2 * ys / b - 1)                      # 0 展中 → 1 翼尖
    if family == "ellipse":
        shape = np.sqrt(np.clip(1 - s**2, 0, 1)); shape = np.maximum(shape, 0.08)
    else:
        shape = 1 - (1 - taper) * s
    c_mid = S_FIX / np.trapezoid(shape, ys)
    c = c_mid * shape
    if family in ("rect", "delta"):
        x_le = np.zeros_like(ys)                   # 前缘直
    elif family in ("fishtail",):
        x_le = c_mid - c                            # 尾缘直（x = c_mid），前缘后掠
    elif family == "ellipse":
        x_le = 0.5 * (c_mid - c)                    # 前后缘对称椭圆
    elif family == "lunate":
        x_le = np.abs(ys - b / 2) * math.tan(math.radians(sweep_deg))
    else:
        raise ValueError(family)
    return dict(stations=[(float(y), float(cc), float(x)) for y, cc, x in zip(ys, c, x_le)],
                airfoil="naca0012"), b, float(c_mid)


def one(job):
    from uvlm_foil import run_case, theta0_for_alpha_max
    g, b, c_mid = planform(job["family"], job["AR"], job.get("taper", 1.0), job.get("sweep", 0.0))
    f, h0, amax, psi, piv = job["f"], job["h0"], job["amax"], job["psi"], job["pivot"]
    try:
        th0 = theta0_for_alpha_max(f, h0, U, psi, amax)   # α_max 是上限
    except ValueError as e:
        if "Pure heave" in str(e):
            th0 = 0.0                                     # 诱导角本来就不到上限 → 纯沉浮
        else:
            return dict(job, error="alpha cap unreachable")   # ψ≠90° 时正弦俯仰压不住攻角 → 违反约束，剔除
    kin = dict(f=f, h0=h0, theta0=th0, psi=psi, pivot=piv)
    t0 = time.time()
    try:
        r = run_case(g, kin, U=U, n_cycles=4, n_chord=6, n_span=12, dt_per_cycle=40)
        if not r.get("ok", True):
            return dict(job, error=str(r.get("error")))
        out = dict(job, theta0=th0, b=b, c_mid=c_mid, St=2 * f * h0 / U,
                   T_mN=r["T_mean"] * 1e3, P_mW=r["P_mean"] * 1e3, CT=r["CT"], CP=r["CP"],
                   eta=r["eta"], AR_eff=r.get("AR", job["AR"]),
                   rate_pk=float(np.degrees(np.radians(th0) * 2 * math.pi * f)),
                   dt=time.time() - t0)
        return out
    except Exception as e:
        return dict(job, error=repr(e))


def run_jobs(jobs, out_csv, nproc=2):
    keys = None
    done = []
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(one, jobs, chunksize=2)):
            done.append(r)
            if i % 10 == 0:
                print(f"{i+1}/{len(jobs)}  " + (f"T={r.get('T_mN', float('nan')):.1f} mN η={r.get('eta', float('nan')):.2f}" if "T_mN" in r else f"ERR {r.get('error')}"), flush=True)
    import csv
    keys = sorted({k for r in done for k in r})
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in done: w.writerow(r)
    print("→", out_csv, len(done))


if __name__ == "__main__":
    stage = sys.argv[1]
    base = dict(family="rect", AR=3.0, taper=1.0, sweep=0.0, pivot=0.25)
    if stage == "A":       # 运动学：矩形 AR 3
        jobs = [dict(base, f=f, h0=h0, amax=a, psi=p)
                for f in (1.0, 1.5, 2.0, 2.5, 3.0) for h0 in (0.02, 0.03, 0.04, 0.05)
                for a in (15, 20, 25, 30) for p in (75, 90, 105)]
    elif stage == "B":     # 几何：在两组运动学上扫平面形
        kins = json.load(open(os.path.join(HERE, "sweep_B_kins.json")))
        jobs = []
        for kin in kins:
            for AR in (2, 3, 4, 5, 6):
                for fam, tapers, sweeps in (("rect", (1.0,), (0,)), ("fishtail", (0.3, 0.5, 0.7), (0,)),
                                            ("delta", (0.3, 0.5, 0.7), (0,)), ("ellipse", (1.0,), (0,)),
                                            ("lunate", (0.4, 0.6), (20, 35))):
                    for tp in tapers:
                        for sw in sweeps:
                            jobs.append(dict(kin, family=fam, AR=AR, taper=tp, sweep=sw, pivot=0.25))
    elif stage == "C":     # 在最优几何上重扫运动学 + 俯仰轴
        geo = json.load(open(os.path.join(HERE, "sweep_C_geo.json")))
        jobs = [dict(geo, f=f, h0=h0, amax=a, psi=p, pivot=pv)
                for f in (1.5, 2.0, 2.5, 3.0) for h0 in (0.03, 0.04, 0.05)
                for a in (20, 25, 30) for p in (75, 90, 105) for pv in (0.20, 0.30, 0.40)]
    else:
        sys.exit("stage A|B|C")
    print(f"stage {stage}: {len(jobs)} runs")
    run_jobs(jobs, os.path.join(HERE, f"sweep_{stage}.csv"))
