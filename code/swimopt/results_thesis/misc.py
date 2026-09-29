# Remaining analyses: E1 model checks, lift/drag split, duty sweep, fluke time history
import sys
from common import *
import foilvec, mujoco
from kin import decompose
W = sys.argv[1]
def bestx(t): return np.array(json.load(open(os.path.join(HERE, f"run_{t}.json")))["best"]["x"])

if W == "liftdrag":
    cfg = toy_cfg(w_yaw=1.0); sw = Swimmer(cfg["model"], cfg)
    from view import demo_params
    out = {"demo": decompose(sw, demo_params(sw.gait))}
    for t in ("toy_s1", "toy_s2", "g_diag", "g2_diag", "g_fb", "g2_fb", "g_inphase", "g2_inphase"):
        out[t] = decompose(sw, bestx(t))
    out = {k: {a: float(b) for a, b in v.items()} for k, v in out.items()}
    save("E_liftdrag_toy.json", out)
    for k, v in out.items(): print(k, {a: round(b, 3) for a, b in v.items()})

elif W == "duty":
    cfg = toy_cfg(w_yaw=1.0); sw = Swimmer(cfg["model"], cfg); B = sw.gait.blocks()
    D = np.linspace(0, 1, 11); res = {}
    for t in ("g_diag", "g2_fb"):
        x = bestx(t); rows = []
        for dd in D:
            xx = x.copy(); xx[B["duty"][0]] = dd; r = sw.rollout(xx)
            rows.append(dict(duty=sw.gait.duty_range[0] + dd*(sw.gait.duty_range[1]-sw.gait.duty_range[0]), speed=r["speed"], power=r["power"], yaw=r["yaw"], ok=r["ok"]))
        res[t] = dict(opt_duty=float(sw.gait.decode(x)["duty"]), rows=rows)
        print(t, res[t]["opt_duty"], [(round(r["duty"], 2), round(r["speed"], 4)) for r in rows])
    save("E_duty.json", res)

elif W == "checks":
    out = {}
    # (a) vectorised vs loop hydro, determinism
    cfg = toy_cfg(); sw = Swimmer(cfg["model"], cfg); m, d = sw.model, sw.data
    from view import demo_params
    x = demo_params(sw.gait); p = sw.gait.decode(x); sw.reset_posture(p)
    mx = 0.0
    for k in range(1500):
        d.ctrl[:] = sw.gait.ctrl(p, k*0.002)
        sw.hydro.apply(d); A = d.xfrc_applied.copy(); sw.hydro._apply_slow(d); Bf = d.xfrc_applied.copy()
        mx = max(mx, float(np.max(np.abs(A-Bf)))); d.xfrc_applied[:] = A; mujoco.mj_step(m, d)
    t0 = time.time(); [sw.hydro.apply(d) for _ in range(2000)]; tv = (time.time()-t0)/2000
    t0 = time.time(); [sw.hydro._apply_slow(d) for _ in range(300)]; ts = (time.time()-t0)/300
    r1 = sw.rollout(x); r2 = sw.rollout(x)
    out["vec_vs_loop_maxdiff_N"] = mx; out["t_vec_ms"] = tv*1e3; out["t_loop_ms"] = ts*1e3
    out["deterministic"] = bool(r1["speed"] == r2["speed"] and r1["energy"] == r2["energy"])
    # (b) static flotation: buoyancy vs weight after settling, both robots, legs at rest posture
    for kind in ("toy", "fluke"):
        cfg = toy_cfg() if kind == "toy" else fluke_cfg()
        sw = Swimmer(cfg["model"], cfg); foilvec.install(sw.foil); m, d = sw.model, sw.data
        x = sw.gait.x0_opt(); pp = sw.gait.decode(sw.gait.expand(x)); pp["A"] = [a*0 for a in pp["A"]]
        sw.reset_posture(pp); zs = []; Fz = []
        for k in range(int(12.0/0.002)):
            d.ctrl[:] = sw.gait.ctrl(pp, 0.0); sw.hydro.apply(d); sw.foil.apply(d)
            if k*0.002 > 8.0: Fz.append(float(d.xfrc_applied[:, 2].sum()))
            if k % 10 == 0: zs.append([k*0.002, float(d.xpos[sw.trunk_id, 2]), float(np.degrees(np.arcsin(np.clip(-d.xmat[sw.trunk_id].reshape(3,3)[2,0], -1, 1))))])
            mujoco.mj_step(m, d)
        ma_all = float(np.sum(sw.hydro.ma_body)) + sum(i["ma"] for i in sw.foil.items)
        mass_true = float(np.sum(m.body_mass[1:])) - ma_all
        buoy = float(np.mean(Fz)) - ma_all*9.81          # remove the constant added-mass weight compensation
        zz = np.array(zs)
        out[f"float_{kind}"] = dict(mass_kg=mass_true, weight_N=mass_true*9.81, buoyancy_mean_N=buoy,
                                    z_eq=float(zz[zz[:,0] > 8, 1].mean()), z_p2p_late=float(np.ptp(zz[zz[:,0] > 8, 1])),
                                    pitch_eq_deg=float(zz[zz[:,0] > 8, 2].mean()), z_t=zs)
    # (c) virtual coast-down on the fluke robot (legs held at 40 deg): identify k = 0.5 rho CdA_eff / m_eff and c
    cfg = fluke_cfg(); sw = Swimmer(cfg["model"], cfg); foilvec.install(sw.foil); m, d = sw.model, sw.data
    pp = sw.gait.decode(sw.gait.expand(sw.gait.x0_opt())); pp["A"] = [a*0 for a in pp["A"]]
    sw.reset_posture(pp)
    for k in range(int(3.0/0.002)):
        sw.hydro.apply(d); sw.foil.apply(d); mujoco.mj_step(m, d)
    d.qvel[0] = 0.25; T = []; V = []
    for k in range(int(12.0/0.002)):
        sw.hydro.apply(d); sw.foil.apply(d); mujoco.mj_step(m, d)
        if k % 5 == 0: T.append(k*0.002); V.append(float(d.qvel[0]))
    T = np.array(T); V = np.array(V)
    # fit m_eff dv/dt = -a v^2 - b v via least squares on the derivative (smoothed)
    dv = np.gradient(V, T); msk = (T > 0.3) & (V > 0.01)
    Amat = np.c_[-V[msk]**2, -V[msk]]; coef, *_ = np.linalg.lstsq(Amat, dv[msk], rcond=None)
    m_eff = float(np.sum(m.body_mass[1:]))       # includes added mass (isotropic trick)
    a_fit, b_fit = coef*m_eff
    # expected hull-only coefficients from config: 0.5*rho*Cd*A*f and cv*f (f = immersion of hull geom)
    hull = [it for it in sw.hydro.items if it["name"] == "trunk_hull"][0]
    out["coast"] = dict(T=T.tolist(), V=V.tolist(), m_eff=m_eff, a_fit=float(a_fit), b_fit=float(b_fit),
                        hull_a_cfg=0.5*1000*hull["cd"][0]*hull["A"][0], hull_b_cfg=hull["cv"])
    save("E_checks.json", out)
    print({k: v for k, v in out.items() if not isinstance(v, dict)})
    for k in ("float_toy", "float_fluke"): print(k, {a: b for a, b in out[k].items() if a != "z_t"})
    print({a: b for a, b in out["coast"].items() if a not in ("T", "V")})

elif W == "hist":     # fluke time histories: reference gait and best feasible grid gait
    res = {}
    for tag, (f, mean, amp) in {"ref": (1.4, np.radians(38.15), np.radians(23.06)), "best": (1.8, np.radians(43.5), np.radians(28.0))}.items():
        cfg = fluke_cfg(); sw = Swimmer(cfg["model"], cfg); foilvec.install(sw.foil); g = sw.gait
        g.freq_range = [0, 4]; g.amp_range = [0, 1]; g.off_range = [0, 1.5]
        from sweeps import fluke_x
        x = fluke_x(sw, f, mean, amp); p = g.decode(g.expand(x)); m, d = sw.model, sw.data
        sw.reset_posture(p)
        for _ in range(int(sw.settle/0.002)): sw.hydro.apply(d); sw.foil.apply(d); mujoco.mj_step(m, d)
        pq = [int(m.jnt_qposadr[j]) for j in range(m.njnt) if (mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, j) or "").startswith("pitch_")]
        sq = [int(m.jnt_qposadr[int(m.actuator_trnid[a, 0])]) for a in range(m.nu)]
        rows = []
        for k in range(int(8.0/0.002)):
            t = k*0.002; d.ctrl[:] = g.ctrl(p, t); sw.hydro.apply(d); sw.foil.apply(d)
            lv = sw.foil.last_vec
            R = d.xmat[sw.trunk_id].reshape(3, 3)
            rows.append([t, float(d.xpos[sw.trunk_id, 0]), float(d.qvel[0]), float(d.xpos[sw.trunk_id, 2]),
                         float(np.degrees(np.arcsin(np.clip(-R[2, 0], -1, 1))))] +
                        [float(np.degrees(d.qpos[i])) for i in sq] + [float(np.degrees(d.qpos[i])) for i in pq] +
                        list(map(float, np.degrees(lv["alpha"]))) + list(map(float, lv["CL"])) + list(map(float, lv["Fh"][:, 0])) + list(map(float, lv["Fh"][:, 2])) + list(map(float, lv["sig"])))
            mujoco.mj_step(m, d)
        res[tag] = dict(f=f, mean_deg=float(np.degrees(mean)), amp_deg=float(np.degrees(amp)), names=g.names,
                        cols="t x vx z pitch shank4 flukepitch4 alpha4 CL4 Fx4 Fz4 sigma4", rows=rows[::2])
        A = np.array(rows); half = A[:, 0] > 4.0
        print(tag, "mean Fx per fluke (mN):", np.round(A[half, 25:29].mean(0)*1e3, 2), "alpha range", np.round(A[half, 13:17].min(), 1), np.round(A[half, 13:17].max(), 1),
              "stall frac", np.round((A[half, 33:37] < 1).mean(), 3))
    save("E_hist.json", res)
