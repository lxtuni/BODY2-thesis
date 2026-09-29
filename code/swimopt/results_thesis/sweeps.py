# Parameter sweeps for the thesis chapter. usage: python sweeps.py NAME
import sys, multiprocessing as mp
from common import *
import foilvec, mujoco
NAME = sys.argv[1]

def mk(kind, **kw):
    cfg = toy_cfg(**kw) if kind == "toy" else fluke_cfg(**kw)
    sw = Swimmer(cfg["model"], cfg); foilvec.install(sw.foil); return sw

def pick(r): return {k: r[k] for k in ("ok","fitness","speed","yaw","pitch_rms","heave_rms","pitch_amp","roll_amp","power","passive_pitch_max","bl_s")}

# ---------------- fluke helpers ----------------
FREF = dict(freq=1.4, qr_mean=np.radians(38.15), qr_amp=np.radians(23.06))   # exported firmware gait (results_fluke_bounded)
def fluke_x(sw, freq, mean, amp, preset="diag"):
    g = sw.gait; B = g.blocks(); x = g.base_x.copy()
    x[B["freq"][0]] = (freq-g.freq_range[0])/(g.freq_range[1]-g.freq_range[0])
    x[B["amp1"][0]] = (amp-g.amp_range[0])/(g.amp_range[1]-g.amp_range[0])
    x[B["offset"][0]] = (mean-g.off_range[0])/(g.off_range[1]-g.off_range[0])
    ph = {"diag":{"FL":0,"BR":0,"FR":.5,"BL":.5},"fb":{"FL":0,"FR":0,"BL":.5,"BR":.5},"lr":{"FL":0,"BL":0,"FR":.5,"BR":.5},
          "wave":{"FL":0,"FR":.25,"BR":.5,"BL":.75},"inphase":{"FL":0,"FR":0,"BL":0,"BR":0}}[preset]
    for i, nm in enumerate(g.names):
        x[B["phase1"][0]+i] = ph[next(L for L in ph if L in nm)]
    return x  # full vector; values may lie outside [0,1] only if outside configured ranges (expand clips)

def fluke_raw(freq, mean, amp, preset="diag", k_mNm=None, rigid=False, extra=None):
    """Evaluate with physical values; bypass the [0,1] range clipping by editing gait ranges."""
    sw = mk("fluke"); g = sw.gait
    g.freq_range = [0.0, 4.0]; g.amp_range = [0.0, 1.0]; g.off_range = [0.0, 1.5]
    m = sw.model
    for j in range(m.njnt):
        nm = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, j) or ""
        if nm.startswith("pitch_"):
            if k_mNm is not None: m.jnt_stiffness[j] = k_mNm*1e-3
            if rigid: m.jnt_range[j] = [-1e-4, 1e-4]
    r = sw.rollout(fluke_x(sw, freq, mean, amp, preset))
    out = pick(r); out.update(freq=freq, mean_deg=np.degrees(mean), amp_deg=np.degrees(amp), preset=preset, k_mNm=k_mNm, rigid=rigid)
    return out

def _fa(args): return fluke_raw(*args[:4], **args[4])

def run_pool(fn, jobs):
    with mp.Pool(2) as P: return P.map(fn, jobs)

if NAME == "fluke_map":        # frequency x shank amplitude, diagonal, mean 38 deg
    # mean shank angle raised where needed so the lower extreme stays at the 15 deg hull-clearance limit
    F = [0.6, 0.9, 1.2, 1.5, 1.8, 2.1]; A = np.radians([6, 12, 17, 23, 28])
    res = run_pool(_fa, [(f, max(FREF["qr_mean"], np.radians(15.5)+a), a, "diag", {}) for f in F for a in A])
    save("S_fluke_map.json", dict(F=F, A_deg=list(np.degrees(A)), res=res))
elif NAME == "fluke_k":        # passive-hinge stiffness
    K = [0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0]
    jobs = [(f, FREF["qr_mean"], FREF["qr_amp"], "diag", dict(k_mNm=k)) for f in (1.0, 1.4) for k in K]
    jobs += [(f, FREF["qr_mean"], FREF["qr_amp"], "diag", dict(rigid=True)) for f in (1.0, 1.4)]
    save("S_fluke_k.json", dict(K=K, res=run_pool(_fa, jobs)))
elif NAME == "fluke_phase":
    jobs = [(FREF["freq"], FREF["qr_mean"], FREF["qr_amp"], p, {}) for p in ("diag", "fb", "lr", "wave", "inphase")]
    save("S_fluke_phase.json", dict(res=run_pool(_fa, jobs)))
elif NAME == "toy_harm":       # 2nd-harmonic flipper amplitude x relative phase on the demo gait
    from view import demo_params
    def job(a):
        A2, psi = a; sw = mk("toy", w_yaw=1.0); g = sw.gait; B = g.blocks(); x = demo_params(g)
        x[B["amp2"][0]+2] = A2/g.amp_range[1]
        x[B["phase2"][0]:B["phase2"][1]] = (psi/360.0) % 1.0
        r = sw.rollout(x); o = pick(r); o.update(A2=A2, psi=psi); return o
    globals()["job"] = job
    A2 = [0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9]; PS = list(range(0, 360, 30))
    save("S_toy_harm.json", dict(A2=A2, PSI=PS, res=run_pool(job, [(a, p) for a in A2 for p in PS])))
