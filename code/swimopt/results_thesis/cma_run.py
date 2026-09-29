# Resumable CMA-ES runner for thesis experiments.
# usage: python cma_run.py TAG KIND MODE SEED BUDGET WALL_S [key=value overrides...]
#   KIND: toy | fluke ; MODE: full | preset:<diag|fb|lr|wave|inphase>
import sys, pickle, os
from common import *
import cma
tag, kind, mode, seed, budget, wall = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
ov = {}
for kv in sys.argv[7:]:
    k, v = kv.split("="); ov[k] = json.loads(v)
cfg = toy_cfg(**ov) if kind == "toy" else fluke_cfg(**ov)
sw = Swimmer(cfg["model"], cfg); g = sw.gait; B = g.blocks()
import foilvec; foilvec.install(sw.foil)

def to_full(z):
    if mode == "full":
        return g.expand(z)
    if mode == "repair":   # fluke: shift the mean shank angle so both stroke extremes stay inside [15.5 deg, 1.24 rad]
        x = g.expand(z).copy(); p = g.decode(x)
        amp = float(p["A"][0][0]); mean = float(p["offs"][0])
        lo, hi = np.radians(15.5) + amp, 1.24 - amp
        mean = float(np.clip(mean, lo, hi)) if lo <= hi else 0.5*(lo+hi)
        x[B["offset"][0]] = (mean - g.off_range[0])/(g.off_range[1]-g.off_range[0])
        return x
    ph1 = g.preset_phases(mode.split(":")[1])          # normalised [0,1]
    z = np.clip(np.asarray(z, float), 0, 1); x = np.full(g.dim, 0.5)
    x[B["freq"][0]] = z[0]; x[B["duty"][0]] = z[1]
    x[B["amp1"][0]:B["amp1"][1]] = z[2:5]; x[B["amp2"][0]:B["amp2"][1]] = z[5:8]
    x[B["phase1"][0]:B["phase1"][1]] = ph1
    x[B["phase2"][0]:B["phase2"][1]] = (z[8] + 2*ph1) % 1.0   # phi2_i = psi + 2 phi1_i : same feathering timing in every leg
    x[B["offset"][0]:B["offset"][1]] = z[9:12]
    return x
dim = g.dim_opt if mode in ("full","repair") else 12
x0 = g.x0_opt() if mode in ("full","repair") else np.full(12, 0.5)
st = os.path.join(HERE, f"run_{tag}.pkl")
if os.path.exists(st):
    S = pickle.load(open(st, "rb"))
else:
    es = cma.CMAEvolutionStrategy(x0, cfg.get("sigma0", 0.25), {"bounds": [0, 1], "popsize": cfg.get("popsize", 10),
                                  "seed": seed, "verbose": -9})
    S = dict(es=es, n=0, gens=[], evals=[], best=dict(fitness=-1e9), done=False, t=0.0)
if S["n"] < budget and not S["es"].stop(): S["done"] = False
t0 = time.time(); t_start = S["t"]; wall_budget = wall
es = S["es"]
while not S["done"]:
    X = es.ask(); F = []; rows = []
    for z in X:
        r = sw.rollout(to_full(z)); S["n"] += 1
        F.append(-r["fitness"])
        row = {k: r[k] for k in ("fitness", "speed", "yaw", "ok", "power", "pitch_rms", "heave_rms", "pitch_amp", "roll_amp", "bl_s", "servo_rate")}
        S["evals"].append(row); rows.append(row)
        if r["fitness"] > S["best"]["fitness"]:
            S["best"] = dict(row, z=list(map(float, z)), x=list(map(float, to_full(z))), eval=S["n"])
    es.tell(X, F)
    fs = [rw["fitness"] for rw in rows if rw["ok"]]
    S["gens"].append(dict(gen=len(S["gens"]) + 1, n=S["n"], sigma=float(es.sigma), best=S["best"]["fitness"],
                          best_speed=S["best"]["speed"], pop_median=float(np.median(fs)) if fs else None,
                          mean=list(map(float, es.mean))))
    if S["n"] >= budget or es.stop(): S["done"] = True
    S["t"] += time.time() - t0; t0 = time.time(); wall -= 0  # accumulate
    pickle.dump(S, open(st + ".tmp", "wb")); os.replace(st + ".tmp", st)
    t_used = S["t"] - t_start
    if t_used > wall_budget: break
save(f"run_{tag}.json", dict(tag=tag, kind=kind, mode=mode, seed=seed, overrides=ov, n=S["n"], done=S["done"],
                            wall_s=S["t"], best=S["best"], gens=S["gens"], evals=S["evals"]))
print(f"{tag}: n={S['n']} done={S['done']} best_speed={S['best'].get('speed',0):.4f} yaw={S['best'].get('yaw',0):.1f} wall={S['t']:.0f}s")
