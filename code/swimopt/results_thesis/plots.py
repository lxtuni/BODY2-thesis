import json, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "..", "..", "thesis_chapter", "figures")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif", "font.size": 8.5,
                     "axes.labelsize": 8.5, "legend.fontsize": 7.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight", "lines.linewidth": 1.2})
C = ["#1f5f8b", "#d1495b", "#edae49", "#00798c", "#66a182", "#8d6a9f", "#555555"]
L = lambda n: json.load(open(os.path.join(H, n)))
def sv(f, n): f.savefig(os.path.join(OUT, n + ".pdf")); f.savefig(os.path.join(OUT, n + ".png"), dpi=200); plt.close(f)
W = 6.0

# ---------- foil coefficients ----------
def coeffs(a, AR=3.125, CLa=None, a_s=np.radians(15), CD0=0.02, e=0.9, CN=1.8):
    CLa = CLa or 2*np.pi*AR/(AR+2); aa = np.abs(a)
    sig = np.clip((1.6*a_s-aa)/(0.6*a_s), 0, 1); sc = np.sin(a)*np.cos(a); CLatt = CLa*sc
    return sig*CLatt+(1-sig)*CN*sc, sig*(CD0+CLatt**2/(np.pi*e*AR))+(1-sig)*CN*np.sin(a)**2, sig
a = np.radians(np.linspace(-90, 90, 721)); CL, CD, sg = coeffs(a)
f, ax = plt.subplots(1, 2, figsize=(W, 2.2), gridspec_kw=dict(width_ratios=[1.7, 1], wspace=0.35))
ax[0].plot(np.degrees(a), CL, c=C[0], label="$C_L$"); ax[0].plot(np.degrees(a), CD, c=C[1], label="$C_D$")
aa_ = np.where(np.abs(a) < np.radians(28), a, np.nan); ax[0].plot(np.degrees(aa_), 3.83*aa_, c=C[0], ls=":", lw=0.8, label="$C_{L\\alpha}\\,\\alpha$")
ax[0].axvspan(15, 24, color="0.9"); ax[0].axvspan(-24, -15, color="0.9")
ax[0].set_xlabel("angle of attack $\\alpha$ [deg]"); ax[0].set_ylabel("coefficient [-]"); ax[0].legend(frameon=False, loc="upper left"); ax[0].set_xlim(-90, 90)
ax[0].text(19.5, 1.55, "blend", ha="center", fontsize=6.5, color="0.35")
ax[1].plot(np.degrees(a), sg, c=C[3]); ax[1].set_xlabel("$\\alpha$ [deg]"); ax[1].set_ylabel("attached-flow weight $\\sigma$ [-]"); ax[1].set_xlim(-90, 90)
sv(f, "fig_foil_coeffs")

# ---------- model checks ----------
E = L("E_checks.json")
f, ax = plt.subplots(1, 2, figsize=(W, 2.2))
for k, c, lab in (("float_toy", C[0], "toy quadruped"), ("float_fluke", C[1], "BODY2 fluke")):
    z = np.array(E[k]["z_t"]); ax[0].plot(z[:, 0], z[:, 1]*1e3, c=c, label=lab)
ax[0].set_xlabel("time [s]"); ax[0].set_ylabel("trunk-origin height $z$ [mm]"); ax[0].legend(frameon=False); ax[0].set_xlim(0, 12)
co = E["coast"]; fh = (0-(E["float_fluke"]["z_eq"]-0.039))/0.05; T = np.array(co["T"]); V = np.array(co["V"])
ax[1].plot(T, V*100, c=C[0], label="simulated coast-down")
def integ(a_, b_, m_):
    v = [V[0]]; dt = T[1]-T[0]
    for _ in T[1:]: v.append(v[-1] - dt*(a_*v[-1]**2 + b_*v[-1])/m_)
    return np.array(v)
ax[1].plot(T, integ(co["a_fit"], co["b_fit"], co["m_eff"])*100, c=C[1], ls="--", label="fitted $m\\dot v=-av^2-bv$")
ax[1].plot(T, integ(co["hull_a_cfg"]*fh, co["hull_b_cfg"]*fh, co["m_eff"])*100, c=C[6], ls=":", label="hull coefficients only")
ax[1].set_xlabel("time [s]"); ax[1].set_ylabel("surge velocity $u$ [cm/s]"); ax[1].legend(frameon=False); ax[1].set_xlim(0, 12)
sv(f, "fig_model_checks")

# ---------- harmonic map ----------
S = L("S_toy_harm.json"); A2 = np.array(S["A2"]); PS = np.array(S["PSI"])
M = np.array([[next(r["speed"] for r in S["res"] if r["A2"] == aa and r["psi"] == pp) for pp in PS] for aa in A2])
f, ax = plt.subplots(1, 2, figsize=(W, 2.3), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.6))
vm = np.abs(M).max()*100
im = ax[0].pcolormesh(np.r_[PS, 360]-15, np.r_[A2, 1.05]-0.075, M*100, cmap="RdBu_r", vmin=-vm, vmax=vm, shading="flat")
cb = f.colorbar(im, ax=ax[0]); cb.set_label("forward speed [cm/s]")
ax[0].set_xlabel("2nd-harmonic phase $\\psi$ [deg]"); ax[0].set_ylabel("flipper $A_2$ [rad]"); ax[0].set_xticks(range(0, 361, 90))
for j, aa in enumerate(A2[[0, 2, 4, 6]]):
    i = list(A2).index(aa); ax[1].plot(PS, M[i]*100, "o-", ms=2.5, c=C[j], label=f"$A_2$={aa:.2f}")
ax[1].axhline(0, c="k", lw=0.5); ax[1].set_xlabel("$\\psi$ [deg]"); ax[1].set_ylabel("speed [cm/s]"); ax[1].legend(frameon=False, fontsize=6.5); ax[1].set_xticks(range(0, 361, 90))
sv(f, "fig_harmonic_map")

# ---------- convergence ----------
def curve(t):
    d = L(f"run_{t}.json"); s = [e["speed"] if e["ok"] else np.nan for e in d["evals"]]
    fit = np.array([e["fitness"] for e in d["evals"]]); best = np.maximum.accumulate(fit)
    # speed of the incumbent
    inc = []; bi = -1e9; bs = 0
    for e in d["evals"]:
        if e["fitness"] > bi: bi, bs = e["fitness"], e["speed"]
        inc.append(bs)
    return np.arange(1, len(fit)+1), best, np.array(inc), d["gens"]
f, ax = plt.subplots(1, 2, figsize=(W, 2.3))
for t, c, ls, lab in (("toy_s1", C[6], "-", "full 35-D, seed 1"), ("toy_s2", C[6], "--", "full 35-D, seed 2"),
                      ("g_diag", C[0], "-", "diagonal 12-D, seed 1"), ("g2_diag", C[0], "--", "diagonal 12-D, seed 2"),
                      ("g_fb", C[1], "-", "front–back 12-D, seed 1"), ("g2_fb", C[1], "--", "front–back 12-D, seed 2")):
    n, b, inc, G = curve(t); ax[0].plot(n, inc*100, c=c, ls=ls, label=lab)
    ax[1].plot([g["n"] for g in G], [g["sigma"] for g in G], c=c, ls=ls)
ax[0].set_xlabel("rollouts"); ax[0].set_ylabel("speed of incumbent [cm/s]"); ax[0].set_xlim(0, 500)
f.legend(*ax[0].get_legend_handles_labels(), frameon=False, fontsize=6.5, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.2))
ax[1].set_xlabel("rollouts"); ax[1].set_ylabel("CMA-ES step size $\\sigma$ [-]"); ax[1].set_ylim(0, 0.4)
sv(f, "fig_convergence")

# ---------- gait patterns (toy) ----------
P = ["diag", "fb", "inphase", "wave", "lr"]; PN = ["diagonal", "front–back", "in-phase", "wave", "left–right"]
def met(t, k): b = L(f"run_{t}.json")["best"]; return b[k]
f, ax = plt.subplots(1, 3, figsize=(W, 2.1), gridspec_kw=dict(wspace=0.5))
for j, (k, sc, lab) in enumerate((("speed", 100, "speed [cm/s]"), ("pitch_rms", 180/np.pi, "pitch std [deg]"), ("power", 1, "mean servo power [W]"))):
    v1 = [met(f"g_{p}", k)*sc for p in P]; v2 = [met(f"g2_{p}", k)*sc for p in P]
    x = np.arange(len(P)); ax[j].bar(x, (np.array(v1)+v2)/2, color=[C[i] for i in range(5)], alpha=0.55, width=0.65)
    ax[j].plot(x-0.12, v1, "k.", ms=4); ax[j].plot(x+0.12, v2, "kx", ms=3.5)
    ax[j].set_xticks(x); ax[j].set_xticklabels(PN, rotation=35, ha="right"); ax[j].set_ylabel(lab)
sv(f, "fig_gait_patterns")

# ---------- lift / drag split ----------
D = L("E_liftdrag_toy.json"); keys = ["g_diag", "g2_diag", "g_fb", "g2_fb", "g_inphase", "g2_inphase", "toy_s1", "toy_s2"]
lab = ["diag s1", "diag s2", "f–b s1", "f–b s2", "in-ph. s1", "in-ph. s2", "35-D s1", "35-D s2"]
f, ax = plt.subplots(1, 2, figsize=(W, 2.2), gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.45))
x = np.arange(len(keys)); lt = np.array([D[k]["lift_thrust_N"] for k in keys])*1e3; dt_ = np.array([D[k]["drag_thrust_N"] for k in keys])*1e3
ax[0].bar(x-0.18, lt, 0.36, color=C[0], label="flow-normal (lift) part"); ax[0].bar(x+0.18, dt_, 0.36, color=C[1], label="flow-aligned (drag) part")
ax[0].axhline(0, c="k", lw=0.5); ax[0].set_xticks(x); ax[0].set_xticklabels(lab, rotation=35, ha="right"); ax[0].set_ylabel("cycle-mean $F_x$, 4 flippers [mN]"); ax[0].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 1.12))
ax[1].scatter([D[k]["path_transverse_ratio"] for k in keys], [D[k]["alpha_thrust_weighted_deg"] for k in keys], c=[C[0]]*2+[C[1]]*2+[C[2]]*2+[C[6]]*2, s=14)
ax[1].set_xlabel("transverse path share [-]"); ax[1].set_ylabel("thrust-weighted incidence [deg]"); ax[1].set_xlim(0.3, 0.7); ax[1].set_ylim(50, 80)
sv(f, "fig_liftdrag")

# ---------- duty ----------
Dd = L("E_duty.json"); f, ax = plt.subplots(figsize=(3.2, 2.1))
for t, c, lab in (("g_diag", C[0], "diagonal optimum"), ("g2_fb", C[1], "front–back optimum")):
    R = Dd[t]["rows"]; du = [r["duty"] for r in R]; sp = [r["speed"]*100 if r["ok"] else np.nan for r in R]
    ax.plot(du, sp, "o-", c=c, ms=3, label=lab); bad = [r["duty"] for r in R if not r["ok"]]
    ax.plot(bad, [0.5]*len(bad), "x", c=c, ms=5); ax.axvline(Dd[t]["opt_duty"], c=c, ls=":", lw=0.8)
ax.set_xlabel("power-stroke duty cycle $d$ [-]"); ax.set_ylabel("speed [cm/s]"); ax.legend(frameon=False, fontsize=6.5)
sv(f, "fig_duty")

# ---------- fluke maps ----------
Mf = L("S_fluke_map.json"); F = Mf["F"]; A = Mf["A_deg"]
sp = np.array([[next(r["speed"] for r in Mf["res"] if r["freq"] == ff and abs(r["amp_deg"]-aa) < 1e-6) for ff in F] for aa in A])
pw = np.array([[next(r["power"] for r in Mf["res"] if r["freq"] == ff and abs(r["amp_deg"]-aa) < 1e-6) for ff in F] for aa in A])
pp = np.array([[next(r["passive_pitch_max"] for r in Mf["res"] if r["freq"] == ff and abs(r["amp_deg"]-aa) < 1e-6) for ff in F] for aa in A])
f, ax = plt.subplots(1, 2, figsize=(W, 2.4), gridspec_kw=dict(wspace=0.55))
cs = ax[0].contourf(F, A, sp*100, levels=12, cmap="viridis"); f.colorbar(cs, ax=ax[0]).set_label("speed [cm/s]")
c2 = ax[0].contour(F, A, pp, levels=[20, 30, 40], colors="w", linewidths=0.7); ax[0].clabel(c2, fmt="%d°", fontsize=6)
ax[0].plot(1.4, 23.06, "r*", ms=7); ax[0].annotate("exported firmware gait", (1.4, 23.06), (1.12, 8.0), color="w", fontsize=6.3, arrowprops=dict(arrowstyle="-", color="w", lw=0.6))
ax[0].plot(1.8, 28.0, "wo", ms=4, mfc="none")
ax[0].axhline(22.65, c="w", ls="--", lw=0.6); ax[0].text(0.62, 24.3, "mean raised to keep $q_r\\geq15°$", color="w", fontsize=6.3)
ax[0].set_xlabel("stroke frequency $f$ [Hz]"); ax[0].set_ylabel("shank amplitude $A$ [deg]")
for i, aa in enumerate(A):
    if aa < 10: continue
    ax[1].plot(pw[i]*1e3, sp[i]*100, "o-", ms=2.5, c=plt.cm.viridis(i/len(A)), label=f"$A$={aa:.0f}°")
ax[1].set_xscale("log"); ax[1].set_xlabel("mean servo power (4 legs) [mW]"); ax[1].set_ylabel("speed [cm/s]"); ax[1].legend(frameon=False, fontsize=6.5)
sv(f, "fig_fluke_map")

# ---------- stiffness ----------
K = L("S_fluke_k.json"); f, ax = plt.subplots(1, 2, figsize=(W, 2.1))
for fr, c in ((1.0, C[0]), (1.4, C[1])):
    R = [r for r in K["res"] if r["freq"] == fr and not r["rigid"]]; rg = [r for r in K["res"] if r["freq"] == fr and r["rigid"]][0]
    ax[0].semilogx([r["k_mNm"] for r in R], [r["speed"]*100 for r in R], "o-", c=c, ms=3, label=f"$f$={fr} Hz"); ax[0].axhline(rg["speed"]*100, c=c, ls=":", lw=0.8)
    ax[1].semilogx([r["k_mNm"] for r in R], [r["passive_pitch_max"] for r in R], "o-", c=c, ms=3)
ax[0].axvline(2.5, c="0.5", lw=0.6, ls="--"); ax[0].text(2.9, 1.2, "current\nleaf spring", fontsize=6.3, color="0.4")
ax[0].text(0.26, 0.15, "dotted: rigid fluke", fontsize=6.3, color="0.4")
ax[0].set_xlabel("hinge stiffness $k$ [mN·m/rad]"); ax[0].set_ylabel("speed [cm/s]"); ax[0].legend(frameon=False, fontsize=6.5)
ax[1].axhline(40, c="k", lw=0.6, ls="--"); ax[1].text(0.3, 41.5, "±40° stop", fontsize=6.3)
ax[1].set_xlabel("hinge stiffness $k$ [mN·m/rad]"); ax[1].set_ylabel("max passive pitch [deg]")
sv(f, "fig_fluke_stiffness")

# ---------- fluke phase ----------
Ph = L("S_fluke_phase.json")["res"]; names = {"diag": "diagonal", "fb": "front–back", "lr": "left–right", "wave": "wave", "inphase": "in-phase"}
f, ax = plt.subplots(1, 3, figsize=(W, 2.0), gridspec_kw=dict(wspace=0.5)); x = np.arange(len(Ph))
for j, (k, sc, lab) in enumerate((("speed", 100, "speed [cm/s]"), ("pitch_rms", 180/np.pi, "pitch std [deg]"), ("heave_rms", 1e3, "heave std [mm]"))):
    ax[j].bar(x, [r[k]*sc for r in Ph], color=[C[i] for i in range(5)], width=0.65); ax[j].set_xticks(x); ax[j].set_xticklabels([names[r["preset"]] for r in Ph], rotation=35, ha="right"); ax[j].set_ylabel(lab)
ax[0].axhline(0, c="k", lw=0.5)
sv(f, "fig_fluke_phase")

# ---------- fluke time history ----------
Hh = L("E_hist.json")["best"]; A_ = np.array(Hh["rows"]); f0 = Hh["f"]; t = A_[:, 0]
m = (t >= 6.0) & (t <= 6.0 + 2/f0); tt = (t[m]-6.0)*f0
f, ax = plt.subplots(4, 1, figsize=(W, 4.6), sharex=True)
ax[0].plot(tt, A_[m, 5], c=C[0], label="shank $q_r$ (FL)"); ax[0].plot(tt, A_[m, 9], c=C[1], label="passive fluke pitch $\\theta$ (FL)"); ax[0].set_ylabel("[deg]"); ax[0].legend(frameon=False, ncol=2)
ax[1].plot(tt, A_[m, 13], c=C[3]); ax[1].set_ylabel("$\\alpha$ FL [deg]"); ax[1].set_ylim(-100, 100); ax[1].axhspan(-15, 15, color="0.92")
ax[2].plot(tt, A_[m, 21]*1e3, c=C[0], label="$F_x$"); ax[2].plot(tt, A_[m, 25]*1e3, c=C[1], label="$F_z$"); ax[2].axhline(0, c="k", lw=0.4); ax[2].set_ylabel("fluke FL [mN]"); ax[2].legend(frameon=False, ncol=2)
ax[3].plot(tt, A_[m, 2]*100, c=C[0]); ax[3].set_ylabel("surge $u$ [cm/s]"); ax[3].set_xlabel("stroke cycles [-]")
sv(f, "fig_fluke_history")
# stats for text
half = t > 4.0
stats = dict(Fx_mean_mN=list(A_[half, 21:25].mean(0)*1e3), Fz_mean_mN=list(A_[half, 25:29].mean(0)*1e3),
             attached_frac=float((A_[half, 29:33] >= 0.999).mean()), stalled_frac=float((A_[half, 29:33] <= 1e-9).mean()),
             pitch_amp=float(np.ptp(A_[half, 9])/2), u_mean=float(A_[half, 2].mean()))
Hr = np.array(L("E_hist.json")["ref"]["rows"]); hr = Hr[:, 0] > 4.0
stats["ref"] = dict(Fx_mean_mN=list(Hr[hr, 21:25].mean(0)*1e3), attached_frac=float((Hr[hr, 29:33] >= 0.999).mean()), stalled_frac=float((Hr[hr, 29:33] <= 1e-9).mean()), u_mean=float(Hr[hr, 2].mean()))
json.dump(stats, open(os.path.join(H, "E_hist_stats.json"), "w"), indent=1)
print(json.dumps(stats, indent=1))
