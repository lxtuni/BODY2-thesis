#!/usr/bin/env python3
"""
用已有的 open / closed 全周期 CFD 算例，扫“收扇时刻 τ_c”与“开扇时刻 τ_o”（两态合成，不需要新算例）
    python3 fold_timing_sweep.py [--out .]
合成规则：τ ∈ [τ_o, τ_c) 用 open 算例的力，其余用 closed 算例（与 analyze_verify.py 的两态合成一致，只是切换时刻可变）。
注意：这是一阶灵敏度——没有模拟折扇过渡本身，open/closed 两条时间历程的尾流历史也各自独立。
输出 fig_fold_timing_sweep.png、fold_timing_sweep.json
"""
import os, sys, json, glob, argparse, importlib.util
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
HERE = os.path.dirname(os.path.abspath(__file__)); LD = os.path.normpath(os.path.join(HERE, "..", ".."))
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=HERE); ap.add_argument("--baseline", default=os.path.join(LD, "results_leg"))
ap.add_argument("--extra", nargs="*", default=[], help="追加算例目录 NAME=DIR（如 opt=~/run/opt，目录里要有 force_open/closed.dat、moment_*.dat、gait_*.json）"); a = ap.parse_args()
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
COL = {"current": "#8a8f96", "V1": "#eb6834", "V2": "#7a5cc0", "V3": "#2a78d6", "V3c12": "#16a085", "V3c06": "#7fb800", "V3trap": "#c2185b"}
NAME_EXTRA = {"V3c12": "V3c12 收拢 12 mm", "V3c06": "V3c06 羽化 5.6 mm", "V3trap": "V3trap 梯形律 + 12 mm"}
NAME = {"current": "现步态 136→76°", "V1": "V1 全行程", "V2": "V2", "V3": "V3 120→70°"}

def load(R):
    def rd(tag):
        d = np.loadtxt(os.path.join(R, f"force_{tag}.dat"), comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
        m = np.loadtxt(os.path.join(R, f"moment_{tag}.dat"), comments="#"); m = m[np.argsort(m[:, 0])]; _, i = np.unique(m[:, 0], return_index=True); m = m[i]
        return dict(t=d[:, 0], F=d[:, 1:4], tm=m[:, 0], M=m[:, 1:4], gait=json.load(open(os.path.join(R, f"gait_{tag}.json"))))
    c = {k: rd(k) for k in ("open", "closed")}; g = c["open"]["gait"]; T, DUTY = g["T"], g["DUTY"]
    ncyc = int(np.floor(min(x["t"][-1] for x in c.values()) / T + 1e-6)); N = 2000
    t0 = (ncyc - 1) * T; tt = np.linspace(t0, t0 + T, N, endpoint=False); tau = (tt - t0) / T
    F = {k: np.stack([np.interp(tt, c[k]["t"], c[k]["F"][:, i]) for i in range(3)], 1) for k in c}
    M = {k: np.stack([np.interp(tt, c[k]["tm"], c[k]["M"][:, i]) for i in range(3)], 1) for k in c}
    for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
        if k in g: os.environ[k] = str(g[k])
    if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
    sp = importlib.util.spec_from_file_location("ld_" + os.path.basename(R), os.path.join(LD, "leg_dynamics.py")); ld = importlib.util.module_from_spec(sp); sp.loader.exec_module(ld)
    psi, phi, _ = ld.gait(tt); PP = [ld.pose(p, q) for p, q in zip(psi, phi)]; E = np.array([p["E"] for p in PP])
    vE = np.gradient(E, tt, axis=0); om = -np.gradient(np.unwrap(psi), tt); rE = E - np.array(ld.O2)
    def comp(to, tc):
        """open 窗口 [to, tc)（按周期取模），其余 closed；返回周期平均 Fx, Fz [mN] 与功率 P [mW]"""
        w = ((tau - to) % 1.0) < ((tc - to) % 1.0)
        Fc = np.where(w[:, None], F["open"], F["closed"]); Mc = np.where(w[:, None], M["open"], M["closed"])
        ME = Mc[:, 1] - (rE[:, 1] * Fc[:, 0] - rE[:, 0] * Fc[:, 2]); P = -(Fc[:, 0] * vE[:, 0] + Fc[:, 2] * vE[:, 1] + ME * om)
        return Fc[:, 0].mean() * 1e3, Fc[:, 2].mean() * 1e3, P.mean() * 1e3
    return dict(T=T, DUTY=DUTY, ncyc=ncyc, tau=tau, comp=comp, Fx_open=F["open"][:, 0] * 1e3, Fx_closed=F["closed"][:, 0] * 1e3,
                tau_open=float(g.get("TAU_OPEN", 0.0)), tau_close=float(g.get("TAU_CLOSE", DUTY)))   # --traj 算例：gait.json 里带开合相位

runs = {"current": a.baseline, **{k: os.path.join(LD, "results_verify", k) for k in ("V1", "V2", "V3")}}
for kv in a.extra:
    k, v = kv.split("=", 1); runs[k] = os.path.expanduser(v); COL.setdefault(k, ["#16a085", "#7fb800", "#c2185b", "#f39c12"][len(runs) % 4]); NAME.setdefault(k, NAME_EXTRA.get(k, k))
tc_grid = np.round(np.arange(0.25, 0.701, 0.01), 3); to_grid = np.round(np.arange(-0.15, 0.151, 0.01), 3)
res = {}
for k, R in runs.items():
    if not os.path.exists(os.path.join(R, "force_open.dat")): continue
    d = load(R); D = d["tau_close"]; TO = d["tau_open"]
    sc = np.array([d["comp"](TO, tc) for tc in tc_grid]); so = np.array([d["comp"](to, D) for to in to_grid])
    base = d["comp"](TO, D)
    ic = int(np.argmax(sc[:, 0])); io = int(np.argmax(so[:, 0]))
    # 划水相内 open 与 closed 的推力差为零的时刻（即“合上也不吃亏”的最早时刻）
    tau = d["tau"]; dF = d["Fx_open"] - d["Fx_closed"]; k5 = max(1, int(0.02 * d["T"] / (d["T"] / len(tau)))); dFs = np.convolve(dF, np.ones(k5) / k5, "same")
    after = (tau > 0.5 * D) & (tau < 0.95); cross = tau[after][np.argmax(dFs[after] <= 0)] if np.any(dFs[after] <= 0) else None
    res[k] = dict(DUTY=D, tau_open=TO, T=d["T"], base=dict(Fx=base[0], Fz=base[1], P=base[2]), tc=tc_grid.tolist(), Fx_tc=sc[:, 0].tolist(), Fz_tc=sc[:, 1].tolist(), P_tc=sc[:, 2].tolist(),
                  to=to_grid.tolist(), Fx_to=so[:, 0].tolist(), P_to=so[:, 2].tolist(), best_tc=float(tc_grid[ic]), best_Fx_tc=float(sc[ic, 0]), best_to=float(to_grid[io]), best_Fx_to=float(so[io, 0]),
                  tau_dF0=None if cross is None else float(cross), dFx_per_0p05_early=float(np.interp(D - 0.05, tc_grid, sc[:, 0]) - base[0]), dFx_per_0p05_late=float(np.interp(D + 0.05, tc_grid, sc[:, 0]) - base[0]))
    print(f"{k:8s} DUTY {D:.2f}: base Fx {base[0]:6.1f} mN, P {base[2]:5.1f} mW | best τ_c {tc_grid[ic]:.2f} → {sc[ic,0]:6.1f} mN | early −0.05: {res[k]['dFx_per_0p05_early']:+5.1f}, late +0.05: {res[k]['dFx_per_0p05_late']:+5.1f} mN | best τ_o {to_grid[io]:+.2f} → {so[io,0]:6.1f} mN | open−closed 推力差过零 τ ≈ {cross}")
json.dump(res, open(os.path.join(a.out, "fold_timing_sweep.json"), "w"), indent=1, ensure_ascii=False)

fig, axs = plt.subplots(1, 3, figsize=(16, 5.4), facecolor=SURF, gridspec_kw=dict(wspace=0.3, left=0.055, right=0.985, top=0.80, bottom=0.14))
def style(ax):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9.5)
ax = axs[0]; style(ax)
for k, r in res.items():
    ax.plot(r["tc"], r["Fx_tc"], color=COL[k], lw=2.2, zorder=3, label=NAME[k]); ax.scatter([r["DUTY"]], [r["base"]["Fx"]], s=60, color=COL[k], zorder=5, edgecolors=SURF, linewidths=1.4)
    ax.scatter([r["best_tc"]], [r["best_Fx_tc"]], s=60, marker="D", color=COL[k], zorder=5, edgecolors=SURF, linewidths=1.4)
ax.set_xlabel("收扇时刻 τ_c（开扇固定在 τ = 0）", color=INK2, fontsize=10); ax.set_ylabel("周期平均推力 F̄x [mN]", color=INK2, fontsize=10)
ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)
ax.set_title("A  收扇：比划水结束（DUTY）早 ≈ 0.1 T 最好（+2–5 %），晚收代价大\n     ● 现用切换点 = DUTY，◆ 扫描最优", loc="left", fontsize=11, color=INK, pad=9)
ax = axs[1]; style(ax)
for k, r in res.items():
    ax.plot(r["to"], r["Fx_to"], color=COL[k], lw=2.2, zorder=3, label=NAME[k]); ax.scatter([r["tau_open"]], [r["base"]["Fx"]], s=60, color=COL[k], zorder=5, edgecolors=SURF, linewidths=1.4)
ax.axvline(0, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2)
ax.set_xlabel("开扇时刻 τ_o（负 = 回收末段提前张开；收扇固定在 DUTY）", color=INK2, fontsize=10); ax.set_ylabel("周期平均推力 F̄x [mN]", color=INK2, fontsize=10)
ax.set_title("B  开扇：在回收末段速度过零前 ≈ 0.05 T 张开最好（+3–4 %）\n     晚开 0.05 T 损失 3–4 %，晚 0.15 T 损失 ≈ 20 % → 宁早勿晚", loc="left", fontsize=11, color=INK, pad=9)
ax = axs[2]; style(ax)
for k, r in res.items():
    Fx, P = np.array(r["Fx_tc"]), np.array(r["P_tc"]); ok = P > 0.5
    ax.plot(np.array(r["tc"])[ok], (Fx / P)[ok], color=COL[k], lw=2.2, zorder=3, label=NAME[k])
    ax.scatter([r["DUTY"]], [r["base"]["Fx"] / r["base"]["P"]], s=60, color=COL[k], zorder=5, edgecolors=SURF, linewidths=1.4)
ax.set_xlabel("收扇时刻 τ_c", color=INK2, fontsize=10); ax.set_ylabel("推力 / 机构输入功率 F̄x / P̄ [mN/mW = N/W]", color=INK2, fontsize=10)
ax.set_title("C  单位功率推力：峰值也在提前收扇一侧\n     （U = 0 系泊工况，只能比相对值）", loc="left", fontsize=11, color=INK, pad=9)
fig.suptitle("用现有 open/closed 算例扫描脚蹼开合时刻（两态合成，第 3 周期，未模拟折扇过渡本身）", fontsize=13.5, x=0.055, ha="left", y=0.965, color=INK)
fig.text(0.055, 0.905, "合成：τ ∈ [τ_o, τ_c) 取 open 算例的力，其余取 closed；功率 P = −(F·v_E + M_E·ω)，与 analyze_verify.py 相同。数据：results_leg（现步态）、results_verify/V1–V3", fontsize=9, color=MUTED, ha="left")
fn = os.path.join(a.out, "fig_fold_timing_sweep.png"); plt.savefig(fn, dpi=120, facecolor=SURF); print("→", fn)
