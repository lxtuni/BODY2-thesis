#!/usr/bin/env python3
"""
步态验证算例深入分析（在 compare_verify.py 之上）：时间历程、逐周期收敛、两对腿叠加后的垂向力残差（船体升沉激励）、权衡图、与模型前沿的对照
    python3 analyze_verify.py [--root results_verify] [--baseline results] [--pareto results/gait_pareto.json] [--out results_verify]
输出：verify_analysis.json、fig_verify_time.png、fig_verify_tradeoff.png、fig_verify_cycles.png、verify_table.md
"""
import os, sys, json, glob, argparse, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
cv = importlib.import_module("compare_verify")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"
COL = {"current": INK, "V1": GREEN, "V2": ORANGE, "V3": BLUE}
LAB = {"current": "现步态", "V1": "V1 前沿点", "V2": "V2 狗式占空比", "V3": "V3 低中心角"}

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=os.path.join(HERE, "results_verify")); ap.add_argument("--baseline", default=os.path.join(HERE, "results"))
ap.add_argument("--pareto", default=None); ap.add_argument("--out", default=None)
a = ap.parse_args(); OUT = a.out or a.root; os.makedirs(OUT, exist_ok=True)
pareto_fn = a.pareto or next((p for p in (os.path.join(a.baseline, "gait_pareto.json"), os.path.join(HERE, "results_leg", "gait_pareto.json")) if os.path.exists(p)), None)

def load_full(R):
    """读 open/closed，合成两态，返回细时间历程（最后一个完整周期）与逐周期均值"""
    def load(tag):
        fn = os.path.join(R, f"force_{tag}.dat")
        if not os.path.exists(fn): return None
        d = np.loadtxt(fn, comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
        m = np.loadtxt(os.path.join(R, f"moment_{tag}.dat"), comments="#"); m = m[np.argsort(m[:, 0])]; _, i = np.unique(m[:, 0], return_index=True); m = m[i]
        return dict(t=d[:, 0], F=d[:, 1:4], M=m[:, 1:4], tm=m[:, 0], gait=json.load(open(os.path.join(R, f"gait_{tag}.json"))))
    c = {k: load(k) for k in ("open", "closed")}
    if c["open"] is None or c["closed"] is None: return None
    g = c["open"]["gait"]; T = g["T"]; DUTY = g["DUTY"]; ncyc = int(np.floor(min(x["t"][-1] for x in c.values()) / T + 1e-6))
    N = 2000
    def cyc(k):
        t0 = (k - 1) * T; tt = np.linspace(t0, t0 + T, N, endpoint=False); tau = (tt - t0) / T
        F = {tag: np.stack([np.interp(tt, c[tag]["t"], c[tag]["F"][:, i]) for i in range(3)], 1) for tag in c}
        M = {tag: np.stack([np.interp(tt, c[tag]["tm"], c[tag]["M"][:, i]) for i in range(3)], 1) for tag in c}
        fo = (tau < DUTY)[:, None]; Fc = np.where(fo, F["open"], F["closed"]); Mc = np.where(fo, M["open"], M["closed"])
        return tt, tau, Fc, Mc, F, M
    per_cycle = []
    for k in range(1, ncyc + 1):
        tt, tau, Fc, Mc, F, M = cyc(k); per_cycle.append(dict(k=k, Fx=float(Fc[:, 0].mean() * 1e3), Fz=float(Fc[:, 2].mean() * 1e3), Fx_open=float(F["open"][:, 0].mean() * 1e3), Fx_closed=float(F["closed"][:, 0].mean() * 1e3)))
    tt, tau, Fc, Mc, F, M = cyc(ncyc)
    # 功率（运动学）
    for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
        if k in g: os.environ[k] = str(g[k])
    if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
    sp = importlib.util.spec_from_file_location("ld_" + os.path.basename(R), os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(sp); sp.loader.exec_module(ld)
    psi, phi, _ = ld.gait(tt); PP = [ld.pose(p, q) for p, q in zip(psi, phi)]; E = np.array([p["E"] for p in PP])
    vE = np.gradient(E, tt, axis=0); om = -np.gradient(np.unwrap(psi), tt); rE = E - np.array(ld.O2)
    ME = Mc[:, 1] - (rE[:, 1] * Fc[:, 0] - rE[:, 0] * Fc[:, 2]); P = -(Fc[:, 0] * vE[:, 0] + Fc[:, 2] * vE[:, 1] + ME * om)
    dt = T / N; k = max(1, int(round(0.02 / dt))); ker = np.ones(k) / k
    sm = lambda y: np.convolve(np.concatenate([y[-k:], y, y[:k]]), ker, "same")[k:-k]
    # 两对腿半周期错开后的合成（Cui 的升力残差；机器人靠浮力支撑，F_z 只是桨板力的无用垂向分量）：F_pair(τ) = F(τ) + F(τ + 1/2)
    half = N // 2; Fz_pair = Fc[:, 2] + np.roll(Fc[:, 2], -half); Fx_pair = Fc[:, 0] + np.roll(Fc[:, 0], -half)
    pw = tau < DUTY; Fx_pw, Fz_pw = Fc[pw, 0].mean(), Fc[pw, 2].mean()          # 划水相平均力矢量：大小与倾角（相对水平）
    return dict(T=T, DUTY=DUTY, ncyc=ncyc, tau=tau, Fx=sm(Fc[:, 0]) * 1e3, Fz=sm(Fc[:, 2]) * 1e3, My=sm(Mc[:, 1]) * 1e3, P=sm(P) * 1e3, psi=np.degrees(psi),
                Fx_open=sm(F["open"][:, 0]) * 1e3, Fx_closed=sm(F["closed"][:, 0]) * 1e3,
                Fx_mean=float(Fc[:, 0].mean() * 1e3), Fz_mean=float(Fc[:, 2].mean() * 1e3), P_mean=float(P.mean() * 1e3), My_peak=float(np.abs(sm(Mc[:, 1])).max() * 1e3),
                Fx_peak=float(sm(Fc[:, 0]).max() * 1e3), Fx_min=float(sm(Fc[:, 0]).min() * 1e3), Fz_peak=float(sm(Fc[:, 2]).max() * 1e3), Fz_min=float(sm(Fc[:, 2]).min() * 1e3), P_peak=float(sm(P).max() * 1e3),
                Fz_pair_std=float(np.std(sm(Fz_pair)) * 1e3), Fz_pair_pp=float((sm(Fz_pair).max() - sm(Fz_pair).min()) * 1e3), Fx_pair_std=float(np.std(sm(Fx_pair)) * 1e3),
                Fz_std=float(np.std(sm(Fc[:, 2])) * 1e3), per_cycle=per_cycle, psi_range=[float(np.degrees(psi).min()), float(np.degrees(psi).max())], f=1 / T,
                wmax_deg_s=float(np.degrees(np.abs(om)).max()),
                Fx_pw=float(Fx_pw * 1e3), Fz_pw=float(Fz_pw * 1e3), Fmag_pw=float(np.hypot(Fx_pw, Fz_pw) * 1e3), tilt_pw_deg=float(np.degrees(np.arctan2(Fz_pw, Fx_pw))),
                Fx_rec=float(Fc[~pw, 0].mean() * 1e3), Fz_rec=float(Fc[~pw, 2].mean() * 1e3))

S = {}
if os.path.exists(os.path.join(a.baseline, "force_open.dat")): S["current"] = load_full(a.baseline)
for d in sorted(glob.glob(os.path.join(a.root, "*"))):
    if os.path.isdir(d) and os.path.exists(os.path.join(d, "force_open.dat")):
        r = load_full(d)
        if r: S[os.path.basename(d)] = r
keys = list(S.keys()); MODEL = cv.MODEL
ref = S.get("current")
rows = ["| 方案 | T [s] / f [Hz] | DUTY | ψ 范围 | 峰值髋角速度 [°/s] | F̄_x [mN] | 模型 | 偏差 | F̄_z | F̄_z/F̄_x | F_x 峰 / 谷 | F_z 峰 / 谷 | |M_y| 峰 [mN·m] | P̄ [mW] | P 峰 | F̄_x/P̄ | 两对腿叠加 F_z 峰峰 | 相对现步态 |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for k in keys:
    r = S[k]; m = MODEL.get(k, {})
    rows.append(f"| {LAB.get(k,k)} | {r['T']:.2f} / {r['f']:.2f} | {r['DUTY']:.2f} | {r['psi_range'][1]:.0f}→{r['psi_range'][0]:.0f}° | {r['wmax_deg_s']:.0f} | **{r['Fx_mean']:.1f}** | {m.get('F0','—')} | {('%+.0f %%' % (100*(r['Fx_mean']/m['F0']-1))) if m.get('F0') else '—'} | {r['Fz_mean']:.1f} | {r['Fz_mean']/r['Fx_mean']:.2f} | {r['Fx_peak']:.0f} / {r['Fx_min']:.0f} | {r['Fz_peak']:.0f} / {r['Fz_min']:.0f} | {r['My_peak']:.1f} | {r['P_mean']:.1f} | {r['P_peak']:.0f} | {r['Fx_mean']/r['P_mean']:.2f} | {r['Fz_pair_pp']:.0f} | {('%+.0f %%' % (100*(r['Fx_mean']/ref['Fx_mean']-1))) if (ref and k!='current') else '—'} |")
cyc_rows = ["", "| 方案 | 周期 1 | 周期 2 | 周期 3 | 2→3 变化 | open 单独 | closed 单独 |", "|---|---|---|---|---|---|---|"]
for k in keys:
    pc = S[k]["per_cycle"]; c = [p["Fx"] for p in pc]
    cyc_rows.append(f"| {LAB.get(k,k)} | " + " | ".join(f"{x:.1f}" for x in c) + f" | {100*(c[-1]/c[-2]-1):+.1f} % | {pc[-1]['Fx_open']:.1f} | {pc[-1]['Fx_closed']:.1f} |")
md = "# 步态验证算例：深入对比（两态合成，最后一个完整周期）\n\n" + "\n".join(rows) + "\n\n逐周期平均推力 [mN]（收敛性）：\n" + "\n".join(cyc_rows) + "\n"
open(os.path.join(OUT, "verify_table.md"), "w").write(md); print(md)
json.dump({k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in r.items()} for k, r in S.items()}, open(os.path.join(OUT, "verify_analysis.json"), "w"), indent=1, ensure_ascii=False)

# ---------------- 图 1：时间历程（相位锁定） ----------------
fig, axs = plt.subplots(2, 2, figsize=(13, 7.6), facecolor=SURF)
panels = [("Fx", "单腿推力 F_x [mN]"), ("Fz", "桨板垂向力 F_z [mN]（无用分量，浮力已支撑重量）"), ("My", "髋轴力矩 M_y [mN·m]（∝ 舵机负荷）"), ("P", "水动力功率 P [mW]")]
for ax, (key, lab) in zip(axs.ravel(), panels):
    ax.set_facecolor(SURF)
    for k in keys:
        r = S[k]; ax.plot(r["tau"], r[key], color=COL.get(k, GRAY), lw=1.6, label=f"{LAB.get(k,k)}（DUTY {r['DUTY']:.2f}）")
        ax.axvline(r["DUTY"], color=COL.get(k, GRAY), lw=0.7, ls=":", alpha=0.7)
    ax.axhline(0, color=GRAY, lw=0.8); ax.set_xlabel("周期相位 τ（点线 = 各方案划水相结束）"); ax.set_ylabel(lab); ax.grid(color=GRID, lw=0.8)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
axs[0, 0].legend(frameon=False, fontsize=8.5, loc="upper right"); axs[0, 0].set_title("CFD 两态合成，最后一个周期，0.02 s 平滑；τ 锁相（各方案周期不同）", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_verify_time.png"), dpi=115, facecolor=SURF); plt.close()

# ---------------- 图 2：权衡图（与模型前沿） ----------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), facecolor=SURF)
ax = axs[0]; ax.set_facecolor(SURF)
if pareto_fn:
    PJ = json.load(open(pareto_fn)); fr = PJ["pareto_F0_M"]
    ax.plot([f["M"] for f in fr], [f["F0"] for f in fr], color=GRAY, lw=1.4, label="行程模型前沿（同力矩最大推力）")
for k in keys:
    r = S[k]; m = MODEL.get(k, {})
    ax.scatter([r["My_peak"]], [r["Fx_mean"]], s=110, color=COL.get(k, GRAY), zorder=5, marker="o", label=f"{LAB.get(k,k)} CFD")
    if m.get("M"): ax.scatter([m["M"]], [m["F0"]], s=60, facecolors="none", edgecolors=COL.get(k, GRAY), zorder=4, marker="o"); ax.plot([m["M"], r["My_peak"]], [m["F0"], r["Fx_mean"]], color=COL.get(k, GRAY), lw=0.8, ls="--")
    ax.annotate(LAB.get(k, k), (r["My_peak"], r["Fx_mean"]), xytext=(6, 4), textcoords="offset points", fontsize=8.5, color=COL.get(k, GRAY))
ax.plot([], [], "o", mfc="none", mec=GRAY, label="模型预测（空心）"); ax.set_xlabel("髋轴峰值力矩 |M_y| [mN·m]"); ax.set_ylabel("周期平均推力 F̄_x [mN]"); ax.set_title("推力 vs 舵机负荷：CFD 落在模型前沿附近", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8, loc="upper left"); ax.set_xlim(0, 30); ax.set_ylim(0, 70); ax.grid(color=GRID, lw=0.8)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax = axs[1]; ax.set_facecolor(SURF)
for k in keys:
    r = S[k]; ax.scatter([r["Fx_mean"]], [r["Fz_mean"]], s=110, color=COL.get(k, GRAY), zorder=5); ax.annotate(f"{LAB.get(k,k)}  F̄z/F̄x = {r['Fz_mean']/r['Fx_mean']:.2f}", (r["Fx_mean"], r["Fz_mean"]), xytext=(6, 4), textcoords="offset points", fontsize=8.5, color=COL.get(k, GRAY))
xx = np.linspace(0, 60, 5); ax.plot(xx, 0.6 * xx, color=GRAY, lw=0.9, ls="--"); ax.text(45, 29, "F̄_z = 0.6 F̄_x", fontsize=8.5, color=INK2)
ax.set_xlabel("周期平均推力 F̄_x [mN]"); ax.set_ylabel("周期平均垂向力 F̄_z [mN]（无用分量）"); ax.set_title("推力 vs 垂向力：V3 把 F̄_z/F̄_x 压到 0.3", loc="left", fontsize=10.5); ax.set_xlim(0, 60); ax.set_ylim(0, 32); ax.grid(color=GRID, lw=0.8)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax = axs[2]; ax.set_facecolor(SURF)
for k in keys:
    r = S[k]; ax.scatter([r["P_mean"]], [r["Fx_mean"]], s=110, color=COL.get(k, GRAY), zorder=5); ax.annotate(f"{LAB.get(k,k)}  {r['Fx_mean']/r['P_mean']:.2f} mN/mW", (r["P_mean"], r["Fx_mean"]), xytext=(6, 6 if k == "V1" else -12), textcoords="offset points", fontsize=8.5, color=COL.get(k, GRAY))
for e in (2.5, 3.0, 3.5): xx = np.linspace(0, 20, 5); ax.plot(xx, e * xx, color=GRID, lw=0.9); xe = min(19.5, 57 / e); ax.text(xe, e * xe + 0.5, f"{e} mN/mW", fontsize=7.5, color=GRAY, ha="right")
ax.set_xlabel("周期平均水动力功率 P̄ [mW]"); ax.set_ylabel("周期平均推力 F̄_x [mN]"); ax.set_title("推力 vs 功率：快步态效率略降", loc="left", fontsize=10.5); ax.set_xlim(0, 20); ax.set_ylim(0, 60); ax.grid(color=GRID, lw=0.8)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_verify_tradeoff.png"), dpi=115, facecolor=SURF); plt.close()

# ---------------- 图 3：逐周期收敛 + open/closed 分解 + 两对腿叠加垂向力 ----------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), facecolor=SURF)
ax = axs[0]; ax.set_facecolor(SURF)
for k in keys:
    pc = S[k]["per_cycle"]; ax.plot([p["k"] for p in pc], [p["Fx"] for p in pc], "o-", color=COL.get(k, GRAY), lw=1.6, label=LAB.get(k, k))
ax.set_xlabel("周期序号"); ax.set_ylabel("周期平均推力 F̄_x [mN]"); ax.set_title("逐周期收敛：第 2→3 周期变化 ≤ 3 %", loc="left", fontsize=10.5); ax.set_xticks([1, 2, 3]); ax.legend(frameon=False, fontsize=8.5); ax.grid(color=GRID, lw=0.8); ax.set_ylim(0, None)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax = axs[1]; ax.set_facecolor(SURF); x = np.arange(len(keys)); w = 0.26
ax.bar(x - w, [S[k]["per_cycle"][-1]["Fx_open"] for k in keys], w, color=[COL.get(k, GRAY) for k in keys], label="open 桨全周期")
ax.bar(x, [S[k]["per_cycle"][-1]["Fx_closed"] for k in keys], w, color=[COL.get(k, GRAY) for k in keys], alpha=0.45, hatch="//", edgecolor=SURF, label="closed 桨全周期")
ax.bar(x + w, [S[k]["Fx_mean"] for k in keys], w, color=[COL.get(k, GRAY) for k in keys], alpha=0.8, edgecolor=INK, linewidth=1.0, label="两态合成")
for i, k in enumerate(keys): ax.text(i + w, S[k]["Fx_mean"] + 1, f"{S[k]['Fx_mean']:.1f}", ha="center", fontsize=8.5)
ax.set_xticks(x); ax.set_xticklabels([LAB.get(k, k) for k in keys], fontsize=9); ax.set_ylabel("周期平均推力 [mN]"); ax.set_title("两态合成 = 划水相取 open + 回收相取 closed", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5); ax.grid(axis="y", color=GRID, lw=0.8)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax = axs[2]; ax.set_facecolor(SURF)
for k in keys:
    r = S[k]; N = len(r["tau"]); half = N // 2; pair = r["Fz"] + np.roll(r["Fz"], -half)
    ax.plot(r["tau"], pair, color=COL.get(k, GRAY), lw=1.6, label=f"{LAB.get(k,k)}  峰峰 {r['Fz_pair_pp']:.0f} mN")
ax.axhline(0, color=GRAY, lw=0.8); ax.set_xlabel("周期相位 τ"); ax.set_ylabel("两对腿（半周期错开）合成垂向力 [mN]（每对 1 条腿）"); ax.set_title("对角步态叠加后的垂向力残差 = 船体升沉激励（Cui 的约束量）", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5); ax.grid(color=GRID, lw=0.8)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_verify_cycles.png"), dpi=115, facecolor=SURF); plt.close()
print("→", OUT)
