#!/usr/bin/env python3
"""
在 CFD 标定的行程模型上扫描步态参数，画“推力 – 舵机峰值力矩 – 升力”的权衡前沿（Pareto），回答“方案 A 是不是最优”。
    python3 scan_gait_pareto.py [results_dir]
变量：幅度 A、中心角 ψc（现机构范围 70–150°）、速率倍数 rf（峰值髋角速度 / 现值）、回收时长 t_rec（慢蜷缩摊到整个回收相）、蜷缩深度 Δφ
输出：results/gait_pareto.json、fig_gait_pareto.png
"""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results_leg")
os.chdir(HERE); sys.argv = ["analyze_stroke.py", R]
src = open("analyze_stroke.py").read(); cut = src.index("# ---------------- 扫描")
g = {"__file__": os.path.abspath("analyze_stroke.py"), "__name__": "__main__"}
import matplotlib; matplotlib.use("Agg")
exec(compile(src[:cut], "analyze_stroke.py", "exec"), g)
cand = g["candidate"]; T = g["T"]; DUTY = g["DUTY"]; SWEEP = g["SWEEP"]; PSI0 = g["PSI0"]; DPHI0 = g["DPHI0"]; WMAX = g["WMAX"]
J_SW0, V_SW0, T_SW0, J_back, V_BACK0, P_rec0, J_rec0 = g["J_SW0"], g["V_SW0"], g["T_SW0"], g["J_back"], g["V_BACK0"], g["P_rec0"], g["J_rec0"]
import matplotlib.pyplot as plt
from matplotlib import font_manager
INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"

def with_recovery(c, t_rec, dphi, U):
    """强制回收相时长 t_rec（慢蜷缩摊满），重算损失与周期平均；t_rec 不能短于两次切换的最短时间"""
    t_min = 2 * T_SW0 * (dphi / DPHI0) + 0.05
    t_rec = max(t_rec, t_min)
    t_sw_eff = t_rec / 2; v_sw = V_SW0 * (dphi / DPHI0) * (T_SW0 / t_sw_eff)
    J_sw = J_SW0 * (dphi / DPHI0) ** 2 * (T_SW0 / t_sw_eff) * ((v_sw + U) / v_sw) ** 2 if v_sw > 0 else 0.0
    v_back = V_BACK0 * (c["A"] / SWEEP) * (DUTY * T / t_rec)
    J_bk = J_back * (c["A"] / SWEEP) * (v_back / V_BACK0) * ((v_back + U) / v_back) ** 2 if v_back > 0 else 0.0
    J_rec = J_sw + J_bk; T_ = c["t_p"] + t_rec
    P_rec = P_rec0 * T * abs(J_rec) / abs(J_rec0); P_pow = c["Pmean"] * c["T"] - P_rec0 * T * abs(c["J_rec"]) / abs(J_rec0)
    o = dict(c); o.update(dphi=dphi, t_rec=t_rec, T=T_, f=1 / T_, duty=c["t_p"] / T_, J_sw=J_sw, J_bk=J_bk, J_rec=J_rec,
                          Fmean=(c["J_pow"] + J_rec) / T_, Fz_mean=c["Fz_mean"] * c["T"] / T_, Pmean=(P_pow + P_rec) / T_)
    o["eff"] = o["Fmean"] / o["Pmean"] if o["Pmean"] > 0 else None; o.pop("t", None); o.pop("Fx", None); return o

PSI_MIN, PSI_MAX = 70.0, 150.0
As = [40, 45, 50, 55, 60, 70, 80]; rfs = [1.0, 1.25, 1.43, 1.75, 2.0]; trecs = [0.45, 0.6, 0.8, 1.0, 1.3]
dphis = [float(x) for x in os.environ.get("DPHIS", "65").split(",")]           # 默认只用 CFD 验证过的蜷缩深度 65°；DPHIS=65,45,30 可加浅蜷缩
rows = []
for A in As:
    for pc in range(80, 131, 5):
        if pc + A / 2 > PSI_MAX or pc - A / 2 < PSI_MIN: continue
        for rf in rfs:
            base = {U: cand(A, pc, DPHI0, U, "rate", "slowtuck", n_t=160, rf=rf) for U in (0.0, 0.05, 0.1)}
            for tr in trecs:
                for dp in dphis:
                    o0 = with_recovery(base[0.0], tr, dp, 0.0); o5 = with_recovery(base[0.05], tr, dp, 0.05); o10 = with_recovery(base[0.1], tr, dp, 0.1)
                    rows.append(dict(A=A, psi_c=pc, rf=rf, t_p=o0["t_p"], t_rec=o0["t_rec"], T=o0["T"], f=o0["f"], duty=o0["duty"], dphi=dp,
                                     F0=o0["Fmean"], F5=o5["Fmean"], F10=o10["Fmean"], Fz0=o0["Fz_mean"], P0=o0["Pmean"], eff0=o0["eff"], M=o0["M_peak"], vtip=o0["v_tip_max"],
                                     inrange=(pc + A / 2 <= PSI0 + 1e-6 and pc - A / 2 >= PSI0 - SWEEP - 1e-6)))
print(len(rows), "candidates")
import pandas as pd
df = pd.DataFrame(rows)
cur = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "tuck"); cur5 = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.05, "rate", "tuck")
A_ = with_recovery(cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", rf=0.625 / 0.4375), 0.8125, DPHI0, 0.0)
A5 = with_recovery(cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.05, "rate", "slowtuck", rf=0.625 / 0.4375), 0.8125, DPHI0, 0.05)
B_ = with_recovery(cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", rf=1.0), 0.625 / 0.35 * 0.65, DPHI0, 0.0)

# ---- 力矩预算下的最优（U=0 与 U=0.05 两种目标），只在现髋角范围内（间隙已核对）与全机构范围各取一个
def best_under(df_, Mlim, key, inrange_only):
    d = df_[(df_.M <= Mlim)]
    if inrange_only: d = d[d.inrange]
    if len(d) == 0: return None
    return d.loc[d[key].idxmax()].to_dict()
M0 = cur["M_peak"]; out = dict(current=dict(F0=cur["Fmean"], F5=cur5["Fmean"], M=cur["M_peak"], Fz0=cur["Fz_mean"], P0=cur["Pmean"]),
                               schemeA=dict(F0=A_["Fmean"], F5=A5["Fmean"], M=A_["M_peak"], Fz0=A_["Fz_mean"], P0=A_["Pmean"], duty=A_["duty"], T=A_["T"]),
                               schemeB=dict(F0=B_["Fmean"], M=B_["M_peak"], Fz0=B_["Fz_mean"], P0=B_["Pmean"], duty=B_["duty"], T=B_["T"]), M0=M0, best={})
for k in (1.0, 1.5, 2.0, 3.0, 4.0):
    for key in ("F0", "F5"):
        for inr in (True, False):
            b = best_under(df, k * M0 + 1e-9, key, inr)
            if b: out["best"][f"M{k}x_{key}_{'inrange' if inr else 'full'}"] = b
            if b: print(f"M ≤ {k}×  目标 {key} {'现范围' if inr else '全范围'}: A={b['A']} ψc={b['psi_c']} rf={b['rf']} t_p={b['t_p']:.2f} t_rec={b['t_rec']:.2f} duty={b['duty']:.2f} Δφ={b['dphi']}  F0={b['F0']:.1f} F5={b['F5']:.1f} F10={b['F10']:.1f} Fz={b['Fz0']:.1f} M={b['M']:.1f} P={b['P0']:.1f} eff={b['eff0']:.2f}")
print(f"current F0={cur['Fmean']:.1f} F5={cur5['Fmean']:.1f} M={M0:.1f};  A: F0={A_['Fmean']:.1f} F5={A5['Fmean']:.1f} M={A_['M_peak']:.1f} Fz={A_['Fz_mean']:.1f};  B: F0={B_['Fmean']:.1f} M={B_['M_peak']:.1f}")
# Pareto 前沿（F0 最大 vs M 最小）
d = df.sort_values("M"); front = []; fmax = -1e9
for _, r in d.iterrows():
    if r.F0 > fmax: front.append(r.to_dict()); fmax = r.F0
out["pareto_F0_M"] = front
d = df.sort_values("M"); front5 = []; fmax = -1e9
for _, r in d.iterrows():
    if r.F5 > fmax: front5.append(r.to_dict()); fmax = r.F5
out["pareto_F5_M"] = front5
def conv(o):
    if isinstance(o, dict): return {k: conv(v) for k, v in o.items()}
    if isinstance(o, list): return [conv(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    if isinstance(o, (np.bool_,)): return bool(o)
    return o
TAG = "" if dphis == [65.0] else "_dphi" + "-".join(f"{int(x)}" for x in dphis)
json.dump(conv(out), open(os.path.join(R, f"gait_pareto{TAG}.json"), "w"), indent=1, ensure_ascii=False)
df.to_csv(os.path.join(R, f"gait_scan{TAG}.csv"), index=False)

# ---- 图
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), facecolor=SURF)
sc_kw = dict(s=9, alpha=0.35, linewidths=0)
for ax, key, ttl in ((axs[0], "F0", "U = 0：周期平均推力 vs 舵机峰值力矩"), (axs[1], "F5", "U = 0.05 m/s：推力 vs 舵机峰值力矩")):
    ax.set_facecolor(SURF)
    sca = ax.scatter(df.M, df[key], c=df.duty, cmap="viridis", vmin=0.2, vmax=0.6, **sc_kw)
    fr = out["pareto_F0_M"] if key == "F0" else out["pareto_F5_M"]
    ax.plot([f["M"] for f in fr], [f[key] for f in fr], color=INK, lw=1.6, label="前沿（同力矩下的最大推力）")
    ax.scatter([cur["M_peak"]], [cur["Fmean"] if key == "F0" else cur5["Fmean"]], marker="*", s=180, color=RED, zorder=5, label="现步态（模型）")
    ax.scatter([A_["M_peak"]], [A_["Fmean"] if key == "F0" else A5["Fmean"]], marker="D", s=70, color=ORANGE, zorder=5, label="方案 A（狗式占空比 0.35）")
    if key == "F0": ax.scatter([B_["M_peak"]], [B_["Fmean"]], marker="s", s=60, color=BLUE, zorder=5, label="方案 B（回收拉长）")
    for k in (1.0, 2.0, 3.0): ax.axvline(k * M0, color=GRAY, lw=0.8, ls=":"); ax.text(k * M0, ax.get_ylim()[0] if False else (df[key].max() * 1.02), f"{k:g}× 现力矩", fontsize=8, color=INK2, ha="center")
    ax.set_xlabel("髋轴峰值力矩 |M| [mN·m]（∝ 舵机负荷）"); ax.set_ylabel("周期平均推力 F̄_x [mN]"); ax.set_title(ttl, loc="left", fontsize=10.5)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(color=GRID, lw=0.8); ax.set_xlim(0, 45)
axs[0].legend(frameon=False, fontsize=8.5, loc="lower right")
cax = axs[1].inset_axes([0.60, 0.10, 0.34, 0.035]); cb = fig.colorbar(sca, cax=cax, orientation="horizontal"); cb.set_label("划水相占空比 t_p / T", fontsize=8.5); cb.ax.tick_params(labelsize=8)
ax = axs[2]; ax.set_facecolor(SURF)
sc3 = ax.scatter(df.F0, df.Fz0, c=df.rf, cmap="plasma", vmin=1.0, vmax=2.0, **sc_kw)
cax3 = ax.inset_axes([0.60, 0.10, 0.34, 0.035]); cb3 = fig.colorbar(sc3, cax=cax3, orientation="horizontal"); cb3.set_label("速率倍数 rf（峰值髋角速度 / 现值）", fontsize=8.5); cb3.ax.tick_params(labelsize=8)
ax.scatter([cur["Fmean"]], [cur["Fz_mean"]], marker="*", s=180, color=RED, zorder=5); ax.scatter([A_["Fmean"]], [A_["Fz_mean"]], marker="D", s=70, color=ORANGE, zorder=5)
xx = np.linspace(0, df.F0.max(), 10); ax.plot(xx, 0.6 * xx, color=GRAY, lw=0.9, ls="--"); ax.text(df.F0.max() * 0.45, 0.6 * df.F0.max() * 0.45 + 4, "F̄_z = 0.6 F̄_x（现步态的升力比）", fontsize=8, color=INK2)
ax.set_xlabel("周期平均推力 F̄_x [mN]（U = 0）"); ax.set_ylabel("周期平均升力 F̄_z [mN]"); ax.set_title(f"推力越大升力越大（颜色 = 速率倍数；Δφ = {'/'.join(f'{int(x)}' for x in dphis)}°）", loc="left", fontsize=10.5)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.grid(color=GRID, lw=0.8)
plt.tight_layout(); plt.savefig(os.path.join(R, f"fig_gait_pareto{TAG}.png"), dpi=115, facecolor=SURF); print(f"→ fig_gait_pareto{TAG}.png")
