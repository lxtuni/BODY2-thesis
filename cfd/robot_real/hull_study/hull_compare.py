#!/usr/bin/env python3
"""
船体形状比选汇总：读 ~/run/hull_<name>_U<xx>/postProcessing/forces 的 force.dat，尾段平均，出表、图、JSON
    python3 hull_compare.py --run ~/run --out results --hydro hull_hydrostatics.json
输出：results/hull_resistance.json、hull_resistance.md、fig_hull_resistance.png、fig_hull_convergence.png，并拷贝各算例的 force.dat
"""
import os, sys, glob, json, argparse, shutil
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

def _cjk():
    import glob as g
    c = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]
    c += g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + g.glob("/mnt/c/Windows/Fonts/simhei.ttf")
    for f in c:
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return
        except Exception: pass
_cjk(); plt.rcParams["axes.unicode_minus"] = False
INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"

ap = argparse.ArgumentParser()
ap.add_argument("--run", default=os.path.expanduser("~/run")); ap.add_argument("--out", default="results")
ap.add_argument("--hydro", default="hull_hydrostatics.json"); ap.add_argument("--tail", type=float, default=0.2)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
HULLS = [("M1_current", "M1 现船体（平底）"), ("M2_roundbottom", "M2 圆滑底"), ("M3_streamlined", "M3 流线型圆舭")]
COL = {"M1_current": INK, "M2_roundbottom": BLUE, "M3_streamlined": GREEN}
RHO = 998.8
hyd = json.load(open(a.hydro)) if os.path.exists(a.hydro) else {}

def read_forces(case):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", "**", "force.dat"), recursive=True))
    if not files: return None
    rows = []
    for fn in files:
        d = np.loadtxt(fn, comments="#")
        if d.ndim == 1: d = d[None]
        rows.append(d)
    d = np.concatenate(rows); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); return d[i]

res = {}
cases = sorted(glob.glob(os.path.join(a.run, "hull_*_U??")))
for case in cases:
    name = os.path.basename(case); hull = name[5:name.rfind("_U")]; U = int(name[-2:]) / 10.0
    d = read_forces(case)
    if d is None: print("跳过（无 force.dat）:", case); continue
    n = len(d); k = int(round(n * (1 - a.tail))); tail = d[k:]
    Fx = -tail[:, 1].mean(); Px = -tail[:, 4].mean(); Vx = -tail[:, 7].mean()            # 阻力取 +（来流 −x，阻力沿 −x）
    Fz = tail[:, 3].mean(); Fy = tail[:, 2].mean()
    spread = float((tail[:, 1].max() - tail[:, 1].min()) / max(abs(tail[:, 1].mean()), 1e-12))
    h1, h2 = tail[: len(tail) // 2, 1].mean(), tail[len(tail) // 2:, 1].mean(); drift = float((h2 - h1) / max(abs(tail[:, 1].mean()), 1e-12))
    hh = hyd.get(hull, {})
    S = hh.get("S_cm2"); Af = hh.get("Amax_cm2")
    q = 0.5 * RHO * U ** 2
    res.setdefault(hull, {})[f"{U:.1f}"] = dict(case=case, U=U, steps=int(d[-1, 0]), Fx_mN=Fx * 1e3, Px_mN=Px * 1e3, Vx_mN=Vx * 1e3, Fz_N=Fz, Fy_mN=Fy * 1e3,
                                                tail_peak2peak=spread, tail_drift=drift, disp_g_cfd=Fz / 9.81 * 1e3,
                                                CD_front=(Fx / (q * Af * 1e-4)) if Af else None, CT_wet=(Fx / (q * S * 1e-4)) if S else None, CF_wet=(Vx / (q * S * 1e-4)) if S else None,
                                                P_mW=Fx * U * 1e3, hist=d[:, [0, 1, 4, 7]].tolist())
    shutil.copyfile(sorted(glob.glob(os.path.join(case, "postProcessing", "**", "force.dat"), recursive=True))[-1], os.path.join(a.out, f"force_{name}.dat"))
    print(f"[{name:<26}] 阻力 {Fx*1e3:6.1f} mN（压差 {Px*1e3:5.1f} / 摩擦 {Vx*1e3:4.1f}）  浮力 {Fz:.2f} N ≈ {Fz/9.81*1e3:.0f} g  尾段峰峰 {100*spread:.1f}% 漂移 {100*drift:+.1f}%  步数 {int(d[-1,0])}")
if not res: sys.exit("没有找到任何算例结果")
speeds = sorted({U for h in res.values() for U in h})
json.dump(dict(results={h: {U: {k: v for k, v in r.items() if k != "hist"} for U, r in hh.items()} for h, hh in res.items()}, hydro=hyd, tail_frac=a.tail), open(os.path.join(a.out, "hull_resistance.json"), "w"), indent=1, ensure_ascii=False)

# ---- 表（markdown）
lines = ["# 船体形状比选：裸船体阻力（interFoam LTS，水线 z = −10 mm，尾段 %d %% 平均）\n" % int(100 * a.tail)]
for U in speeds:
    lines.append(f"\n## U = {U} m/s\n")
    lines.append("| 船体 | 总阻力 [mN] | 压差 [mN] | 摩擦 [mN] | 相对 M1 | C_D(迎流面) | C_T(湿表面) | 排水量 CFD [g] | 有效功率 [mW] | 尾段峰峰 / 漂移 |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    ref = res.get("M1_current", {}).get(U, {}).get("Fx_mN")
    for h, lab in HULLS:
        r = res.get(h, {}).get(U)
        if not r: continue
        rel = f"{100 * (r['Fx_mN'] / ref - 1):+.1f} %" if ref else "—"
        lines.append(f"| {lab} | **{r['Fx_mN']:.1f}** | {r['Px_mN']:.1f} | {r['Vx_mN']:.1f} | {rel} | {r['CD_front']:.3f} | {r['CT_wet']:.4f} | {r['disp_g_cfd']:.0f} | {r['P_mW']:.1f} | {100*r['tail_peak2peak']:.1f} % / {100*r['tail_drift']:+.1f} % |" if r["CD_front"] else
                     f"| {lab} | **{r['Fx_mN']:.1f}** | {r['Px_mN']:.1f} | {r['Vx_mN']:.1f} | {rel} | — | — | {r['disp_g_cfd']:.0f} | {r['P_mW']:.1f} | {100*r['tail_peak2peak']:.1f} % / {100*r['tail_drift']:+.1f} % |")
if hyd:
    lines.append("\n## 静水力学（同水线）\n")
    lines.append("| 船体 | 排水量 [g] | 湿表面 [cm²] | 水线 L×B [mm] | 吃水 [mm] | 最大剖面 [cm²] | C_B | C_P | BM [mm] |"); lines.append("|---|---|---|---|---|---|---|---|---|")
    for h, lab in HULLS:
        q = hyd.get(h)
        if q: lines.append(f"| {lab} | {q['disp_g']:.0f} | {q['S_cm2']:.0f} | {q['LWL_mm']:.0f} × {q['BWL_mm']:.0f} | {q['T_mm']:.0f} | {q['Amax_cm2']:.1f} | {q['Cb']:.2f} | {q['Cp']:.2f} | {q['BM_mm']:.0f} |")
open(os.path.join(a.out, "hull_resistance.md"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))

# ---- 图 1：阻力分解柱状图
fig, axs = plt.subplots(1, len(speeds), figsize=(6.5 * len(speeds), 5), facecolor=SURF, squeeze=False)
for ax, U in zip(axs[0], speeds):
    ax.set_facecolor(SURF); xs = np.arange(len(HULLS)); w = 0.6
    for i, (h, lab) in enumerate(HULLS):
        r = res.get(h, {}).get(U)
        if not r: continue
        ax.bar(i, r["Px_mN"], w, color=COL[h], label="压差 (Druck)" if i == 0 else None)
        ax.bar(i, r["Vx_mN"], w, bottom=r["Px_mN"], color=COL[h], alpha=0.45, hatch="//", edgecolor=SURF, label="摩擦 (Reibung)" if i == 0 else None)
    ymax = max(res[h][U]["Fx_mN"] for h, _ in HULLS if h in res and U in res[h]); ax.set_ylim(0, 1.45 * ymax)
    for i, (h, lab) in enumerate(HULLS):
        r = res.get(h, {}).get(U)
        if r: ax.text(i, r["Fx_mN"] + 0.02 * ymax, f"{r['Fx_mN']:.1f} mN\n压差 {r['Px_mN']:.1f} / 摩擦 {r['Vx_mN']:.1f}", ha="center", va="bottom", fontsize=9, color=INK)
    ref = res.get("M1_current", {}).get(U, {}).get("Fx_mN")
    ax.set_xticks(xs); ax.set_xticklabels([lab + (f"\n{100*(res[h][U]['Fx_mN']/ref-1):+.0f} %" if (ref and h in res and U in res[h] and h != "M1_current") else "") for h, lab in HULLS], fontsize=9.5)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_ylabel("阻力 [mN]", color=INK); ax.set_title(f"裸船体阻力  U = {U} m/s（尾段 {int(100*a.tail)} % 平均）", loc="left", fontsize=11); ax.legend(frameon=False, fontsize=9, loc="upper right")
plt.tight_layout(); plt.savefig(os.path.join(a.out, "fig_hull_resistance.png"), dpi=110, facecolor=SURF); plt.close()

# ---- 图 2：收敛历程
fig, axs = plt.subplots(1, len(speeds), figsize=(6.5 * len(speeds), 4.6), facecolor=SURF, squeeze=False)
for ax, U in zip(axs[0], speeds):
    ax.set_facecolor(SURF)
    for h, lab in HULLS:
        r = res.get(h, {}).get(U)
        if not r: continue
        hst = np.array(r["hist"]); ax.plot(hst[:, 0], -hst[:, 1] * 1e3, color=COL[h], lw=1.6, label=f"{lab}  {r['Fx_mN']:.1f} mN")
        ax.axvspan(hst[-1, 0] * (1 - a.tail), hst[-1, 0], color=GRAY, alpha=0.08)
    ymax = max(res[h][U]["Fx_mN"] for h, _ in HULLS if h in res and U in res[h]) * 2.0
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.set_ylim(0, ymax); ax.grid(color=GRID, lw=0.8); ax.set_xlabel("LTS 伪时间步", color=INK); ax.set_ylabel("总阻力 [mN]", color=INK); ax.set_title(f"收敛历程  U = {U} m/s（灰色 = 平均窗口）", loc="left", fontsize=11); ax.legend(frameon=False, fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(a.out, "fig_hull_convergence.png"), dpi=110, facecolor=SURF); plt.close()
print("→", a.out)
