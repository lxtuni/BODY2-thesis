#!/usr/bin/env python3
"""
步态验证算例的收集与对比（V1 前沿点 / V2 狗式占空比 / V3 低中心角 … vs 现步态基准 vs 行程模型预测）
    python3 compare_verify.py collect <算例目录> <方案名> <fan>           # 把 postProcessing 的 force/moment 与 gait.json 拷到 results_verify/<方案名>/
    python3 compare_verify.py analyze [--root results_verify] [--baseline results] [--out results_verify]
analyze：每个方案目录需有 force_open.dat / force_closed.dat / moment_*.dat / gait_*.json；两态按 τ < DUTY 取 open、其余取 closed 合成，
        取最后一个完整周期，算周期平均 F_x / F_z / M_y、峰值、各相位冲量、功率，与模型预测（MODEL 表）并排；输出 verify_summary.json/.md 与 fig_verify.png。
"""
import os, sys, glob, json, argparse, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
# 行程模型预测（scan_gait_pareto.py / eval_dog_schemes.py，U = 0；current 的 CFD 实测 23.4 mN）
MODEL = {"current": dict(F0=21.0, F5=2.9, Fz=12.1, M=9.0, P=7.3, label="现步态（T 1.25, 60°, 136→76°, DUTY 0.5, 快蜷缩）"),
         "V1": dict(F0=47.7, F5=20.7, Fz=23.3, M=18.5, P=17.3, label="V1 前沿点（T 0.89, 60°, 135→75°, DUTY 0.49, 慢蜷缩 0.45 s）"),
         "V2": dict(F0=36.5, F5=17.0, Fz=16.7, M=18.5, P=11.6, label="V2 狗式占空比（T 1.25, 60°, 136→76°, DUTY 0.35, 慢蜷缩 0.81 s）"),
         "V3": dict(F0=42.6, F5=17.5, Fz=13.0, M=17.9, P=15.4, label="V3 低中心角（T 0.81, 50°, 120→70°, DUTY 0.45, 慢蜷缩 0.45 s）")}

def cat_dat(case, name):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", "forces", "*", name)))
    if not files: return None, None
    rows = []
    for fn in files:
        d = np.loadtxt(fn, comments="#")
        if d.ndim == 1: d = d[None]
        rows.append(d)
    d = np.concatenate(rows); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True)
    return d[i], open(files[0]).read().split("\n")[:4]

def collect(case, scheme, fan, root):
    out = os.path.join(root, scheme); os.makedirs(out, exist_ok=True); case = os.path.expanduser(case)
    for name, tag in (("force.dat", "force"), ("moment.dat", "moment")):
        d, hdr = cat_dat(case, name)
        if d is None: print(f"!! {case}: 没有 postProcessing/forces/*/{name}"); continue
        with open(os.path.join(out, f"{tag}_{fan}.dat"), "w") as fh:
            fh.write("\n".join(h for h in hdr if h.startswith("#")) + "\n"); np.savetxt(fh, d, fmt="%.9e")
        print(f"{scheme}/{tag}_{fan}.dat  ← {len(d)} 行，t = {d[0,0]:.4f} … {d[-1,0]:.4f} s")
    gj = os.path.join(case, "gait.json")
    if os.path.exists(gj):
        import shutil; shutil.copyfile(gj, os.path.join(out, f"gait_{fan}.json"))
    lg = os.path.join(case, "log.overPimpleDyMFoam")
    if os.path.exists(lg):
        txt = open(lg, errors="ignore").read(); import re
        m = re.findall(r"ClockTime = ([\d.]+) s", txt); n = len(re.findall(r"^Time = ", txt, re.M))
        json.dump(dict(steps=n, clock_s=float(m[-1]) if m else None, end="End" in txt[-2000:]), open(os.path.join(out, f"run_{fan}.json"), "w"))
        print(f"   {n} 步，用时 {float(m[-1])/3600:.2f} h" if m else "")

def analyze_scheme(R):
    def load(tag):
        fn = os.path.join(R, f"force_{tag}.dat")
        if not os.path.exists(fn): return None
        d = np.loadtxt(fn, comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
        mfn = os.path.join(R, f"moment_{tag}.dat")
        m = np.loadtxt(mfn, comments="#") if os.path.exists(mfn) else None
        if m is not None: m = m[np.argsort(m[:, 0])]; _, i = np.unique(m[:, 0], return_index=True); m = m[i]
        g = json.load(open(os.path.join(R, f"gait_{tag}.json")))
        return dict(t=d[:, 0], F=d[:, 1:4], P=d[:, 4:7], V=d[:, 7:10], M=(m[:, 1:4] if m is not None else np.zeros_like(d[:, 1:4])), tm=(m[:, 0] if m is not None else d[:, 0]), gait=g)
    cases = {k: v for k, v in ((k, load(k)) for k in ("open", "closed", "plate")) if v is not None}   # 无海象运算符：兼容 Ubuntu 18.04 的 python 3.6
    if not cases: return None
    g = (cases.get("open") or cases.get("closed") or cases.get("plate"))["gait"]; T = g["T"]; DUTY = g["DUTY"]; FC = g.get("FAN_CLOSE", -1); tc = DUTY if (FC is None or FC < 0) else FC
    ncyc = int(np.floor(min(c["t"][-1] for c in cases.values()) / T + 1e-6))
    if ncyc < 1: return dict(partial=True, t_end=float(min(c["t"][-1] for c in cases.values())), T=T)
    t0, t1 = (ncyc - 1) * T, ncyc * T; N = 2500; tt = np.linspace(t0, t1, N, endpoint=False); tau = (tt - t0) / T; dt = T / N
    def smooth(y, w=0.02):
        k = max(1, int(round(w / dt))); ker = np.ones(k) / k; yp = np.concatenate([y[-k:], y, y[:k]]); return np.convolve(yp, ker, "same")[k:-k]
    def cols(c, key, tkey="t"): return np.stack([np.interp(tt, c[tkey], c[key][:, i]) for i in range(3)], 1)
    res = {tag: dict(F=cols(c, "F"), P=cols(c, "P"), V=cols(c, "V"), M=cols(c, "M", "tm")) for tag, c in cases.items()}
    if "plate" in res:                                                   # 第二套：刚性板 + 脚踝羽化，一次算例含整个周期，不合成
        comp = res["plate"]; mode = "plate_single"
    elif "open" in res and "closed" in res:
        fo = (tau < tc)[:, None]; comp = {k: np.where(fo, res["open"][k], res["closed"][k]) for k in res["open"]}; mode = "composite"
    else:
        comp = res.get("open") or res.get("closed"); mode = "open_only" if "open" in res else "closed_only"
    sf = g.get("SWITCH_FRAC", 0.3); tr = min(max(sf, 0.05), 0.5) * (1 - DUTY)
    phases = [("划水相", 0.0, DUTY), ("收拢切换", DUTY, DUTY + tr), ("蜷缩前扫", DUTY + tr, 1 - tr), ("伸展切换", 1 - tr, 1.0)]
    def contrib(y): return [float(np.trapezoid(y[(tau >= a_) & (tau < b_)], tt[(tau >= a_) & (tau < b_)]) / T * 1e3) if ((tau >= a_) & (tau < b_)).any() else 0.0 for _, a_, b_ in phases]
    Fx, Fz, My = comp["F"][:, 0], comp["F"][:, 2], comp["M"][:, 1]
    # 功率：P = −(F·v_E + M_E·ω)，运动学来自 leg_dynamics（与 analyze_leg_xyz.py 相同）
    power = None
    try:
        for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
            if k in g: os.environ[k] = str(g[k])
        if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
        sp = importlib.util.spec_from_file_location("ld_" + os.path.basename(R), os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(sp); sp.loader.exec_module(ld)
        psi, phi, _ = ld.gait(tt); PP = [ld.pose(p, q) for p, q in zip(psi, phi)]; E = np.array([p["E"] for p in PP])
        vE = np.gradient(E, tt, axis=0); om = -np.gradient(np.unwrap(psi), tt); O2 = np.array(ld.O2)
        rE = E - O2; ME_y = My - (rE[:, 1] * comp["F"][:, 0] - rE[:, 0] * comp["F"][:, 2])      # 力矩从 O2 搬到 E：M_E = M_O2 − r × F（绕 y）
        power = -(comp["F"][:, 0] * vE[:, 0] + comp["F"][:, 2] * vE[:, 1] + ME_y * om)
    except Exception as e:
        print("  功率计算跳过：", e)
    out = dict(mode=mode, ncyc=ncyc, T=T, DUTY=DUTY, f=1 / T, gait=dict(SWEEP=g.get("SWEEP"), psi0=g.get("psi0_deg"), SWITCH_FRAC=sf, PHI_RET=g.get("PHI_RET"), U=g.get("U", 0.0)),
               Fx_mean=float(Fx.mean() * 1e3), Fz_mean=float(Fz.mean() * 1e3), My_mean=float(My.mean() * 1e3),
               Fx_peak_s=float(smooth(Fx).max() * 1e3), Fx_min_s=float(smooth(Fx).min() * 1e3), Fz_peak_s=float(np.abs(smooth(Fz)).max() * 1e3), My_peak_s=float(np.abs(smooth(My)).max() * 1e3),
               Fx_phases=contrib(Fx), Fz_phases=contrib(Fz), phases=[p[0] for p in phases], phase_ranges=[(p[1], p[2]) for p in phases],
               P_mean_mW=float(power.mean() * 1e3) if power is not None else None, P_peak_mW=float(smooth(power).max() * 1e3) if power is not None else None,
               tau=tau[::10].tolist(), Fx_s=(smooth(Fx) * 1e3)[::10].tolist(), Fz_s=(smooth(Fz) * 1e3)[::10].tolist(),
               open_mean=float(res["open"]["F"][:, 0].mean() * 1e3) if "open" in res else None, closed_mean=float(res["closed"]["F"][:, 0].mean() * 1e3) if "closed" in res else None)
    for tag in ("open", "closed", "plate"):
        rj = os.path.join(R, f"run_{tag}.json")
        if os.path.exists(rj): out[f"run_{tag}"] = json.load(open(rj))
    return out

def analyze(root, baseline, outdir):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc"):
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
        except Exception: pass
    plt.rcParams["axes.unicode_minus"] = False
    INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"
    COL = {"current": INK, "V1": GREEN, "V2": ORANGE, "V3": BLUE}
    schemes = {}
    if baseline and os.path.isdir(baseline) and os.path.exists(os.path.join(baseline, "force_open.dat")): schemes["current"] = baseline
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if os.path.isdir(d) and glob.glob(os.path.join(d, "force_*.dat")): schemes[os.path.basename(d)] = d
    if not schemes: sys.exit("没有任何方案结果")
    S = {}
    for k, d in schemes.items():
        r = analyze_scheme(d)
        if r: S[k] = r; print(f"[{k}] " + (f"未满一个周期（t={r['t_end']:.2f}/{r['T']} s）" if r.get("partial") else f"{r['mode']}  第 {r['ncyc']} 周期  F̄x {r['Fx_mean']:.1f} mN  F̄z {r['Fz_mean']:.1f}  |My|峰 {r['My_peak_s']:.1f} mN·m  P̄ {r['P_mean_mW'] if r['P_mean_mW'] is None else round(r['P_mean_mW'],1)} mW  相位 {np.round(r['Fx_phases'],1).tolist()}"))
    os.makedirs(outdir, exist_ok=True)
    json.dump(dict(schemes=S, model=MODEL), open(os.path.join(outdir, "verify_summary.json"), "w"), indent=1, ensure_ascii=False)
    # ---- 表
    L = ["# 步态验证算例 vs 行程模型\n", "| 方案 | 步态 | 状态 | CFD F̄_x [mN] | 模型 F̄_x | 偏差 | CFD F̄_z | 模型 F̄_z | |M_y| 峰 [mN·m] | 模型 |M| | P̄ [mW] | 划水相 / 回收相贡献 [mN] | F̄_x 相对现步态 |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    ref = S.get("current", {}).get("Fx_mean")
    for k, r in S.items():
        m = MODEL.get(k, {}); lab = m.get("label", k)
        if r.get("partial"): L.append(f"| {k} | {lab} | 未满 1 周期 | — | {m.get('F0','—')} | — | — | {m.get('Fz','—')} | — | {m.get('M','—')} | — | — | — |"); continue
        dev = f"{100*(r['Fx_mean']/m['F0']-1):+.0f} %" if m.get("F0") else "—"
        rel = f"{100*(r['Fx_mean']/ref-1):+.0f} %" if (ref and k != "current") else "—"
        L.append(f"| {k} | {lab} | {r['mode']}，第 {r['ncyc']} 周期 | **{r['Fx_mean']:.1f}** | {m.get('F0','—')} | {dev} | {r['Fz_mean']:.1f} | {m.get('Fz','—')} | {r['My_peak_s']:.1f} | {m.get('M','—')} | {r['P_mean_mW'] if r['P_mean_mW'] is None else round(r['P_mean_mW'],1)} | {r['Fx_phases'][0]:+.1f} / {sum(r['Fx_phases'][1:]):+.1f} | {rel} |")
    open(os.path.join(outdir, "verify_summary.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    # ---- 图：左 F̄x/F̄z CFD vs 模型；右 合成 F_x(τ)
    full = {k: r for k, r in S.items() if not r.get("partial")}
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.4), facecolor=SURF)
    ax = axs[0]; ax.set_facecolor(SURF); keys = list(full.keys()); x = np.arange(len(keys)); w = 0.2
    ax.bar(x - 1.5 * w, [full[k]["Fx_mean"] for k in keys], w, color=[COL.get(k, GRAY) for k in keys], label="CFD F̄_x")
    ax.bar(x - 0.5 * w, [MODEL.get(k, {}).get("F0", 0) for k in keys], w, color=[COL.get(k, GRAY) for k in keys], alpha=0.4, hatch="//", edgecolor=SURF, label="模型 F̄_x")
    ax.bar(x + 0.5 * w, [full[k]["Fz_mean"] for k in keys], w, color=RED, label="CFD F̄_z")
    ax.bar(x + 1.5 * w, [MODEL.get(k, {}).get("Fz", 0) for k in keys], w, color=RED, alpha=0.4, hatch="//", edgecolor=SURF, label="模型 F̄_z")
    for i, k in enumerate(keys): ax.text(i - 1.5 * w, full[k]["Fx_mean"] + 0.8, f"{full[k]['Fx_mean']:.1f}", ha="center", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(keys); ax.set_ylabel("周期平均力 [mN]"); ax.set_title("周期平均推力 / 垂向力：CFD（实色）vs 行程模型（斜纹）", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5, ncol=2)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax = axs[1]; ax.set_facecolor(SURF)
    for k, r in full.items(): ax.plot(r["tau"], r["Fx_s"], color=COL.get(k, GRAY), lw=1.6, label=f"{k}  F̄x {r['Fx_mean']:.1f} mN（T {r['T']} s）")
    ax.axhline(0, color=GRAY, lw=0.8); ax.set_xlabel("周期相位 τ（各方案 DUTY 不同：划水相结束在 τ = DUTY）"); ax.set_ylabel("单腿 F_x [mN]（0.02 s 平滑）")
    ax.set_title("两态合成推力随相位（最后一个完整周期）", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(color=GRID, lw=0.8)
    plt.tight_layout(); plt.savefig(os.path.join(outdir, "fig_verify.png"), dpi=115, facecolor=SURF); print("→", os.path.join(outdir, "fig_verify.png"))

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "collect":
        ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("case"); ap.add_argument("scheme"); ap.add_argument("fan", choices=["open", "closed", "plate"]); ap.add_argument("--root", default=os.path.join(HERE, "results_verify"))
        a = ap.parse_args(); collect(a.case, a.scheme, a.fan, a.root)
    else:
        ap = argparse.ArgumentParser(); ap.add_argument("cmd", nargs="?", default="analyze"); ap.add_argument("--root", default=os.path.join(HERE, "results_verify"))
        ap.add_argument("--baseline", default=os.path.join(HERE, "results")); ap.add_argument("--out", default=None)
        a = ap.parse_args(); analyze(a.root, a.baseline, a.out or a.root)
