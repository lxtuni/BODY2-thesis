#!/usr/bin/env python3
"""
整机 overset 算例（make_full_cfd.py）的受力后处理：
    python3 full_forces.py --open ~/run/full_open_U01 --closed ~/run/full_closed_U01
    python3 full_forces.py --closed ~/run/full_closed_U01            # 只有一个也行
    python3 full_forces.py --sweep results_full/full_U*.json         # 多个航速 → 净推力–航速曲线、自航点
每条腿按自己的相位拼接两态（tau_leg < FAN_CLOSE 用 open，其余 closed），船体阻力两算例都报；
输出：results_full/full_forces_<tag>.png、full_U<U>.json，并拷 force_*.dat。
"""
import argparse, glob, json, os, sys
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
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
LEGCOL = {"FL": "#2a78d6", "FR": "#eb6834", "RL": "#3a9c6e", "RR": "#8e6bd6", "hull": "#a8a7a1", "total": "#0b0b0b"}

ap = argparse.ArgumentParser()
ap.add_argument("--open"); ap.add_argument("--closed")
ap.add_argument("--sweep", nargs="*", help="多个 full_U*.json → 推力–航速曲线")
ap.add_argument("--out", default=os.environ.get("COPY_TO", "/mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results_full"))
ap.add_argument("--no-linkage", action="store_true", help="不把连杆（曲柄/链接杆/摇臂/曲柄2/从动杆）水下阻力（条带法）加进腿的力里")
ap.add_argument("--wl", type=float, default=None, help="连杆湿长用的水线 z [m]（默认取算例 gait.json 的 WL）")
a = ap.parse_args()
OUT = os.path.abspath(os.path.expanduser(a.out))
if not os.path.isdir(os.path.dirname(OUT.rstrip("/"))): OUT = os.path.abspath("results_full")
os.makedirs(OUT, exist_ok=True)

# ---------------- 多航速汇总
if a.sweep:
    rows = []
    for f in a.sweep:
        for fn in glob.glob(f):
            d = json.load(open(fn)); rows.append(d)
    rows.sort(key=lambda d: d["U"])
    if len(rows) < 2: sys.exit("--sweep 至少要两个航速的 json")
    U = np.array([d["U"] for d in rows]); Fnet = np.array([d["net"]["mean_Fx"] for d in rows]); Fleg = np.array([d["legs_total"]["mean_Fx"] for d in rows]); Fh = np.array([d["hull"]["mean_Fx"] for d in rows])
    print(f"{'U [m/s]':>8}{'四腿推力':>10}{'船体阻力':>10}{'净力':>8}  [mN]")
    for u, fl, fh, fn in zip(U, Fleg, Fh, Fnet): print(f"{u:8.2f}{fl:10.1f}{fh:10.1f}{fn:8.1f}")
    Ueq = None
    for i in range(len(U) - 1):
        if Fnet[i] > 0 >= Fnet[i + 1]: Ueq = U[i] + (U[i + 1] - U[i]) * Fnet[i] / (Fnet[i] - Fnet[i + 1]); break
    print("自航点（净力 = 0）：", f"U ≈ {Ueq:.3f} m/s" if Ueq else "扫掠范围内没有过零（推力始终 %s 阻力）" % ("大于" if Fnet[-1] > 0 else "小于"))
    fig, ax = plt.subplots(figsize=(7.5, 4.6), facecolor=SURF); ax.set_facecolor(SURF)
    ax.plot(U, Fleg, "o-", color="#2a78d6", lw=2, label="四腿周期平均推力")
    ax.plot(U, -Fh, "s-", color="#eb6834", lw=2, label="船体阻力（−F_x）")
    ax.plot(U, Fnet, "^-", color=INK, lw=2.4, label="净力 = 推力 − 阻力")
    ax.axhline(0, color=INK, lw=0.8)
    if Ueq: ax.axvline(Ueq, color=INK2, ls=":", lw=1); ax.text(Ueq, ax.get_ylim()[1] * 0.9, f" 自航点 ≈ {Ueq:.2f} m/s", color=INK2)
    ax.set_xlabel("航速 U [m/s]"); ax.set_ylabel("[mN]"); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.grid(axis="y", color=GRID)
    ax.legend(frameon=False); ax.set_title("整机：推力、阻力与净力随航速", loc="left", fontsize=11)
    plt.tight_layout(); fn = os.path.join(OUT, "fig_full_sweep.png"); plt.savefig(fn, dpi=110, facecolor=SURF); print("图:", fn)
    json.dump(dict(U=U.tolist(), Fleg=Fleg.tolist(), Fhull=Fh.tolist(), Fnet=Fnet.tolist(), Ueq=Ueq), open(os.path.join(OUT, "full_sweep.json"), "w"), indent=1)
    sys.exit(0)

if not a.open and not a.closed: sys.exit("给 --open / --closed 算例目录，或 --sweep")

def read_fo(case, name):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", f"forces_{name}", "*", "force.dat")))
    if not files: return None
    rows = []
    for fn in files:
        d = np.loadtxt(fn, comments="#")
        if d.ndim == 1: d = d[None]
        rows.append(d)
    d = np.concatenate(rows); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
    m = None
    mf = sorted(glob.glob(os.path.join(case, "postProcessing", f"forces_{name}", "*", "moment.dat")))
    if mf:
        mm = np.concatenate([np.loadtxt(f, comments="#").reshape(-1, 10) for f in mf]); mm = mm[np.argsort(mm[:, 0])]; _, i = np.unique(mm[:, 0], return_index=True); m = mm[i]
    return dict(t=d[:, 0], F=d[:, 1:4], P=d[:, 4:7], V=d[:, 7:10], M=m)

cases = {}
for tag, path in (("open", a.open), ("closed", a.closed)):
    if not path: continue
    path = os.path.expanduser(path); g = json.load(open(os.path.join(path, "gait.json")))
    legs = list(g["legs"].keys()); fo = {L: read_fo(path, L) for L in legs}
    fo = {k: v for k, v in fo.items() if v is not None}
    if not fo: print(f"!! {path}: 没有 postProcessing/forces_*/force.dat，跳过"); continue
    hull = read_fo(path, "hull")
    cases[tag] = dict(path=path, gait=g, legs=fo, hull=hull)
    print(f"[{tag}] {path}: 腿 {list(fo)}，船体 {'有' if hull is not None else '无'}，t = {list(fo.values())[0]['t'][-1]:.3f} s")
if not cases: sys.exit("没有可用结果")
g = next(iter(cases.values()))["gait"]; T = g["T"]; DUTY = g["DUTY"]; U = g["U"]
FC = g.get("FAN_CLOSE", -1); tc = DUTY if (FC is None or FC < 0) else FC
t_end = min(min(v["t"][-1] for v in c["legs"].values()) for c in cases.values())
ncyc = int(np.floor(t_end / T + 1e-6))
if ncyc < 1: sys.exit("还没跑满一个周期")
t0, t1 = (ncyc - 1) * T, ncyc * T; tt = np.linspace(t0, t1, 1500); tau = (tt - t0) / T

res = {}   # res[tag][leg] = dict(Fx, Fz, My)
for tag, c in cases.items():
    res[tag] = {}
    for L, d in c["legs"].items():
        r = dict(Fx=np.interp(tt, d["t"], d["F"][:, 0]), Fz=np.interp(tt, d["t"], d["F"][:, 2]), Fy=np.interp(tt, d["t"], d["F"][:, 1]))
        if d["M"] is not None: r["My"] = np.interp(tt, d["M"][:, 0], d["M"][:, 2]); r["Mz"] = np.interp(tt, d["M"][:, 0], d["M"][:, 3])
        res[tag][L] = r
    if c["hull"] is not None:
        d = c["hull"]; r = dict(Fx=np.interp(tt, d["t"], d["F"][:, 0]), Fz=np.interp(tt, d["t"], d["F"][:, 2]), Fy=np.interp(tt, d["t"], d["F"][:, 1]),
                                Px=np.interp(tt, d["t"], d["P"][:, 0]), Vx=np.interp(tt, d["t"], d["V"][:, 0]))
        if d["M"] is not None: r["My"] = np.interp(tt, d["M"][:, 0], d["M"][:, 2]); r["Mz"] = np.interp(tt, d["M"][:, 0], d["M"][:, 3])
        res[tag]["hull"] = r

# ---- 每条腿按自己的相位拼两态
legs = list(next(iter(cases.values()))["legs"].keys())
phase = {L: g["legs"][L]["phase"] for L in legs}
comp = {}
if "open" in res and "closed" in res:
    for L in legs:
        tl = (tau + phase[L]) % 1.0; use_open = tl < tc
        comp[L] = {k: np.where(use_open, res["open"][L][k], res["closed"][L][k]) for k in res["open"][L] if k in res["closed"][L]}
    # 船体：两算例平均（差值就是两态几何对船体阻力的影响范围）
    if "hull" in res["open"] and "hull" in res["closed"]:
        comp["hull"] = {k: 0.5 * (res["open"]["hull"][k] + res["closed"]["hull"][k]) for k in res["open"]["hull"] if k in res["closed"]["hull"]}
        comp["hull_range"] = (res["open"]["hull"]["Fx"].mean() * 1e3, res["closed"]["hull"]["Fx"].mean() * 1e3)
main_tag = "composite" if comp else next(iter(res))
main = comp if comp else res[main_tag]

# ---- 连杆水下阻力（条带法，含前进速度 U 的相对速度），按每条腿的相位/位置算，加到该腿的 F_x、F_z 上
link = {}
if not a.no_linkage:
    try:
        import importlib.util, copy
        HERE = os.path.dirname(os.path.abspath(__file__))
        for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W"]:
            if k in g: os.environ[k] = str(g[k])
        if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
        spec = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(spec); spec.loader.exec_module(ld)
        spec2 = importlib.util.spec_from_file_location("lk", os.path.join(HERE, "linkage_drag.py")); lk = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(lk)
        wl = a.wl if a.wl is not None else g.get("WL", -0.010)
        legs_all = list(next(iter(cases.values()))["legs"].keys())
        main = {k: (dict(v) if isinstance(v, dict) else v) for k, v in main.items()}      # 拷贝一份再加
        for L in legs_all:
            lg = g["legs"][L]
            Lx, Lz, _ = lk.leg_forces(ld, tt, phase=lg["phase"], dx=lg["dx"], wl=wl, U=U)
            link[L] = dict(Fx=Lx, Fz=Lz)
            main[L]["Fx"] = main[L]["Fx"] + Lx; main[L]["Fz"] = main[L]["Fz"] + Lz
        print(f"连杆水下阻力（条带法，水线 {wl:+.3f} m，U = {U} m/s）已加入各腿：" + "  ".join(f"{L} {link[L]['Fx'].mean()*1e3:+.2f} mN" for L in legs_all))
    except Exception as e:
        print("连杆条带法加载失败（腿的力里不含连杆）:", e); link = {}

def tot(d, key):
    return sum(d[L][key] for L in legs if L in d)
summary = dict(U=U, T=T, ncycles=ncyc, tag=main_tag, legs={}, phase=phase, tc=tc, fan_clip=g.get("fan_clip", {}), gaps=g.get("gaps", {}))
print(f"\n第 {ncyc} 周期（t = {t0:.2f}–{t1:.2f} s），U = {U} m/s，+x = 推力，主结果 = {main_tag}")
print(f"{'':<10}{'F̄_x [mN]':>11}{'峰值':>8}{'谷值':>8}{'F̄_z':>8}{'F̄_y':>8}")
for L in legs:
    r = main[L]; d = dict(mean_Fx=float(r["Fx"].mean() * 1e3), peak=float(r["Fx"].max() * 1e3), trough=float(r["Fx"].min() * 1e3), mean_Fz=float(r["Fz"].mean() * 1e3), mean_Fy=float(r["Fy"].mean() * 1e3))
    for tag in res:
        if L in res[tag]: d[f"mean_Fx_{tag}"] = float(res[tag][L]["Fx"].mean() * 1e3)
    if L in link: d["linkage_mean_Fx"] = float(link[L]["Fx"].mean() * 1e3); d["linkage_peak"] = float(np.abs(link[L]["Fx"]).max() * 1e3); d["linkage_mean_Fz"] = float(link[L]["Fz"].mean() * 1e3)
    summary["legs"][L] = d
    print(f"{L:<10}{d['mean_Fx']:11.1f}{d['peak']:8.0f}{d['trough']:8.0f}{d['mean_Fz']:8.1f}{d['mean_Fy']:8.1f}")
Fleg = tot(main, "Fx"); Fzleg = tot(main, "Fz"); Fyleg = tot(main, "Fy")
summary["legs_total"] = dict(mean_Fx=float(Fleg.mean() * 1e3), mean_Fz=float(Fzleg.mean() * 1e3), mean_Fy=float(Fyleg.mean() * 1e3), peak=float(Fleg.max() * 1e3), trough=float(Fleg.min() * 1e3))
if link:
    Lxt = sum(link[L]["Fx"] for L in link); summary["linkage_total"] = dict(mean_Fx=float(Lxt.mean() * 1e3), peak=float(np.abs(Lxt).max() * 1e3), mean_Fz=float(sum(link[L]["Fz"] for L in link).mean() * 1e3), included_in_legs=True)
print(f"{'四腿合计':<10}{Fleg.mean()*1e3:11.1f}{Fleg.max()*1e3:8.0f}{Fleg.min()*1e3:8.0f}{Fzleg.mean()*1e3:8.1f}{Fyleg.mean()*1e3:8.1f}" + (f"   （其中连杆 {summary['linkage_total']['mean_Fx']:+.2f} mN）" if link else ""))
if "hull" in main:
    h = main["hull"]; summary["hull"] = dict(mean_Fx=float(h["Fx"].mean() * 1e3), mean_Fz=float(h["Fz"].mean() * 1e3), peak=float(h["Fx"].max() * 1e3), trough=float(h["Fx"].min() * 1e3))
    if "Px" in h: summary["hull"]["mean_Px"] = float(h["Px"].mean() * 1e3); summary["hull"]["mean_Vx"] = float(h["Vx"].mean() * 1e3)
    if "hull_range" in comp: summary["hull"]["range_open_closed"] = [float(x) for x in comp["hull_range"]]
    print(f"{'船体':<10}{h['Fx'].mean()*1e3:11.1f}{h['Fx'].max()*1e3:8.0f}{h['Fx'].min()*1e3:8.0f}{h['Fz'].mean()*1e3:8.1f}" + (f"   （open 算例 {comp['hull_range'][0]:+.1f} / closed {comp['hull_range'][1]:+.1f}）" if "hull_range" in comp else ""))
    Fnet = Fleg + h["Fx"]; summary["net"] = dict(mean_Fx=float(Fnet.mean() * 1e3), peak=float(Fnet.max() * 1e3), trough=float(Fnet.min() * 1e3))
    print(f"{'净力':<10}{Fnet.mean()*1e3:11.1f}{Fnet.max()*1e3:8.0f}{Fnet.min()*1e3:8.0f}   → {'净推进' if Fnet.mean() > 0 else '净减速'}（U = {U} m/s）")
else:
    Fnet = Fleg; summary["net"] = dict(mean_Fx=float(Fnet.mean() * 1e3))
My = Mz = None
if all("My" in main[L] for L in legs):
    My = tot(main, "My") + (main["hull"]["My"] if "hull" in main and "My" in main["hull"] else 0)
    Mz = tot(main, "Mz") + (main["hull"]["Mz"] if "hull" in main and "Mz" in main["hull"] else 0)
    summary["My"] = dict(mean=float(My.mean() * 1e3), absmax=float(np.abs(My).max() * 1e3)); summary["Mz"] = dict(mean=float(Mz.mean() * 1e3), absmax=float(np.abs(Mz).max() * 1e3))
    print(f"绕船体中点的俯仰力矩 M_y：均值 {My.mean()*1e3:+.2f} mN·m，峰值 {np.abs(My).max()*1e3:.2f} mN·m；偏航力矩 M_z：均值 {Mz.mean()*1e3:+.2f}，峰值 {np.abs(Mz).max()*1e3:.2f} mN·m")
Fy_hull = main["hull"]["Fy"] if "hull" in main else 0 * Fyleg; Fynet = Fyleg + Fy_hull; Fznet = Fzleg + (main["hull"]["Fz"] if "hull" in main else 0)
summary["net"].update(mean_Fy=float(Fynet.mean() * 1e3), peak_Fy=float(np.abs(Fynet).max() * 1e3), mean_Fz=float(Fznet.mean() * 1e3), peak_Fz=float(np.abs(Fznet).max() * 1e3))
print(f"整机三向合力（周期平均 / 峰值）：F_x {summary['net']['mean_Fx']:+.1f} / {max(abs(summary['net'].get('peak',0)), abs(summary['net'].get('trough',0))):.0f}   F_y {Fynet.mean()*1e3:+.1f} / {np.abs(Fynet).max()*1e3:.0f}   F_z {Fznet.mean()*1e3:+.1f} / {np.abs(Fznet).max()*1e3:.0f}  [mN]")
# 逐周期收敛
cyc = []
for n in range(1, ncyc + 1):
    vals = {}
    for tag, c in cases.items():
        k = None
        for L, d in c["legs"].items():
            kk = (d["t"] >= (n - 1) * T) & (d["t"] < n * T)
            vals[f"{tag}_{L}"] = float(np.trapezoid(d["F"][kk, 0], d["t"][kk]) / T * 1e3)
        if c["hull"] is not None:
            d = c["hull"]; kk = (d["t"] >= (n - 1) * T) & (d["t"] < n * T); vals[f"{tag}_hull"] = float(np.trapezoid(d["F"][kk, 0], d["t"][kk]) / T * 1e3)
    cyc.append(dict(n=n, **vals))
summary["cycles"] = cyc
tagname = "_".join(k for k in ("open", "closed") if k in cases)
json.dump(summary, open(os.path.join(OUT, f"full_U{U:.2f}_{tagname}.json"), "w"), indent=1, ensure_ascii=False)

# ---- 图
fig, axs = plt.subplots(2, 3, figsize=(21, 9.5), facecolor=SURF)
for ax in axs.flat:
    ax.set_facecolor(SURF); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.grid(axis="y", color=GRID, lw=0.8)
ax = axs[0, 2]
for L in legs: ax.plot(tau, main[L]["Fz"] * 1e3, color=LEGCOL[L], lw=1.4, label=L)
if "hull" in main: ax.plot(tau, main["hull"]["Fz"] * 1e3, color=LEGCOL["hull"], lw=1.6, label="船体")
ax.plot(tau, Fznet * 1e3, color=INK, lw=2.2, label=f"整机 F_z（均值 {Fznet.mean()*1e3:+.1f}）")
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.set_ylabel("F_z [mN]  (+ 向上)"); ax.legend(frameon=False, fontsize=8.5)
ax.set_title("竖向力 F_z（各腿、船体、整机）", loc="left", fontsize=11)
ax = axs[1, 0]
for L in legs: ax.plot(tau, main[L]["Fy"] * 1e3, color=LEGCOL[L], lw=1.4, label=L)
if "hull" in main: ax.plot(tau, main["hull"]["Fy"] * 1e3, color=LEGCOL["hull"], lw=1.6, label="船体")
ax.plot(tau, Fynet * 1e3, color=INK, lw=2.2, label=f"整机 F_y（均值 {Fynet.mean()*1e3:+.1f}，峰值 {np.abs(Fynet).max()*1e3:.0f}）")
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.set_ylabel("F_y [mN]  (+ 左舷)"); ax.legend(frameon=False, fontsize=8.5)
ax.set_title("横向力 F_y（对角步态的左右交替）", loc="left", fontsize=11)
ax = axs[0, 0]
for L in legs: ax.plot(tau, main[L]["Fx"] * 1e3, color=LEGCOL[L], lw=1.6, label=f"{L}（相位 {phase[L]}，均值 {summary['legs'][L]['mean_Fx']:+.1f}）")
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.set_ylabel("F_x [mN]  (+ 推力)"); ax.legend(frameon=False, fontsize=8.5)
ax.set_title(f"各腿 x 向力（{main_tag}，第 {ncyc} 周期，U = {U} m/s）", loc="left", fontsize=11)
ax = axs[0, 1]
ax.plot(tau, Fleg * 1e3, color="#2a78d6", lw=2, label=f"四腿合计（均值 {Fleg.mean()*1e3:+.1f}）")
if "hull" in main:
    ax.plot(tau, main["hull"]["Fx"] * 1e3, color=LEGCOL["hull"], lw=1.8, label=f"船体（均值 {main['hull']['Fx'].mean()*1e3:+.1f}）")
    ax.plot(tau, Fnet * 1e3, color=INK, lw=2.4, label=f"净力（均值 {Fnet.mean()*1e3:+.1f}）")
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.set_ylabel("F_x [mN]"); ax.legend(frameon=False, fontsize=9)
ax.set_title("整机 x 向力", loc="left", fontsize=11)
ax = axs[1, 1]
if My is not None:
    ax.plot(tau, My * 1e3, color="#eb6834", lw=1.8, label=f"俯仰 M_y（均值 {My.mean()*1e3:+.2f}，峰值 {np.abs(My).max()*1e3:.2f}）")
    ax.plot(tau, Mz * 1e3, color="#2a78d6", lw=1.8, label=f"偏航 M_z（均值 {Mz.mean()*1e3:+.2f}，峰值 {np.abs(Mz).max()*1e3:.2f}）")
    ax.set_ylabel("[mN·m]（参考点 = 船体中点、水线）")
else: ax.text(0.5, 0.5, "没有 moment.dat", ha="center", transform=ax.transAxes)
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.legend(frameon=False, fontsize=9)
ax.set_title("整机力矩（纵摇 / 偏航激励）", loc="left", fontsize=11)
ax = axs[1, 2]
names = legs + (["船体"] if "hull" in main else []) + ["净力"]
vals = [summary["legs"][L]["mean_Fx"] for L in legs] + ([summary["hull"]["mean_Fx"]] if "hull" in main else []) + [summary["net"]["mean_Fx"]]
cols = [LEGCOL[L] for L in legs] + ([LEGCOL["hull"]] if "hull" in main else []) + [INK]
y = np.arange(len(names))[::-1]; ax.barh(y, vals, color=cols, height=0.6)
for yi, v in zip(y, vals): ax.text(v + (1 if v >= 0 else -1), yi, f"{v:+.1f}", va="center", ha="left" if v >= 0 else "right", fontsize=9)
ax.set_yticks(y); ax.set_yticklabels(names); ax.axvline(0, color=INK, lw=0.8); ax.set_xlabel("周期平均 F_x [mN]"); ax.grid(axis="x", color=GRID); ax.grid(axis="y", visible=False)
ax.set_title("周期平均 F_x 分解" + ("（腿含连杆阻力）" if link else ""), loc="left", fontsize=11)
plt.tight_layout(); fn = os.path.join(OUT, f"full_forces_U{U:.2f}_{tagname}.png"); plt.savefig(fn, dpi=110, facecolor=SURF); print("图:", fn)
# 拷原始数据
import shutil
for tag, c in cases.items():
    for L in list(c["legs"]) + (["hull"] if c["hull"] is not None else []):
        src = sorted(glob.glob(os.path.join(c["path"], "postProcessing", f"forces_{L}", "*", "force.dat")))
        if src: shutil.copyfile(src[0], os.path.join(OUT, f"full_force_{L}_{tag}_U{U:.2f}.dat"))
    shutil.copyfile(os.path.join(c["path"], "gait.json"), os.path.join(OUT, f"full_gait_{tag}_U{U:.2f}.json"))
print("结果已拷到:", OUT)
