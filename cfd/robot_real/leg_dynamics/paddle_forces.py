#!/usr/bin/env python3
"""
单腿 CFD 后处理：读 overPimpleDyMFoam 的 forces 输出，做两态脚蹼合成，和准定常模型对比
    python3 paddle_forces.py --open ~/run/leg_open --closed ~/run/leg_closed
    python3 paddle_forces.py --open ~/run/leg_open            # 只有一个算例也行
输出：paddle_forces.png（最后一个周期的 F_x、F_z 历程 + 合成 + 准定常叠加）、控制台里的周期平均表。
合成规则：脚蹼张开的相位（tau < FAN_CLOSE，默认 = DUTY）取 open 算例，其余取 closed 算例；切换瞬态忽略。
"""
import argparse, glob, json, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
def _cjk_font():
    """找一个有中文的字体：系统 Noto CJK，或 WSL 里借 Windows 的微软雅黑/黑体"""
    import glob as _g, warnings
    cands = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]
    cands += _g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + _g.glob("/mnt/c/Windows/Fonts/simhei.ttf") + _g.glob("/mnt/c/Windows/Fonts/simsun.ttc")
    cands += [f for f in font_manager.findSystemFonts() if any(k in f for k in ("wqy", "WenQuanYi", "SourceHanSans", "NotoSansSC"))]
    for f in cands:
        try:
            font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return True
        except Exception: pass
    warnings.filterwarnings("ignore", message="Glyph .* missing"); return False
_cjk_font()
plt.rcParams["axes.unicode_minus"] = False

ap = argparse.ArgumentParser()
ap.add_argument("--open"); ap.add_argument("--closed")
ap.add_argument("--out", default=None)
ap.add_argument("--no-qs", action="store_true", help="不叠加准定常模型")
a = ap.parse_args()
if not a.open and not a.closed: sys.exit("至少给一个算例：--open 或 --closed")

COPY_TO = os.environ.get("COPY_TO", "/mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results")

def read_forces(case):
    files = sorted(glob.glob(os.path.join(case, "postProcessing", "forces", "*", "force.dat")))
    if not files:
        print(f"!! {case}: 还没有 postProcessing/forces/*/force.dat（没跑完或没跑起来），先跳过这个算例"); return None
    rows = []
    for fn in files:                     # 多个启动目录时按时间拼接
        d = np.loadtxt(fn, comments="#")
        if d.ndim == 1: d = d[None]
        rows.append(d)
    d = np.concatenate(rows); d = d[np.argsort(d[:, 0])]
    # 去重（重启会重复时间）
    _, idx = np.unique(d[:, 0], return_index=True); d = d[idx]
    hdr = open(files[0]).read(2000)
    # 列名：v2606 "total_x total_y total_z pressure_x ..."；旧版 "(total) (pressure) (viscous)"
    return dict(t=d[:, 0], total=d[:, 1:4], pressure=d[:, 4:7], viscous=d[:, 7:10])

cases = {}
for tag, path in (("open", a.open), ("closed", a.closed)):
    if path:
        path = os.path.expanduser(path); r_ = read_forces(path)
        if r_ is None: continue
        cases[tag] = r_
        gj = os.path.join(path, "gait.json")
        cases[tag]["gait"] = json.load(open(gj)) if os.path.exists(gj) else {}
        print(f"[{tag}] {path}: {len(cases[tag]['t'])} 个时间点，t = {cases[tag]['t'][0]:.3f} … {cases[tag]['t'][-1]:.3f} s")
if not cases: sys.exit("两个算例都还没有结果。")
g = next(iter(cases.values()))["gait"]
T = g.get("T", 1.25); DUTY = g.get("DUTY", 0.5); FAN_CLOSE = g.get("FAN_CLOSE", -1); U = g.get("U", 0.0)
tc = DUTY if FAN_CLOSE is None or FAN_CLOSE < 0 else FAN_CLOSE
print(f"步态：T = {T} s, DUTY = {DUTY}, 脚蹼张开相位 tau ∈ [0, {tc}), U = {U} m/s")

# 最后一个完整周期
t_end = min(c["t"][-1] for c in cases.values())
n_full = int(np.floor(t_end / T + 1e-6))
if n_full < 1: sys.exit("还没跑满一个周期")
t0, t1 = (n_full - 1) * T, n_full * T
tt = np.linspace(t0, t1, 1000)
def interp(c, key, comp):
    return np.interp(tt, c["t"], c[key][:, comp])
res = {}
for tag, c in cases.items():
    res[tag] = dict(Fx=interp(c, "total", 0), Fz=interp(c, "total", 2), Px=interp(c, "pressure", 0), Vx=interp(c, "viscous", 0))
tau = (tt - t0) / T
fan_open = tau < tc
if "open" in res and "closed" in res:
    comp = {k: np.where(fan_open, res["open"][k], res["closed"][k]) for k in res["open"]}
    res["composite"] = comp

print(f"\n第 {n_full} 个周期（t = {t0:.2f}–{t1:.2f} s）平均，+x = 推力：")
print(f"{'':<12}{'F_x mean [mN]':>14}{'  压差':>8}{'  黏性':>8}{'划水相峰 [mN]':>14}{'回收相谷 [mN]':>14}{'F_z mean [mN]':>14}")
pw = tau < DUTY
summary = {}
for tag, r in res.items():
    m = r["Fx"].mean() * 1e3; summary[tag] = m
    print(f"{tag:<12}{m:14.1f}{r['Px'].mean()*1e3:8.1f}{r['Vx'].mean()*1e3:8.1f}{r['Fx'][pw].max()*1e3:14.1f}{r['Fx'][~pw].min()*1e3:14.1f}{r['Fz'].mean()*1e3:14.1f}")
if "open" in summary and "closed" in summary:
    print(f"\n脚蹼张开带来的收益（open − closed，全程）：{summary['open']-summary['closed']:+.1f} mN；两态合成净推力：{summary.get('composite', float('nan')):.1f} mN")

# 准定常叠加
qs = None
if not a.no_qs:
    try:
        for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W"]:
            if k in g: os.environ[k] = str(g[k])
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import leg_dynamics as ld
        r = ld.simulate(U); m = r["m"]; tq = (r["t"][m] - r["t"][m][0]) / ld.T_PER
        qs = dict(tau=tq, Fx=r["F"][m, 0] * 1e3, Fd=(r["F"][m, 0] - r["Fam"][m, 0]) * 1e3)
        print(f"准定常模型（同参数）：净推力 含附加质量 {qs['Fx'].mean():.1f} mN / 仅压差 {qs['Fd'].mean():.1f} mN")
    except Exception as e:
        print("（准定常叠加失败：", e, "）")

# ---- 图
fig, axs = plt.subplots(1, 2, figsize=(16, 5.5))
ax = axs[0]
cols = {"open": "#c0392b", "closed": "#2980b9", "composite": "k"}
for tag, r in res.items():
    ax.plot(tau, r["Fx"] * 1e3, color=cols[tag], lw=2.5 if tag == "composite" else 1.6, label=f"CFD {tag}  (mean {summary[tag]:+.1f} mN)")
if qs is not None:
    ax.plot(qs["tau"], qs["Fx"], "--", color="gray", lw=1.2, label=f"准定常 含附加质量 (mean {qs['Fx'].mean():+.1f})")
    ax.plot(qs["tau"], qs["Fd"], ":", color="gray", lw=1.2, label=f"准定常 仅压差 (mean {qs['Fd'].mean():+.1f})")
ax.axvspan(0, DUTY, color="#c0392b", alpha=.07); ax.axvspan(0, tc, color="#f39c12", alpha=.07)
ax.axhline(0, color="k", lw=.8); ax.set_xlabel("t / T（最后一个周期）"); ax.set_ylabel("F_x [mN]  (+ 推力)")
ax.set_title(f"单腿 x 向力 · CFD (overset, {g.get('turb','laminar')})  U = {U} m/s"); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.text(DUTY / 2, ax.get_ylim()[1] * .92, "划水相", ha="center", color="#c0392b")
ax = axs[1]
for tag, r in res.items():
    ax.plot(tau, r["Fz"] * 1e3, color=cols[tag], lw=1.6, label=f"{tag} F_z")
    if tag != "composite":
        ax.plot(tau, r["Vx"] * 1e3, ":", color=cols[tag], lw=1, label=f"{tag} 黏性 F_x")
ax.axvspan(0, DUTY, color="#c0392b", alpha=.07); ax.axhline(0, color="k", lw=.8)
ax.set_xlabel("t / T"); ax.set_ylabel("[mN]"); ax.set_title("竖向力 F_z 与黏性分量"); ax.grid(alpha=.3); ax.legend(fontsize=8)
plt.tight_layout()
first = os.path.expanduser(a.open or a.closed)
out = a.out or os.path.join(os.path.dirname(first), "paddle_forces.png")
plt.savefig(out, dpi=100); print("图已存:", out)
# 顺手拷到共享文件夹（Windows 里能直接看到）
if COPY_TO and os.path.isdir(os.path.dirname(COPY_TO.rstrip("/"))):
    import shutil
    os.makedirs(COPY_TO, exist_ok=True)
    tagname = "_".join(k for k in ("open", "closed") if k in cases)
    dst = os.path.join(COPY_TO, f"paddle_forces_{tagname}.png"); shutil.copyfile(out, dst)
    for tag, path in (("open", a.open), ("closed", a.closed)):
        if tag in cases:
            src = glob.glob(os.path.join(os.path.expanduser(path), "postProcessing", "forces", "*", "force.dat"))[0]
            shutil.copyfile(src, os.path.join(COPY_TO, f"force_{tag}.dat"))
            gj = os.path.join(os.path.expanduser(path), "gait.json")
            if os.path.exists(gj): shutil.copyfile(gj, os.path.join(COPY_TO, f"gait_{tag}.json"))
            for fo, nm in (("forcesLinks", "force_links"), ):          # --leg 算例：连杆段单独的力
                fl = glob.glob(os.path.join(os.path.expanduser(path), "postProcessing", fo, "*", "force.dat"))
                if fl: shutil.copyfile(fl[0], os.path.join(COPY_TO, f"{nm}_{tag}.dat"))
            ml = glob.glob(os.path.join(os.path.expanduser(path), "postProcessing", "forces", "*", "moment.dat"))
            if ml: shutil.copyfile(ml[0], os.path.join(COPY_TO, f"moment_{tag}.dat"))
    print(f"已拷到共享文件夹: {COPY_TO}  （图 + force_*.dat + gait_*.json）")
