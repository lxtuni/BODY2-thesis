#!/usr/bin/env python3
"""
单腿动态 CFD 结果分析（报告用）：读 results/force_{open,closed}.dat + gait_*.json，
输出 report_data.json 和一组图（时间历程、分相位分解、逐周期收敛、F_z、准定常对照）。
    python3 analyze_leg.py <results_dir> [out_dir]
"""
import os, sys, json, glob, importlib.util
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
BLUE, RED, GRAY, INK, INK2, SURF, ORANGE = "#2a78d6", "#e34948", "#a8a7a1", "#0b0b0b", "#52514e", "#fcfcfb", "#eb6834"

R = os.path.abspath(sys.argv[1]); OUT = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else R
os.makedirs(OUT, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))

def load(tag):
    fn = os.path.join(R, f"force_{tag}.dat")
    if not os.path.exists(fn): return None
    d = np.loadtxt(fn, comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
    g = json.load(open(os.path.join(R, f"gait_{tag}.json"))) if os.path.exists(os.path.join(R, f"gait_{tag}.json")) else {}
    mfn = os.path.join(R, f"moment_{tag}.dat"); M = None
    if os.path.exists(mfn):
        m = np.loadtxt(mfn, comments="#"); m = m[np.argsort(m[:, 0])]; _, i = np.unique(m[:, 0], return_index=True); M = m[i]
    return dict(t=d[:, 0], F=d[:, 1:4], P=d[:, 4:7], V=d[:, 7:10], gait=g, M=M)

cases = {k: v for k in ("open", "closed") if (v := load(k)) is not None}
if not cases: sys.exit("results 里没有 force_open.dat / force_closed.dat")
g = next(iter(cases.values()))["gait"]; T = g.get("T", 1.25); DUTY = g.get("DUTY", 0.5); U = g.get("U", 0.0)
FC = g.get("FAN_CLOSE", -1); tc = DUTY if (FC is None or FC < 0) else FC
ncyc = int(np.floor(min(c["t"][-1] for c in cases.values()) / T + 1e-6))
t0, t1 = (ncyc - 1) * T, ncyc * T
tt = np.linspace(t0, t1, 2000); tau = (tt - t0) / T
def smooth(y, x, w=0.02):
    out = np.empty_like(y)
    for i in range(len(y)):
        k = (x >= x[i] - w / 2) & (x <= x[i] + w / 2); out[i] = y[k].mean()
    return out
res = {}
for tag, c in cases.items():
    r = {k: np.interp(tt, c["t"], c[a][:, i]) for k, (a, i) in {"Fx": ("F", 0), "Fy": ("F", 1), "Fz": ("F", 2), "Px": ("P", 0), "Vx": ("V", 0)}.items()}
    r["Fx_s"] = smooth(r["Fx"], tt); r["Fz_s"] = smooth(r["Fz"], tt)
    if c["M"] is not None: r["My"] = np.interp(tt, c["M"][:, 0], c["M"][:, 2]); r["My_s"] = smooth(r["My"], tt)
    res[tag] = r
fan_open = tau < tc
if "open" in res and "closed" in res:
    res["composite"] = {k: np.where(fan_open, res["open"][k], res["closed"][k]) for k in res["open"]}

# ---- phases
rec_a = DUTY; tr = g.get("SWITCH_FRAC", 0.3) * (1 - DUTY)
phases = [("划水相（伸展，扇张开）", 0.0, DUTY), ("收拢切换 φ 105→40°", DUTY, DUTY + tr), ("蜷缩前扫", DUTY + tr, 1 - tr), ("伸展切换 φ 40→105°", 1 - tr, 1.0)]
def contrib(fx):
    out = []
    for name, a, b in phases:
        k = (tau >= a) & (tau < b); out.append(float(np.trapezoid(fx[k], tt[k]) / T * 1e3))
    return out

data = dict(T=T, DUTY=DUTY, U=U, tc=tc, ncycles=ncyc, gait=g, cases={}, phases=[p[0] for p in phases])
for tag, r in res.items():
    pw = tau < DUTY
    d = dict(mean_Fx=float(r["Fx"].mean() * 1e3), mean_Px=float(r["Px"].mean() * 1e3), mean_Vx=float(r["Vx"].mean() * 1e3),
             mean_Fz=float(r["Fz"].mean() * 1e3), mean_Fy=float(r["Fy"].mean() * 1e3),
             peak_power_raw=float(r["Fx"][pw].max() * 1e3), peak_power_s=float(r["Fx_s"][pw].max() * 1e3),
             trough_raw=float(r["Fx"].min() * 1e3), trough_s=float(r["Fx_s"].min() * 1e3),
             peak_rec_s=float(r["Fx_s"][~pw].max() * 1e3), Fz_absmax_s=float(np.abs(r["Fz_s"]).max() * 1e3),
             contrib=contrib(r["Fx"]))
    if "My_s" in r: d["My_absmax"] = float(np.abs(r["My_s"]).max() * 1e3); d["My_mean"] = float(r["My"].mean() * 1e3)
    data["cases"][tag] = d
# per-cycle convergence (raw cases)
for tag, c in cases.items():
    cyc = []
    for n in range(1, ncyc + 1):
        k = (c["t"] >= (n - 1) * T) & (c["t"] < n * T)
        cyc.append(dict(n=n, Fx=float(np.trapezoid(c["F"][k, 0], c["t"][k]) / T * 1e3), Fz=float(np.trapezoid(c["F"][k, 2], c["t"][k]) / T * 1e3)))
    data["cases"][tag]["cycles"] = cyc
    data["cases"][tag]["nsteps"] = int(len(c["t"])); data["cases"][tag]["dt_mean_ms"] = float(np.mean(np.diff(c["t"])) * 1e3)
    k = (tau > 0.1) & (tau < 0.3); data["cases"][tag]["jitter_std"] = float(np.std(res[tag]["Fx"][k] - res[tag]["Fx_s"][k]) * 1e3)

# ---- quasi-steady
qs = None
try:
    for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
        if k in g: os.environ[k] = str(g[k])
    if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
    spec = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(spec); spec.loader.exec_module(ld)
    r = ld.simulate(U); m = r["m"]; tq = (r["t"][m] - r["t"][m][0]) / T
    qs = dict(tau=tq, Fx=r["F"][m, 0] * 1e3, Fd=(r["F"][m, 0] - r["Fam"][m, 0]) * 1e3)
    qsc = []
    for name, a, b in phases:
        k = (tq >= a) & (tq < b); qsc.append(float(np.trapezoid(qs["Fd"][k], r["t"][m][k]) / T))
    data["qs"] = dict(mean_total=float(qs["Fx"].mean()), mean_drag=float(qs["Fd"].mean()), contrib_drag=qsc, peak_drag=float(qs["Fd"].max()))
    # servo demand from kinematics
    tipv = np.linalg.norm(np.gradient(r["tip"][m], r["t"][m], axis=0), axis=1)
    data["kin"] = dict(tip_speed_max=float(tipv.max()), tip_speed_power_mean=float(tipv[r["power"][m]].mean()),
                       tip_x=[float(r["tip"][m][:, 0].min()), float(r["tip"][m][:, 0].max())], tip_z=[float(r["tip"][m][:, 1].min()), float(r["tip"][m][:, 1].max())],
                       th2_range=[float(np.degrees(r["th2"][m].min())), float(np.degrees(r["th2"][m].max()))])
except Exception as e:
    print("准定常对照失败:", e)

# ---- --leg 算例：连杆段单独的 CFD 力（force_links_<tag>.dat），桨 = 总力 − 连杆
for tag, c in cases.items():
    fl = os.path.join(R, f"force_links_{tag}.dat")
    if not os.path.exists(fl): continue
    d = np.loadtxt(fl, comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
    Lx = np.interp(tt, d[:, 0], d[:, 1]); Lz = np.interp(tt, d[:, 0], d[:, 3]); Lx_s = smooth(Lx, tt)
    data["cases"][tag]["links_cfd"] = dict(mean_Fx=float(Lx.mean() * 1e3), mean_Fz=float(Lz.mean() * 1e3), peak_s=float(np.abs(Lx_s).max() * 1e3),
                                           peak_raw=float(np.abs(Lx).max() * 1e3), contrib=contrib(Lx),
                                           share_mean=float(Lx.mean() / res[tag]["Fx"].mean()) if res[tag]["Fx"].mean() else None)
    res[tag]["Lx_s"] = Lx_s
    print(f"[{tag}] CFD 连杆段（--leg）：均值 {Lx.mean()*1e3:+.2f} mN（占总力 {100*data['cases'][tag]['links_cfd']['share_mean'] or 0:.1f}%），平滑峰值 {np.abs(Lx_s).max()*1e3:.1f} mN")
json.dump(data, open(os.path.join(OUT, "report_data.json"), "w"), indent=1, ensure_ascii=False)

# ================= figures =================
cols = {"open": ORANGE, "closed": BLUE, "composite": INK}
# 1) time histories
fig, axs = plt.subplots(1, 2, figsize=(16, 5.8), facecolor=SURF)
ax = axs[0]; ax.set_facecolor(SURF)
for tag in ("open", "closed", "composite"):
    if tag not in res: continue
    r = res[tag]; lw = 2.6 if tag == "composite" else 1.5
    ax.plot(tau, r["Fx"] * 1e3, color=cols[tag], lw=0.5, alpha=0.2)
    ax.plot(tau, r["Fx_s"] * 1e3, color=cols[tag], lw=lw, label={"open": "全程张开 open", "closed": "全程收拢 closed", "composite": "两态合成 composite"}[tag] + f"（均值 {data['cases'][tag]['mean_Fx']:+.1f} mN）")
if qs is not None: ax.plot(qs["tau"], qs["Fd"], "--", color=GRAY, lw=1.3, label=f"准定常叶素 仅压差（{qs['Fd'].mean():+.1f}）")
for tag in ("open", "closed"):
    if tag in res and "Lx_s" in res[tag]:
        ax.plot(tau, res[tag]["Lx_s"] * 1e3, ":", color=cols[tag], lw=1.6, label=f"{tag} 其中连杆段 CFD（{data['cases'][tag]['links_cfd']['mean_Fx']:+.2f}）")
for name, a, b in phases: ax.axvline(a, color=INK2, lw=0.6, alpha=0.35); ax.text((a + b) / 2, 0.97, name, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8.5, color=INK2)
ax.axvspan(0, tc, color=ORANGE, alpha=0.06)
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel(f"t / T（第 {ncyc} 周期）", color=INK); ax.set_ylabel("单腿 x 向力 F_x [mN]  (+ = 推力)", color=INK)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e6e5e1", lw=0.8); ax.legend(frameon=False, fontsize=8.5, loc="lower right"); ax.set_title(f"x 向力时间历程（20 ms 平滑；淡线为原始）  U = {U} m/s", loc="left", fontsize=11)
ax = axs[1]; ax.set_facecolor(SURF)
for tag in ("open", "closed", "composite"):
    if tag not in res: continue
    ax.plot(tau, res[tag]["Fz_s"] * 1e3, color=cols[tag], lw=2.2 if tag == "composite" else 1.4, label=tag)
for name, a, b in phases: ax.axvline(a, color=INK2, lw=0.6, alpha=0.35)
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T", color=INK); ax.set_ylabel("竖向力 F_z [mN]  (+ = 向上)", color=INK)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e6e5e1", lw=0.8); ax.legend(frameon=False, fontsize=9); ax.set_title("竖向力（20 ms 平滑）", loc="left", fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_forces.png"), dpi=110, facecolor=SURF); plt.close()

# 2) phase contribution bars
tags = [t for t in ("open", "closed", "composite") if t in res]
fig, ax = plt.subplots(figsize=(11, 5.2), facecolor=SURF); ax.set_facecolor(SURF)
y = np.arange(len(phases))[::-1]; h = 0.8 / len(tags)
for j, tag in enumerate(tags):
    c = data["cases"][tag]["contrib"]; yy = y + (j - (len(tags) - 1) / 2) * h
    ax.barh(yy, c, height=h * 0.92, color=cols[tag], edgecolor=SURF, linewidth=1.5, label=tag)
    for yi, v in zip(yy, c): ax.text(v + (0.6 if v >= 0 else -0.6), yi, f"{v:+.1f}", va="center", ha="left" if v >= 0 else "right", fontsize=8.5, color=INK)
if qs is not None: ax.scatter(data["qs"]["contrib_drag"], y, marker="|", s=260, color=INK2, zorder=3, label="准定常同相位（仅压差）")
ax.axvline(0, color=INK, lw=0.8); ax.set_yticks(y); ax.set_yticklabels([p[0] for p in phases], fontsize=9.5)
ax.set_xlabel("对周期平均推力的贡献 [mN]", color=INK)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="x", color="#e6e5e1", lw=0.8); ax.legend(frameon=False, fontsize=9, loc="lower right")
tot = "  ".join(f"{t} {data['cases'][t]['mean_Fx']:+.1f}" for t in tags)
ax.set_title(f"分相位冲量 → 周期平均推力 [mN]   合计：{tot}", loc="left", fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_phases.png"), dpi=110, facecolor=SURF); plt.close()

# 3) cycle convergence
fig, ax = plt.subplots(figsize=(7, 4), facecolor=SURF); ax.set_facecolor(SURF)
for tag, c in cases.items():
    cyc = data["cases"][tag]["cycles"]; ax.plot([q["n"] for q in cyc], [q["Fx"] for q in cyc], "o-", color=cols[tag], lw=1.8, ms=7, label=f"{tag} F_x")
    for q in cyc: ax.text(q["n"], q["Fx"] + 0.8, f"{q['Fx']:.1f}", ha="center", fontsize=8.5, color=INK)
ax.set_xlabel("周期序号", color=INK); ax.set_ylabel("周期平均 F_x [mN]", color=INK); ax.set_xticks(range(1, ncyc + 1))
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e6e5e1", lw=0.8); ax.legend(frameon=False, fontsize=9); ax.set_title("逐周期平均推力（周期性检验）", loc="left", fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_cycles.png"), dpi=110, facecolor=SURF); plt.close()


# 4) paddle two-state outlines (to scale)
try:
    S, H, SW, HW, R_, ANG = ld.STEM, ld.HEAD, ld.STEM_W, ld.HEAD_W, ld.FAN_R, np.deg2rad(ld.FAN_DEG)
    fig, axs = plt.subplots(1, 2, figsize=(7.5, 4.2), facecolor=SURF, sharey=True)
    kk = np.linspace(0, 1, 24); sh = [(HW / 2 * np.sin(np.pi / 2 * q), S + 0.025 * (1 - np.cos(np.pi / 2 * q))) for q in kk]
    right = [(w, z) for w, z in sh if w > SW / 2]
    closed = [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + right + [(HW / 2, S + H), (-HW / 2, S + H)] + [(-w, z) for w, z in right[::-1]] + [(-SW / 2, S)]
    th = np.linspace(ANG / 2, -ANG / 2, 48); fan = [(R_ * np.sin(x), S + R_ * np.cos(x)) for x in th]
    opn = [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + fan + [(-SW / 2, S)]
    for ax, poly, col, name in ((axs[0], closed, BLUE, f"收拢态 closed（CAD 桨头，宽 {HW*1000:.0f} mm）"), (axs[1], opn, ORANGE, f"张开态 open（扇 {ld.FAN_DEG:.0f}°，r = {R_*1000:.0f} mm）")):
        P = np.array(poly) * 1000; ax.fill(P[:, 0], -P[:, 1], color=col, alpha=0.75, lw=0)
        ax.set_aspect("equal"); ax.set_facecolor(SURF); ax.set_title(name, fontsize=10, loc="left"); ax.set_xlabel("横向 y [mm]")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.grid(color="#e6e5e1", lw=0.8); ax.plot(0, 0, "o", color=INK, ms=5); ax.text(2, 2, "销孔 E", fontsize=8.5, color=INK2)
    axs[0].set_ylabel("沿桨轴距销孔 [mm]（向下为负）")
    A_c = S * SW + HW * H - (1 - np.pi / 4) * HW * 0.025 / 2; A_o = S * SW + 0.5 * R_**2 * ANG
    fig.suptitle(f"脚蹼两态几何（面积 {A_c*1e4:.1f} → {A_o*1e4:.1f} cm²，×{A_o/A_c:.2f}；厚 {ld.THK*1000:.1f} mm）", fontsize=11)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_paddles.png"), dpi=110, facecolor=SURF); plt.close()
    data["paddle_area_cm2"] = [float(A_c * 1e4), float(A_o * 1e4)]
except Exception as e:
    print("桨板轮廓图失败:", e)

# 5) CFD frame montages from the mp4s (if present): 4 instants of the last cycle
try:
    import subprocess, shutil
    from PIL import Image
    for view in ("side", "3d"):
        cands = glob.glob(os.path.join(R, f"paddle_{view}_leg_open*.mp4")) + glob.glob(os.path.join(R, f"paddle_{view}.mp4"))
        if not cands: continue
        mp4 = cands[0]
        dur = float(subprocess.check_output(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "default=nw=1:nk=1", mp4]).decode().strip())
        nfr = int(round(dur * 12))                                     # FPS=12 → 帧数；每周期 writes=20 帧 → 3 周期 60 帧
        per = nfr / 3.0 if nfr >= 40 else nfr
        picks = [(0.2, "划水相中段"), (0.55, "收拢切换"), (0.75, "蜷缩前扫"), (0.92, "伸展切换")]
        ims = []
        for tau_, lab in picks:
            fi = int(min(nfr - 1, round(2 * per + tau_ * per))) if nfr >= 40 else int(round(tau_ * per))
            fn = os.path.join(OUT, f"_fr_{view}_{int(tau_*100)}.png")
            subprocess.run(["ffmpeg", "-v", "error", "-y", "-i", mp4, "-vf", f"select=eq(n\\,{fi})", "-vframes", "1", fn], check=True)
            im = Image.open(fn).convert("RGB"); ims.append((im, f"t/T = {tau_:.2f}  {lab}"))
        w, h = ims[0][0].size; sc = 0.5; w, h = int(w * sc), int(h * sc)
        from PIL import ImageDraw, ImageFont
        fnt = None
        for fp in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f] + glob.glob("/mnt/c/Windows/Fonts/msyh.ttc") + glob.glob("/mnt/c/Windows/Fonts/simhei.ttf"):
            try: fnt = ImageFont.truetype(fp, 18); break
            except Exception: pass
        sheet = Image.new("RGB", (w * 2, h * 2), "white")
        for i, (im, lab) in enumerate(ims):
            im2 = im.resize((w, h)); d = ImageDraw.Draw(im2); d.rectangle((0, h - 28, 300, h), fill=(252, 252, 251))
            d.text((8, h - 25), lab, fill=(11, 11, 11), font=fnt)
            sheet.paste(im2, ((i % 2) * w, (i // 2) * h))
        sheet.save(os.path.join(OUT, f"fig_cfd_{view}.png")); data[f"fig_cfd_{view}"] = True
        for f in glob.glob(os.path.join(OUT, f"_fr_{view}_*.png")): os.remove(f)
except Exception as e:
    print("CFD 截帧失败:", e)
json.dump(data, open(os.path.join(OUT, "report_data.json"), "w"), indent=1, ensure_ascii=False)

# ---- 连杆条带法估算（linkage_drag.py）→ linkage.json + fig_linkage.png
try:
    import subprocess
    wl = str(g.get("WL", -0.010))
    subprocess.run([sys.executable, os.path.join(HERE, "linkage_drag.py"), R, "--out", OUT, "--wl", wl], check=True)
except Exception as e:
    print("连杆估算失败:", e)

# ---- console summary
print(f"周期 T={T}s DUTY={DUTY} U={U}  用第 {ncyc} 周期")
for tag in tags:
    d = data["cases"][tag]
    print(f"[{tag:<9}] 均值 F_x {d['mean_Fx']:+6.1f} mN (压差 {d['mean_Px']:+.1f} 黏性 {d['mean_Vx']:+.2f})  F_z {d['mean_Fz']:+.1f}  划水峰(平滑) {d['peak_power_s']:.0f}  谷(平滑) {d['trough_s']:.0f}  分相位 {[round(x,1) for x in d['contrib']]}")
if qs is not None: print(f"[准定常   ] 含附加质量 {data['qs']['mean_total']:+.1f}  仅压差 {data['qs']['mean_drag']:+.1f}")
print("图:", ", ".join(f for f in ("fig_forces.png", "fig_phases.png", "fig_cycles.png")), "→", OUT)
