#!/usr/bin/env python3
"""
船体形状比选——在 hull_compare.py 汇总的基础上做二次分析，供报告使用
    python3 analyze_hulls.py [--res results] [--hydro hull_hydrostatics.json] [--eqdisp m2_eqdisp.json]
输入：results/hull_resistance.json、results/force_hull_*.dat、hull_hydrostatics.json、（可选）results/force_hull_M2_roundbottom_U02_eqdisp.dat
输出：results/hull_analysis.json、fig_hull_components.png、fig_hull_speed.png、fig_hull_budget.png、fig_hull_fairness.png
"""
import os, sys, glob, json, argparse
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
HULLS = [("M1_current", "M1 现船体（平底）"), ("M2_roundbottom", "M2 圆滑底"), ("M3_streamlined", "M3 流线型圆舭")]
SHORT = {"M1_current": "M1", "M2_roundbottom": "M2", "M3_streamlined": "M3"}
COL = {"M1_current": INK, "M2_roundbottom": BLUE, "M3_streamlined": GREEN}
RHO, NU, G = 998.8, 1.09e-6, 9.81           # DTCHull 模板物性
R_STATIC_02 = 77.6                          # 静止浮态报告：四腿在水下、U=0.2 m/s 的整机阻力 [mN]
R_HULL_ASSUMED_V11 = 25.0                   # v1.1 单腿报告中按迎流面积比例估算的船体阻力 [mN]

ap = argparse.ArgumentParser()
ap.add_argument("--res", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"))
ap.add_argument("--hydro", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "hull_hydrostatics.json"))
ap.add_argument("--eqdisp", default=None, help="M2 同排水量水线的静水力学 json（eqdisp.py 输出）")
a = ap.parse_args(); RES = a.res
R = json.load(open(os.path.join(RES, "hull_resistance.json"))); hyd = json.load(open(a.hydro))
res = R["results"]; speeds = sorted({float(U) for h in res.values() for U in h})

def hist(name):
    fn = os.path.join(RES, f"force_{name}.dat")
    if not os.path.exists(fn): return None
    d = np.loadtxt(fn, comments="#"); return dict(t=d[:, 0], Fx=-d[:, 1] * 1e3, Px=-d[:, 4] * 1e3, Vx=-d[:, 7] * 1e3, Fz=d[:, 3])

def windows(h):
    n = len(h["t"]); out = {}
    for f in (0.2, 0.4):
        k = int(round(n * (1 - f)))
        out[f"w{int(f*100)}"] = dict(Fx=float(h["Fx"][k:].mean()), Px=float(h["Px"][k:].mean()), Vx=float(h["Vx"][k:].mean()), std=float(h["Fx"][k:].std()),
                                     std_P=float(h["Px"][k:].std()), std_V=float(h["Vx"][k:].std()), lo=float(h["Fx"][k:].min()), hi=float(h["Fx"][k:].max()))
    k = int(round(n * 0.6)); x = h["Fx"][k:] - h["Fx"][k:].mean(); zc = np.where(np.diff(np.sign(x)) != 0)[0]
    out["period_steps"] = float(2 * np.diff(h["t"][k:][zc]).mean()) if len(zc) > 2 else None
    out["trend_per100"] = float(np.polyfit(h["t"][k:], h["Fx"][k:], 1)[0] * 100)
    return out

A = {"hulls": {}, "speeds": speeds, "const": dict(rho=RHO, nu=NU, g=G, R_static_02=R_STATIC_02, R_hull_assumed_v11=R_HULL_ASSUMED_V11)}
for h, lab in HULLS:
    A["hulls"][h] = dict(label=lab, hydro=hyd[h], U={})
    for U in speeds:
        r = res[h].get(f"{U:.1f}")
        if not r: continue
        H = hist(os.path.basename(r["case"])); w = windows(H) if H else {}
        L = hyd[h]["LWL_mm"] * 1e-3; S = hyd[h]["S_cm2"] * 1e-4; Af = hyd[h]["Amax_cm2"] * 1e-4; q = 0.5 * RHO * U ** 2
        Re = U * L / NU; Fr = U / np.sqrt(G * L)
        CF_ittc = 0.075 / (np.log10(Re) - 2) ** 2; CF_blasius = 1.328 / np.sqrt(Re)
        A["hulls"][h]["U"][f"{U:.1f}"] = dict(Fx=r["Fx_mN"], Px=r["Px_mN"], Vx=r["Vx_mN"], Fz=r["Fz_N"], disp_g_cfd=r["disp_g_cfd"], P_mW=r["P_mW"],
                                             Re=Re, Fr=Fr, CT=r["Fx_mN"] * 1e-3 / (q * S), CP=r["Px_mN"] * 1e-3 / (q * S), CF=r["Vx_mN"] * 1e-3 / (q * S),
                                             CD_front=r["Fx_mN"] * 1e-3 / (q * Af), CF_ittc=CF_ittc, CF_blasius=CF_blasius, CF_over_ittc=r["Vx_mN"] * 1e-3 / (q * S) / CF_ittc,
                                             R_visc_ittc_mN=CF_ittc * q * S * 1e3, R_visc_blasius_mN=CF_blasius * q * S * 1e3, win=w,
                                             share_of_static=r["Fx_mN"] / R_STATIC_02 * (U / 0.2) ** (-1.95) if U == 0.2 else None)
    # 速度指数（两点）
    u = sorted(A["hulls"][h]["U"].keys())
    if len(u) == 2:
        r1, r2 = A["hulls"][h]["U"][u[0]], A["hulls"][h]["U"][u[1]]; du = np.log(float(u[1]) / float(u[0]))
        A["hulls"][h]["n_total"] = float(np.log(r2["Fx"] / r1["Fx"]) / du); A["hulls"][h]["n_press"] = float(np.log(r2["Px"] / r1["Px"]) / du); A["hulls"][h]["n_visc"] = float(np.log(r2["Vx"] / r1["Vx"]) / du)
        A["hulls"][h]["n_total_w40"] = float(np.log(r2["win"]["w40"]["Fx"] / r1["win"]["w40"]["Fx"]) / du) if r1.get("win") and r2.get("win") else None

# 相对 M1
for h, _ in HULLS:
    for U in A["hulls"][h]["U"]:
        ref = A["hulls"]["M1_current"]["U"][U]
        A["hulls"][h]["U"][U]["rel_M1"] = A["hulls"][h]["U"][U]["Fx"] / ref["Fx"] - 1
        A["hulls"][h]["U"][U]["rel_M1_w40"] = A["hulls"][h]["U"][U]["win"]["w40"]["Fx"] / ref["win"]["w40"]["Fx"] - 1

# ---- M2 同排水量修正（一阶估算：摩擦 ∝ 湿表面，压差 ∝ 最大剖面）+ 若有 CFD 校核算例则读入
eq = {}
if a.eqdisp and os.path.exists(a.eqdisp):
    E = json.load(open(a.eqdisp)); keys = sorted(E.keys(), key=float); e0, e1 = E[keys[0]], E[keys[-1]]    # -10.00 与 -7.77
    kS = e1["S_cm2"] / e0["S_cm2"]; kA = e1["Amax_cm2"] / e0["Amax_cm2"]
    eq = dict(WL_mm=float(keys[-1]), hydro_eq=e1, k_S=kS, k_Amax=kA, est={})
    for U in A["hulls"]["M2_roundbottom"]["U"]:
        r = A["hulls"]["M2_roundbottom"]["U"][U]; ref = A["hulls"]["M1_current"]["U"][U]
        est = r["Vx"] * kS + r["Px"] * kA
        eq["est"][U] = dict(Fx=est, Vx=r["Vx"] * kS, Px=r["Px"] * kA, rel_M1=est / ref["Fx"] - 1)
    H = hist("hull_M2_roundbottom_U02_eqdisp")
    if H is not None:
        w = windows(H); eq["cfd"] = dict(U=0.2, Fx=w["w20"]["Fx"], Px=w["w20"]["Px"], Vx=w["w20"]["Vx"], Fz=float(H["Fz"][int(0.8 * len(H["Fz"])):].mean()), win=w,
                                          rel_M1=w["w20"]["Fx"] / A["hulls"]["M1_current"]["U"]["0.2"]["Fx"] - 1, steps=int(H["t"][-1]))
        eq["cfd"]["disp_g_cfd"] = eq["cfd"]["Fz"] / G * 1e3
A["M2_eqdisp"] = eq

def conv(o):
    if isinstance(o, dict): return {k: conv(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [conv(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return float(o)
    return o
json.dump(conv(A), open(os.path.join(RES, "hull_analysis.json"), "w"), indent=1, ensure_ascii=False)

# ---------------- 图 1：U=0.2 分量收敛历程（压差在摆、摩擦已收敛） ----------------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.4), facecolor=SURF, sharey=True)
for ax, (h, lab) in zip(axs, HULLS):
    ax.set_facecolor(SURF); r = A["hulls"][h]["U"]["0.2"]; H = hist(os.path.basename(res[h]["0.2"]["case"]))
    m = H["t"] >= 300
    ax.plot(H["t"][m], H["Fx"][m], color=COL[h], lw=1.8, label="总阻力 (gesamt)")
    ax.plot(H["t"][m], H["Px"][m], color=ORANGE, lw=1.4, label="压差 (Druck)")
    ax.plot(H["t"][m], H["Vx"][m], color=GRAY, lw=1.4, label="摩擦 (Reibung)")
    ax.axvspan(1200, 1500, color=GRAY, alpha=0.10, lw=0); ax.axvspan(900, 1200, color=GRAY, alpha=0.05, lw=0)
    w = r["win"]; ax.hlines(w["w20"]["Fx"], 1200, 1500, color=COL[h], lw=1.0, ls="--"); ax.hlines(w["w40"]["Fx"], 900, 1500, color=COL[h], lw=0.8, ls=":")
    ax.text(1500, w["w20"]["Fx"], f"  20%: {w['w20']['Fx']:.1f}", va="center", ha="left", fontsize=8.5, color=INK, clip_on=False)
    ax.text(1500, w["w40"]["Fx"] + (0.9 if w["w40"]["Fx"] > w["w20"]["Fx"] else -0.9), f"  40%: {w['w40']['Fx']:.1f}", va="center", ha="left", fontsize=8.5, color=INK2, clip_on=False)
    ax.set_title(f"{lab}   U = 0.2 m/s", loc="left", fontsize=10.5); ax.set_xlabel("LTS 伪时间步"); ax.set_xlim(300, 1560)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(color=GRID, lw=0.8)
    ax.text(320, 0.6, f"±1σ(尾段 20%) = {w['w20']['std']:.2f} mN   摆动周期 ≈ {w['period_steps']:.0f} 步" if w["period_steps"] else "", fontsize=8.5, color=INK2)
axs[0].set_ylabel("阻力 [mN]"); axs[0].set_ylim(0, 16.5); axs[0].legend(frameon=False, fontsize=9, loc="lower right", bbox_to_anchor=(1.0, 0.08))
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig_hull_components.png"), dpi=110, facecolor=SURF); plt.close()

# ---------------- 图 2：阻力随速度 + 摩擦系数与经验线对比 ----------------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURF)
ax = axs[0]; ax.set_facecolor(SURF); Ug = np.linspace(0.04, 0.26, 60)
for h, lab in HULLS:
    hh = A["hulls"][h]; us = sorted(float(u) for u in hh["U"]); Fx = [hh["U"][f"{u:.1f}"]["Fx"] for u in us]; sd = [hh["U"][f"{u:.1f}"]["win"]["w20"]["std"] for u in us]
    n = hh["n_total"]; k = Fx[-1] / us[-1] ** n
    ax.plot(Ug, k * Ug ** n, color=COL[h], lw=1.4, ls="--", alpha=0.8)
    ax.errorbar(us, Fx, yerr=sd, fmt="o", color=COL[h], ms=6, capsize=3, lw=1.4, label=f"{lab}   R ∝ U^{n:.2f}")
    ax.annotate(f"{Fx[-1]:.1f}", (us[-1], Fx[-1]), xytext=(6, -4 if h == 'M3_streamlined' else 4), textcoords="offset points", fontsize=9, color=COL[h])
ax.set_xlabel("航速 U [m/s]"); ax.set_ylabel("裸船体阻力 [mN]"); ax.set_title("阻力随航速（点 = CFD ±1σ，虚线 = 两点幂律外推）", loc="left", fontsize=10.5)
ax.set_xlim(0, 0.27); ax.set_ylim(0, 16); ax.grid(color=GRID, lw=0.8); ax.legend(frameon=False, fontsize=9, loc="upper left")
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax2 = ax.secondary_xaxis("top", functions=(lambda u: u / np.sqrt(G * 0.2325), lambda f: f * np.sqrt(G * 0.2325))); ax2.set_xlabel("Froude 数 Fr = U/√(gL)", fontsize=9, color=INK2); ax2.tick_params(labelsize=8, colors=INK2)
ax = axs[1]; ax.set_facecolor(SURF); Re = np.logspace(4.1, 4.9, 60)
ax.plot(Re, 0.075 / (np.log10(Re) - 2) ** 2, color=INK2, lw=1.2, label="ITTC-57 湍流平板 0.075/(lg Re−2)²")
ax.plot(Re, 1.328 / np.sqrt(Re), color=INK2, lw=1.2, ls=":", label="Blasius 层流平板 1.328/√Re")
for h, lab in HULLS:
    hh = A["hulls"][h]; us = sorted(float(u) for u in hh["U"])
    ax.plot([hh["U"][f"{u:.1f}"]["Re"] for u in us], [hh["U"][f"{u:.1f}"]["CF"] for u in us], "o-", color=COL[h], ms=6, lw=1.2, label=f"{lab} C_F（CFD 摩擦 / ½ρU²S）")
    ax.plot([hh["U"][f"{u:.1f}"]["Re"] for u in us], [hh["U"][f"{u:.1f}"]["CP"] for u in us], "s--", color=COL[h], ms=5, lw=1.0, alpha=0.7, mfc="none")
ax.plot([], [], "s--", color=GRAY, mfc="none", label="空心方块：C_P（压差 / ½ρU²S）")
ax.set_xscale("log"); ax.set_xlabel("Reynolds 数 Re = U·L_WL/ν"); ax.set_ylabel("系数 [-]"); ax.set_title("摩擦系数落在层流与湍流平板之间；压差系数是差异所在", loc="left", fontsize=10.5)
ax.set_ylim(0, 0.032); ax.grid(color=GRID, lw=0.8, which="both"); ax.legend(frameon=False, fontsize=8.2, loc="upper right")
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig_hull_speed.png"), dpi=110, facecolor=SURF); plt.close()

# ---------------- 图 3：在整机阻力中的位置 ----------------
fig, ax = plt.subplots(figsize=(11, 4.2), facecolor=SURF); ax.set_facecolor(SURF)
rows = [("整机静止浮态（M1 船体 + 四腿在水下）\n静止浮态报告", R_STATIC_02, INK2, None),
        ("v1.1 报告中船体阻力的估算值（迎流面积比例）", R_HULL_ASSUMED_V11, GRAY, "//")]
for h, lab in HULLS:
    rows.append((f"裸船体 CFD  {lab}", A["hulls"][h]["U"]["0.2"]["Fx"], COL[h], None))
    if h == "M2_roundbottom":
        if eq.get("est"): rows.append(("裸船体 M2 圆滑底  同排水量修正（估算）", eq["est"]["0.2"]["Fx"], BLUE, "//"))
        if eq.get("cfd"): rows.append(("裸船体 M2 圆滑底  同排水量水线 CFD 校核", eq["cfd"]["Fx"], BLUE, "xx"))
ys = np.arange(len(rows))[::-1]
for y, (lab, v, c, hch) in zip(ys, rows):
    ax.barh(y, v, 0.62, color=c, alpha=0.55 if hch else 0.95, hatch=hch, edgecolor=SURF if hch else c)
    ax.text(v + 0.8, y, f"{v:.1f} mN", va="center", fontsize=9.5, color=INK)
ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=9); ax.set_xlabel("U = 0.2 m/s 阻力 [mN]")
ax.set_xlim(0, 90); ax.grid(axis="x", color=GRID, lw=0.8)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
sh = A["hulls"]["M1_current"]["U"]["0.2"]["Fx"] / R_STATIC_02 * 100
ax.set_title(f"裸船体只占四腿在水下时整机阻力的 {sh:.0f} %：船型优化的上限约 3 mN；v1.1 的 25 mN 估算偏高 {R_HULL_ASSUMED_V11/A['hulls']['M1_current']['U']['0.2']['Fx']:.1f} 倍", loc="left", fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig_hull_budget.png"), dpi=110, facecolor=SURF); plt.close()

# ---------------- 图 4：同水线 vs 同排水量的公平性 ----------------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.4), facecolor=SURF)
ax = axs[0]; ax.set_facecolor(SURF); xs = np.arange(3); w = 0.38
disp = [hyd[h]["disp_g"] for h, _ in HULLS]; cfd = [A["hulls"][h]["U"]["0.2"]["disp_g_cfd"] for h, _ in HULLS]
ax.bar(xs - w / 2, disp, w, color=[COL[h] for h, _ in HULLS], label="静水力学（剖面积分）")
ax.bar(xs + w / 2, cfd, w, color=[COL[h] for h, _ in HULLS], alpha=0.45, hatch="//", edgecolor=SURF, label="CFD 浮力反算 F_z/g")
for i, (d1, d2) in enumerate(zip(disp, cfd)): ax.text(i - w / 2, d1 + 4, f"{d1:.0f}", ha="center", fontsize=9); ax.text(i + w / 2, d2 + 4, f"{d2:.0f}", ha="center", fontsize=9, color=INK2)
ax.axhline(hyd["M1_current"]["disp_g"], color=RED, lw=1.0, ls="--"); ax.text(1.0, hyd["M1_current"]["disp_g"] + 6, "机器人重量目标 = M1 排水量 385 g", ha="center", fontsize=8.5, color=RED)
ax.set_xticks(xs); ax.set_xticklabels([lab for _, lab in HULLS], fontsize=9.5); ax.set_ylabel("排水量 [g]  (水线 z = −10 mm)"); ax.set_ylim(0, 530)
ax.set_title("同水线下 M2 少排水 41 g：它其实是一条“更小”的船", loc="left", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.grid(axis="y", color=GRID, lw=0.8)
ax = axs[1]; ax.set_facecolor(SURF)
if eq.get("est"):
    labels = ["M1 现船体", "M2 同水线\n(−10 mm, 344 g)", f"M2 同排水量\n({eq['WL_mm']:.1f} mm, 385 g) 估算", "M3 流线型\n(同水线且同排水量)"]
    vals = [A["hulls"]["M1_current"]["U"]["0.2"]["Fx"], A["hulls"]["M2_roundbottom"]["U"]["0.2"]["Fx"], eq["est"]["0.2"]["Fx"], A["hulls"]["M3_streamlined"]["U"]["0.2"]["Fx"]]
    cols = [INK, BLUE, BLUE, GREEN]; hch = [None, None, "//", None]
    if eq.get("cfd"): labels.insert(3, f"M2 同排水量\nCFD 校核"); vals.insert(3, eq["cfd"]["Fx"]); cols.insert(3, BLUE); hch.insert(3, "xx")
    for i, (v, c, hh_) in enumerate(zip(vals, cols, hch)):
        ax.bar(i, v, 0.6, color=c, alpha=0.5 if hh_ else 0.95, hatch=hh_, edgecolor=SURF if hh_ else c)
        ax.text(i, v + 0.25, f"{v:.1f} mN\n{100*(v/vals[0]-1):+.0f} %" if i else f"{v:.1f} mN", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(labels, fontsize=8.8); ax.set_ylabel("U = 0.2 m/s 阻力 [mN]"); ax.set_ylim(0, 14.5)
    ax.set_title("按同排水量修正后，M2 的优势大半消失；M3 的优势不受影响", loc="left", fontsize=10.5)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8)
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig_hull_fairness.png"), dpi=110, facecolor=SURF); plt.close()

# ---------------- 终端摘要 ----------------
for h, lab in HULLS:
    hh = A["hulls"][h]
    for U in sorted(hh["U"]):
        r = hh["U"][U]; w = r["win"]
        print(f"{lab:<14} U={U}  R={r['Fx']:5.2f} (P {r['Px']:.2f} / V {r['Vx']:.2f})  40%窗 {w['w40']['Fx']:.2f}  ±1σ {w['w20']['std']:.2f}  Re {r['Re']:.0f} Fr {r['Fr']:.3f}  "
              f"C_T {r['CT']:.4f} C_P {r['CP']:.4f} C_F {r['CF']:.4f} (ITTC {r['CF_ittc']:.4f}, Blasius {r['CF_blasius']:.4f}, C_F/ITTC {r['CF_over_ittc']:.2f})  rel M1 {100*r['rel_M1']:+.1f}% (40%: {100*r['rel_M1_w40']:+.1f}%)")
    print(f"   速度指数 n: 总 {hh['n_total']:.2f}  压差 {hh['n_press']:.2f}  摩擦 {hh['n_visc']:.2f}   (40% 窗 总 {hh['n_total_w40']:.2f})")
if eq: print("M2 同排水量:", json.dumps(conv({k: v for k, v in eq.items() if k != 'hydro_eq'}), ensure_ascii=False))
print("→", RES)
