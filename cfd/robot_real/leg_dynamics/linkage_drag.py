#!/usr/bin/env python3
"""
连杆（曲柄1、链接杆、摇臂、曲柄2、从动杆）的水动力估算：CFD 里只有脚蹼受力，这里用条带法把整套连杆
在水线以下部分的阻力算出来，和脚蹼 CFD 力放在一起比较，并给出水线位置的敏感性。
    python3 linkage_drag.py [results_dir] [--wl -0.010] [--cd-scale 1.0] [--out results_dir]
输出：results/fig_linkage.png, results/linkage.json；控制台里的逐杆表。
杆件截面从 STEP 零件 STL 量得（横向厚 t_y × 面内宽 w）；杆在面内垂直于杆轴运动时迎流面 = 湿长 × t_y，
顺流深度 = w，取矩形柱 C_d（近方形 2.0，w/t≈1.7 取 1.8，w/t≈2.7 取 1.4）。
"""
import os, sys, json, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

# 零件实测（parts/*.stl，FL 腿）：销距 [m]、零件全长 [m]、横向厚 t_y [m]、面内宽 w [m]、C_d
LINKS = {
    "crank1":   dict(label="曲柄1 O1–A",  pins=0.028,  part=0.037, ty=0.0053, w=0.0093, cd=2.0),
    "rod":      dict(label="链接杆 A–B",  pins=0.035,  part=0.040, ty=0.0030, w=0.0050, cd=1.8),
    "rocker":   dict(label="摇臂 B–O2–C", pins=0.0448, part=0.050, ty=0.0030, w=0.0080, cd=1.4),
    "crank2":   dict(label="曲柄2 O2–D",  pins=0.039,  part=0.047, ty=0.0078, w=0.0090, cd=2.0),
    "follower": dict(label="从动杆 C–E",  pins=0.039,  part=0.045, ty=0.0030, w=0.0050, cd=1.8),
}
ORDER = ["crank2", "follower", "rocker", "rod", "crank1"]

def joints(ld, t):
    psi, phi, power = ld.gait(t)
    J = {k: np.zeros((len(t), 2)) for k in ("A", "B", "C", "D", "E")}
    for i, (p, q) in enumerate(zip(psi, phi)):
        P = ld.pose(float(p), float(q))
        th1 = ld.crank1_from_rocker(float(p) + np.pi)
        J["A"][i] = ld.O1 + ld.L_CRANK1 * ld.rot(th1)
        for k in ("B", "C", "D", "E"): J[k][i] = P[k]
    J["O1"] = np.tile(ld.O1, (len(t), 1)); J["O2"] = np.tile(ld.O2, (len(t), 1))
    return J, psi, phi, power

SEGS = {"crank1": ("O1", "A"), "rod": ("A", "B"), "rocker": ("B", "C"), "crank2": ("O2", "D"), "follower": ("C", "E")}

def link_forces(ld, t, J, wl, cd_scale=1.0, nstrip=60, rho=998.8, U=0.0):
    """返回 {link: dict(Fx(t), Fz(t), wet(t))}，力为流体作用在杆上的力 [N]；U = 机器人前进速度（机体系里水以 -U 流过）"""
    out = {}
    for name, (a, b) in SEGS.items():
        L = LINKS[name]; P, Q = J[a], J[b]
        e = Q - P; Lp = np.linalg.norm(e, axis=1, keepdims=True); e = e / Lp
        ext = (L["part"] - L["pins"]) / 2                      # 零件比销距长的部分，两端各伸一半
        s = np.linspace(-ext, L["pins"] + ext, nstrip); ds = s[1] - s[0]
        X = P[:, None, :] + s[None, :, None] * e[:, None, :]  # (nt, ns, 2)
        V = np.gradient(X, t, axis=0); V = V + np.array([U, 0.0])[None, None, :]   # 相对水的速度
        n = np.stack([-e[:, 1], e[:, 0]], axis=1)             # 面内法向
        vn = np.einsum("tsk,tk->ts", V, n)
        wet = (X[:, :, 1] < wl).astype(float)
        dF = -0.5 * rho * L["cd"] * cd_scale * L["ty"] * ds * np.abs(vn) * vn * wet   # 沿 n 的力（阻力，与速度反向）
        Fx = (dF * n[:, None, 0]).sum(1); Fz = (dF * n[:, None, 1]).sum(1)
        out[name] = dict(Fx=Fx, Fz=Fz, wet=wet.sum(1) * ds, vmax=np.sqrt((V**2).sum(-1))[wet > 0].max() if wet.any() else 0.0)
    return out

def leg_forces(ld, t, phase=0.0, dx=0.0, wl=-0.010, U=0.0, cd_scale=1.0):
    """整机用：某条腿（相位 phase、机构 x 平移 dx）在时刻数组 t 上整套连杆的水下受力合计 (Fx, Fz) [N] 及各杆明细"""
    T = ld.T_PER
    J, psi, phi, power = joints(ld, t + phase * T)
    J = {k: v + np.array([dx, 0.0]) for k, v in J.items()}
    lf = link_forces(ld, t, J, wl, cd_scale, U=U)
    Fx = sum(v["Fx"] for v in lf.values()); Fz = sum(v["Fz"] for v in lf.values())
    return Fx, Fz, lf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default=os.path.join(HERE, "results"))
    ap.add_argument("--wl", type=float, default=-0.010, help="水线 z [m]（真机约在船体中部 −0.010；船底 −0.035）")
    ap.add_argument("--cd-scale", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    R = os.path.abspath(os.path.expanduser(a.results)); OUT = os.path.abspath(a.out) if a.out else R
    os.makedirs(OUT, exist_ok=True)

    # 步态参数跟 CFD 算例一致
    gj = None
    for tag in ("open", "closed"):
        f = os.path.join(R, f"gait_{tag}.json")
        if os.path.exists(f): gj = json.load(open(f)); break
    if gj:
        for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W"]:
            if k in gj: os.environ[k] = str(gj[k])
        if "psi0_deg" in gj: os.environ["PSI_START"] = str(gj["psi0_deg"])
    import importlib.util
    spec = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(spec); spec.loader.exec_module(ld)
    T = ld.T_PER; DUTY = ld.DUTY
    t = np.linspace(0, T, 2001); tau = t / T
    J, psi, phi, power = joints(ld, t)

    lf = link_forces(ld, t, J, a.wl, a.cd_scale)
    tot_Fx = sum(v["Fx"] for v in lf.values()); tot_Fz = sum(v["Fz"] for v in lf.values())
    tr = getattr(ld, "SWITCH_FRAC", 0.3) * (1 - DUTY)
    phases = [("划水相", 0, DUTY), ("收拢切换", DUTY, DUTY + tr), ("蜷缩前扫", DUTY + tr, 1 - tr), ("伸展切换", 1 - tr, 1.0)]
    def contrib(F):
        return [float(np.trapezoid(F[(tau >= p) & (tau < q)], t[(tau >= p) & (tau < q)]) / T * 1e3) for _, p, q in phases]

    # CFD 脚蹼力（最后一个周期）作对比
    cfd = None
    for tag in ("open", "closed"):
        f = os.path.join(R, f"force_{tag}.dat")
        if os.path.exists(f):
            d = np.loadtxt(f, comments="#"); d = d[np.argsort(d[:, 0])]
            n = int(np.floor(d[-1, 0] / T + 1e-6)); k = (d[:, 0] >= (n - 1) * T) & (d[:, 0] <= n * T)
            Fx = np.interp(t + (n - 1) * T, d[k, 0], d[k, 1])
            cfd = dict(tag=tag, Fx=Fx, mean=float(Fx.mean() * 1e3), peak=float(np.abs(Fx).max() * 1e3), contrib=contrib(Fx)); break

    data = dict(WL=a.wl, cd_scale=a.cd_scale, T=T, DUTY=DUTY, links={}, phases=[p[0] for p in phases])
    print(f"水线 z = {a.wl:+.3f} m，C_d×{a.cd_scale}，周期 T = {T} s\n{'杆件':<14}{'湿长 min–max [mm]':>20}{'湿部最大速度 [m/s]':>20}{'F_x 均值 [mN]':>14}{'|F| 峰值 [mN]':>14}{'F_z 均值 [mN]':>14}")
    for name in ORDER:
        v = lf[name]; L = LINKS[name]
        d = dict(label=L["label"], wet_mm=[float(v["wet"].min() * 1e3), float(v["wet"].max() * 1e3)], vmax=float(v["vmax"]),
                 mean_Fx=float(v["Fx"].mean() * 1e3), peak=float(np.hypot(v["Fx"], v["Fz"]).max() * 1e3), mean_Fz=float(v["Fz"].mean() * 1e3), contrib=contrib(v["Fx"]))
        data["links"][name] = d
        print(f"{L['label']:<14}{d['wet_mm'][0]:9.1f} – {d['wet_mm'][1]:<8.1f}{d['vmax']:>20.2f}{d['mean_Fx']:>14.2f}{d['peak']:>14.1f}{d['mean_Fz']:>14.2f}")
    data["total"] = dict(mean_Fx=float(tot_Fx.mean() * 1e3), peak=float(np.hypot(tot_Fx, tot_Fz).max() * 1e3), peak_Fx=float(np.abs(tot_Fx).max() * 1e3),
                         mean_Fz=float(tot_Fz.mean() * 1e3), contrib=contrib(tot_Fx))
    print(f"{'合计':<14}{'':>20}{'':>20}{data['total']['mean_Fx']:>14.2f}{data['total']['peak']:>14.1f}{data['total']['mean_Fz']:>14.2f}")
    if cfd:
        data["cfd"] = dict(tag=cfd["tag"], mean_Fx=cfd["mean"], peak=cfd["peak"], contrib=cfd["contrib"])
        data["ratio_mean"] = abs(data["total"]["mean_Fx"]) / abs(cfd["mean"]) if cfd["mean"] else None
        data["ratio_peak"] = data["total"]["peak"] / cfd["peak"]
        print(f"脚蹼 CFD（{cfd['tag']}）：均值 {cfd['mean']:+.1f} mN，峰值 {cfd['peak']:.0f} mN → 连杆合计占均值 {100*data['ratio_mean']:.1f}%、占峰值 {100*data['ratio_peak']:.1f}%")

    # 水线敏感性
    sens = []
    for wl in (-0.015, -0.010, -0.005, 0.0, 0.005, 0.010):
        for cs in (1.0, 1.5):
            r = link_forces(ld, t, J, wl, cs); Fx = sum(v["Fx"] for v in r.values()); Fz = sum(v["Fz"] for v in r.values())
            sens.append(dict(WL=wl, cd_scale=cs, mean_Fx=float(Fx.mean() * 1e3), peak=float(np.hypot(Fx, Fz).max() * 1e3),
                             wet_total_max=float(sum(v["wet"].max() for v in r.values()) * 1e3)))
    data["sensitivity"] = sens
    print("水线敏感性（合计 F_x 均值 / 峰值 [mN]）：" + "  ".join(f"WL{s['WL']:+.3f}: {s['mean_Fx']:+.2f}/{s['peak']:.1f}" for s in sens if s["cd_scale"] == 1.0))
    json.dump(data, open(os.path.join(OUT, "linkage.json"), "w"), indent=1, ensure_ascii=False)

    # ---- 图
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    import glob as g
    for f in [x for x in font_manager.findSystemFonts() if "NotoSansCJK" in x and "Regular" in x] + g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + g.glob("/mnt/c/Windows/Fonts/simhei.ttf"):
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
        except Exception: pass
    plt.rcParams["axes.unicode_minus"] = False
    SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    COL = {"crank2": "#2a78d6", "follower": "#eb6834", "rocker": "#3a9c6e", "rod": "#8e6bd6", "crank1": "#a8a7a1"}
    fig, axs = plt.subplots(2, 2, figsize=(15, 9.5), facecolor=SURF)
    for ax in axs.flat:
        ax.set_facecolor(SURF); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.grid(axis="y", color=GRID, lw=0.8)
    ax = axs[0, 0]
    for name in ORDER: ax.plot(tau, lf[name]["wet"] * 1e3, color=COL[name], lw=1.8, label=LINKS[name]["label"])
    for _, p, q in phases: ax.axvline(p, color=INK2, lw=0.6, alpha=0.35)
    ax.axvspan(0, DUTY, color="#eb6834", alpha=0.05); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.set_ylabel("水线以下长度 [mm]")
    ax.set_title(f"各杆没入水中的长度（水线 z = {a.wl:+.3f} m）", loc="left", fontsize=11); ax.legend(frameon=False, fontsize=8.5)
    ax = axs[0, 1]
    for name in ORDER: ax.plot(tau, lf[name]["Fx"] * 1e3, color=COL[name], lw=1.4, label=LINKS[name]["label"])
    ax.plot(tau, tot_Fx * 1e3, color=INK, lw=2.4, label="连杆合计")
    for _, p, q in phases: ax.axvline(p, color=INK2, lw=0.6, alpha=0.35)
    ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_xlabel("t / T"); ax.set_ylabel("F_x [mN]  (+ 推力)")
    ax.set_title("连杆 x 向力（条带法）", loc="left", fontsize=11); ax.legend(frameon=False, fontsize=8.5)
    ax = axs[1, 0]
    labels = ["脚蹼 CFD" if cfd else "脚蹼（无 CFD）", "连杆合计"] + [LINKS[n]["label"] for n in ORDER]
    means = [cfd["mean"] if cfd else 0.0, data["total"]["mean_Fx"]] + [data["links"][n]["mean_Fx"] for n in ORDER]
    peaks = [cfd["peak"] if cfd else 0.0, data["total"]["peak"]] + [data["links"][n]["peak"] for n in ORDER]
    y = np.arange(len(labels))[::-1]
    ax.barh(y + 0.2, means, height=0.38, color="#2a78d6", label="周期平均 F_x")
    ax.barh(y - 0.2, peaks, height=0.38, color="#a8a7a1", label="|F| 峰值")
    for yi, m, pk in zip(y, means, peaks):
        ax.text(max(m, 0) + 1.5, yi + 0.2, f"{m:+.2f}", va="center", fontsize=8.5, color=INK); ax.text(pk + 1.5, yi - 0.2, f"{pk:.1f}", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5); ax.set_xlabel("[mN]"); ax.set_xscale("symlog", linthresh=5); ax.grid(axis="x", color=GRID, lw=0.8); ax.grid(axis="y", visible=False)
    ax.axvline(0, color=INK, lw=0.8); ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title("量级对比：脚蹼 CFD vs 连杆（symlog 轴）", loc="left", fontsize=11)
    ax = axs[1, 1]
    for cs, col, lab in ((1.0, "#2a78d6", "C_d 实测截面"), (1.5, "#eb6834", "C_d × 1.5")):
        ss = [s for s in sens if s["cd_scale"] == cs]
        ax.plot([s["WL"] * 1e3 for s in ss], [s["peak"] for s in ss], "o-", color=col, lw=1.8, label=f"峰值 {lab}")
        ax.plot([s["WL"] * 1e3 for s in ss], [abs(s["mean_Fx"]) for s in ss], "s--", color=col, lw=1.4, label=f"|均值| {lab}")
    ax.axvline(a.wl * 1e3, color=INK2, lw=0.8, ls=":"); ax.text(a.wl * 1e3, ax.get_ylim()[1] * 0.95, " 本报告水线", fontsize=8.5, color=INK2)
    ax.set_xlabel("水线 z [mm]（船底 −35，船体中部 −10）"); ax.set_ylabel("连杆合计 [mN]"); ax.legend(frameon=False, fontsize=8.5)
    ax.set_title("连杆力对水线位置的敏感性", loc="left", fontsize=11)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_linkage.png"), dpi=110, facecolor=SURF); plt.close()
    print("图:", os.path.join(OUT, "fig_linkage.png"), " 数据:", os.path.join(OUT, "linkage.json"))

if __name__ == "__main__":
    main()
