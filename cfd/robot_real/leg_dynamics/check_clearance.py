#!/usr/bin/env python3
"""
步态可行性快查：给定步态（环境变量 T SWEEP PSI_START DUTY SWITCH_FRAC PHI_RET …，同 leg_dynamics.py），
算前腿桨板在划水相（张开扇）与回收相（收拢桨）里离船体的最小间隙，以及桨尖是否出水 / 髋角是否越过机构极限。
    T=0.81 DUTY=0.45 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5 python3 check_clearance.py [--hull ../hull_body_only.stl] [--y 0.075] [--out fig.png]
不需要 OpenFOAM；几秒钟。间隙 < 5 mm 视为真机会碰（板厚已计入），整机 overset 算例还需要 ≥ 10–14 mm（会自动裁扇）。
"""
import os, sys, argparse, struct, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--hull", default=os.path.join(HERE, "..", "hull_body_only.stl"))
ap.add_argument("--y", type=float, default=0.075, help="前腿桨中心 |y| [m]（整机算例 0.075；CAD 0.0688）")
ap.add_argument("--dx", type=float, default=0.0, help="腿机构相对前腿的 x 平移（后腿 -0.1435）")
ap.add_argument("--wl", type=float, default=-0.010)
ap.add_argument("--out", default=None, help="输出图 png（默认不画）")
ap.add_argument("--n", type=int, default=240)
a = ap.parse_args()
sp = importlib.util.spec_from_file_location("ld", os.path.join(HERE, "leg_dynamics.py")); ld = importlib.util.module_from_spec(sp); sp.loader.exec_module(ld)
T = ld.T_PER

def read_stl(fn):
    b = open(fn, "rb").read()
    if b[:5] == b"solid" and b"facet" in b[:400]:
        return np.array([[float(x) for x in ln.split()[1:4]] for ln in b.decode(errors="ignore").splitlines() if ln.strip().startswith("vertex")]).reshape(-1, 3, 3)
    n = struct.unpack("<I", b[80:84])[0]
    return np.frombuffer(b[84:84 + n * 50], dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]))["v"].reshape(-1, 3, 3).astype(float)
tree = None
if os.path.exists(a.hull):
    TT = read_stl(a.hull)
    uu, vv = np.meshgrid(np.linspace(0, 1, 4), np.linspace(0, 1, 4)); uu = uu.ravel(); vv = vv.ravel(); k = uu + vv <= 1; uu, vv = uu[k], vv[k]
    pc = (TT[:, None, 0] + uu[None, :, None] * (TT[:, None, 1] - TT[:, None, 0]) + vv[None, :, None] * (TT[:, None, 2] - TT[:, None, 0])).reshape(-1, 3)
    try:
        from scipy.spatial import cKDTree; tree = cKDTree(pc)
    except Exception:
        pc = pc[:: max(1, len(pc) // 20000)]
        class _B:
            def query(self, S):
                d = np.empty(len(S))
                for i in range(0, len(S), 256):
                    blk = S[i:i + 256]; d[i:i + 256] = np.sqrt(((blk[:, None, :] - pc[None]) ** 2).sum(-1)).min(1)
                return d, None
        tree = _B()
else:
    print(f"!! 找不到船体 STL {a.hull}，只做机构/出水检查")

def outline(fan_open):
    S, H, SW, HW = ld.STEM, ld.HEAD, ld.STEM_W, ld.HEAD_W
    if fan_open:
        ang = np.deg2rad(ld.FAN_DEG); th = np.linspace(ang / 2, -ang / 2, 48)
        return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + [(ld.FAN_R * np.sin(x), S + ld.FAN_R * np.cos(x)) for x in th] + [(-SW / 2, S)]
    kk = np.linspace(0, 1, 24); sh = [(HW / 2 * np.sin(np.pi / 2 * q), S + 0.025 * (1 - np.cos(np.pi / 2 * q))) for q in kk]
    right = [(w, z) for w, z in sh if w > SW / 2]
    return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + right + [(HW / 2, S + H), (-HW / 2, S + H)] + [(-w, z) for w, z in right[::-1]] + [(-SW / 2, S)]
def lat_range(poly, s_):
    xs = []; n = len(poly)
    for i in range(n):
        (l1, s1), (l2, s2) = poly[i], poly[(i + 1) % n]
        if (s1 - s_) * (s2 - s_) <= 0 and s1 != s2: xs.append(l1 + (l2 - l1) * (s_ - s1) / (s2 - s1))
    return (min(xs), max(xs)) if xs else None
def sample(poly, thk, n_s=30, n_lat=7):
    ss = [p[1] for p in poly]; pts = []
    for s_ in np.linspace(min(ss) + 1e-4, max(ss) - 1e-4, n_s):
        r = lat_range(poly, s_)
        if r is None: continue
        for la in np.linspace(r[0], r[1], n_lat):
            for off in (-thk / 2, thk / 2): pts.append((s_, la, off))
    return np.array(pts)

tt = np.linspace(0, T, a.n, endpoint=False); tau = tt / T
psi, phi, power = ld.gait(tt); P = [ld.pose(p, q) for p, q in zip(psi, phi)]
E = np.array([p["E"] for p in P]) + np.array([a.dx, 0.0]); u = np.array([p["u"] for p in P]); nn = np.stack([-u[:, 1], u[:, 0]], 1)
locO, locC = sample(outline(True), ld.THK), sample(outline(False), ld.THK)
gap = np.full(len(tt), np.nan); ztop = np.zeros(len(tt)); ztip = np.zeros(len(tt))
for k in range(len(tt)):
    loc = locO if power[k] else locC
    S = np.column_stack([E[k, 0] + loc[:, 0] * u[k, 0] + loc[:, 2] * nn[k, 0], a.y + loc[:, 1], E[k, 1] + loc[:, 0] * u[k, 1] + loc[:, 2] * nn[k, 1]])
    ztop[k] = S[:, 2].max(); ztip[k] = (E[k] + (ld.STEM + ld.HEAD) * u[k])[1]
    if tree is not None: gap[k] = tree.query(S)[0].min()
pw = power; rc = ~power
psi_deg = np.degrees(psi)
print(f"步态：T={T} s  DUTY={ld.DUTY}  SWEEP={ld.SWEEP}°  PSI_START={np.degrees(ld.PSI0):.0f}° → ψ 范围 {psi_deg.min():.1f}–{psi_deg.max():.1f}°  SWITCH_FRAC={ld.SWITCH_FRAC}  PHI_RET={ld.PHI_RET}°   桨中心 y={a.y} m")
lim = getattr(ld, "PSI_LIMITS", None)
if psi_deg.min() < 70 - 1e-6 or psi_deg.max() > 150 + 1e-6: print(f"!! 髋角超出机构范围 70–150°")
if tree is not None:
    i1 = np.argmin(np.where(pw, gap, np.inf)); i2 = np.argmin(np.where(rc, gap, np.inf))
    print(f"划水相（张开扇）与船体最小间隙 {gap[i1]*1000:5.1f} mm  @ τ={tau[i1]:.2f}, ψ={psi_deg[i1]:.0f}°")
    print(f"回收相（收拢桨）与船体最小间隙 {gap[i2]*1000:5.1f} mm  @ τ={tau[i2]:.2f}, ψ={psi_deg[i2]:.0f}°")
    ok = gap[i1] > 0.005 and gap[i2] > 0.005
    print("结论：" + ("真机不碰船（>5 mm）" if ok else "!! 间隙 <5 mm，真机可能碰船") + ("；整机 overset 算例需 ≥10–14 mm，张开扇内侧会被自动裁掉" if min(gap[i1], gap[i2]) < 0.014 else ""))
print(f"桨顶最高 z = {ztop.max()*1000:.1f} mm（水线 {a.wl*1000:.0f} mm）{'  !! 桨在回收/划水时出水' if ztop.max() > a.wl else '  桨全程在水下'}；桨尖最低 z = {ztip.min()*1000:.1f} mm")
if a.out:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 3.4)); ax.plot(tau, gap * 1e3, color="#0b0b0b"); ax.axvspan(0, ld.DUTY, color="#2a78d6", alpha=0.08)
    ax.axhline(5, color="#e34948", ls=":", lw=1); ax.axhline(14, color="#a8a7a1", ls=":", lw=1); ax.set_xlabel("tau"); ax.set_ylabel("gap to hull [mm]"); ax.set_ylim(0, None)
    ax.set_title(f"SWEEP {ld.SWEEP} PSI_START {np.degrees(ld.PSI0):.0f} DUTY {ld.DUTY} T {T}: power(open) / recovery(closed)", fontsize=9)
    plt.tight_layout(); plt.savefig(a.out, dpi=110); print("→", a.out)
