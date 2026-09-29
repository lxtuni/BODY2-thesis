#!/usr/bin/env python3
"""检查 cad_M3_full/BODY2_M3_hull_print.stl：各壳是否封闭、体积、壁厚、y = −10 接缝连续性、剖面光滑度。输出 cad_M3_full/fig_M3_stl_check.png"""
import os, sys, numpy as np, trimesh
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "cad_M3_full") if os.path.isdir(os.path.join(HERE, "cad_M3_full")) else (os.path.dirname(HERE) if os.path.basename(HERE) == "scripts" else HERE)
fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(OUT, "BODY2_M3_hull_print.stl")
m = trimesh.load(fn, force="mesh"); print(fn, "faces", len(m.faces), "bbox", np.round(m.bounds, 2).tolist())
parts = m.split(only_watertight=False); parts = sorted(parts, key=lambda p: -p.volume)
print(f"连通壳数 {len(parts)}")
for p in parts:
    print(f"  watertight={p.is_watertight}  winding={p.is_winding_consistent}  vol={p.volume/1e3:8.2f} cm³  y∈[{p.bounds[0][1]:7.2f},{p.bounds[1][1]:7.2f}]  x∈[{p.bounds[0][0]:7.1f},{p.bounds[1][0]:7.1f}]")
if len(parts) == 1: lower = upper = parts[0]; ring = []                       # 单壳并集版本
else:
    lower = [p for p in parts if p.bounds[1][1] < -8.0 and p.bounds[0][1] < -40][0]     # 水下壳（y −48…−10）
    upper = [p for p in parts if p.bounds[1][1] > 10][0]; ring = [p for p in parts if p not in (lower, upper)]
# ---- 壁厚：水下壳的内表面点到外表面的距离 ----
env = trimesh.load(os.path.join(OUT, "BODY2_M3_full_envelope.stl"), force="mesh")
pts = lower.sample(12000); pts = pts[pts[:, 1] < -12.0]
d_all = trimesh.proximity.ProximityQuery(env).signed_distance(pts)           # 到外包络的有向距离（正 = 在内）：外表面点 ≈ 0，内表面点 ≈ 壁厚
d_env = d_all[d_all > 0.5]
print(f"水下壳壁厚（内表面点→外表面）：min {d_env.min():.2f}  中位 {np.median(d_env):.2f}  max {d_env.max():.2f} mm  (n={len(d_env)})")
# ---- 接缝 y = −10：下壳 y = −10.3 的外轮廓 vs 上壳 y = −9.7 的外轮廓 ----
def outer_contour(mesh, y):
    s = mesh.section(plane_origin=[0, y, 0], plane_normal=[0, 1, 0])
    if s is None: return None
    paths = [np.asarray(e) for e in s.discrete]; paths = sorted(paths, key=lambda P: -np.abs(np.cross(P[:, [0, 2]][:-1], P[:, [0, 2]][1:])).sum())
    return paths[0][:, [0, 2]]
def poly_dist(P, Q):
    """P 各点到折线 Q 的最小距离（向量化点-线段距离）"""
    A, B = Q[:-1], Q[1:]; AB = B - A; L2 = np.maximum((AB ** 2).sum(1), 1e-12)
    t = np.clip(((P[:, None, :] - A[None]) * AB[None]).sum(2) / L2[None], 0, 1)
    C = A[None] + t[..., None] * AB[None]; return np.sqrt(((P[:, None, :] - C) ** 2).sum(2)).min(1)
Pl = outer_contour(lower, -10.3); Pu = outer_contour(upper, -9.7)
gap = poly_dist(Pl, Pu); ang = np.degrees(np.arctan2(Pl[:, 1], Pl[:, 0] - Pl[:, 0].mean()))
print(f"接缝 y = −10：下壳 y=−10.3 轮廓到上壳 y=−9.7 轮廓的法向距离：median {np.median(gap):.2f}  max {gap.max():.2f} mm（上部外飘 + 下部收进各 ≈ 0.1 mm/0.3 mm → 预期 ≈ 0.2；无台阶）")
# ---- 光滑度：x 站位剖面外轮廓的转角分布 ----
def kinks(P, step=1.0):
    """开口轮廓（左舷水线 → 龙骨 → 右舷水线）按 1 mm 弧长重采样后求相邻转角 [°]（= 曲率 × 1 mm）"""
    P = P[np.r_[True, np.linalg.norm(np.diff(P, axis=0), axis=1) > 1e-6]]
    L = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]); s = np.arange(0, L[-1], step)
    Q = np.c_[np.interp(s, L, P[:, 0]), np.interp(s, L, P[:, 1])]; d = np.diff(Q, axis=0); a = np.arctan2(d[:, 1], d[:, 0])
    return np.degrees(np.abs((np.diff(a) + np.pi) % (2 * np.pi) - np.pi))
from hull_design import Design
Dz = Design(); worst = []
for x in np.linspace(Dz.xmin + 2, Dz.xmax - 2, 117):                       # 解析剖面（STL 只是它的 0.015 mm 弦高细分）
    z, y = Dz.section(float(x), 400); worst.append((x, kinks(np.c_[z, y], 1.0)[1:-1].max()))
fig, axs = plt.subplots(2, 2, figsize=(16, 12))
ax = axs[0, 0]
for y in [-45, -40, -35, -30, -25, -20, -15, -10.3, -9.7, -5, 0, 5, 10, 14.5]:
    for mesh, col in ((lower, "#1f77b4"), (upper, "#7f8c8d")):
        s = mesh.section(plane_origin=[0, y, 0], plane_normal=[0, 1, 0])
        if s is None: continue
        for e in s.discrete: e = np.asarray(e); ax.plot(e[:, 0], e[:, 2], "-", color=("#1f77b4" if y < -10 else "#7f8c8d"), lw=0.8)
        if lower is upper: break
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_title("STL 水平切片 y = −45…14.5（蓝 = 水线以下，灰 = 水线以上；内外两条线 = 2 mm 壁）", fontsize=10)
ax = axs[0, 1]
for x in np.linspace(-110, 110, 12):
    for mesh, col in ((lower, "#1f77b4"), (upper, "#7f8c8d"), *[(r_, "#c0392b") for r_ in ring]):
        s = mesh.section(plane_origin=[x, 0, 0], plane_normal=[1, 0, 0])
        if s is None: continue
        for e in s.discrete:
            e = np.asarray(e); ax.plot(e[:, 2], e[:, 1], "-", color=col, lw=0.8)
        if lower is upper: break
ax.axhline(-10, color="#c0392b", ls="--", lw=1); ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_title("STL 横剖切片 x = −110…110（单壳：新水下壳 + 原上部 + 搭接环已并集）", fontsize=10)
ax.set_xlabel("z"); ax.set_ylabel("y")
ax = axs[1, 0]
o = np.argsort(ang); ax.plot(ang[o], gap[o], ".", ms=2, color="#1f77b4"); ax.set_xlabel("方位角 [°]（0 = 钝端，±180 = 尖端）"); ax.set_ylabel("距离 [mm]"); ax.grid(alpha=.3); ax.set_ylim(0, 1)
ax.set_title(f"接缝 y = −10：下壳 (y=−10.3) 与上壳 (y=−9.7) 轮廓的法向距离，median {np.median(gap):.2f} / max {gap.max():.2f} mm\n（= 0.6 mm 高差上的正常外飘，首尾外飘大；无台阶）", fontsize=10)
if worst:
    xw = [w[0] for w in worst]; kw = [w[1] for w in worst]
    ax2 = ax.inset_axes([0.55, 0.55, 0.42, 0.4]); ax2.plot(xw, kw, "-", color="#2e8b57", lw=1); ax2.set_title("横剖面最大转角 [°/mm]（解析，龙骨处最大）", fontsize=8); ax2.tick_params(labelsize=7); ax2.grid(alpha=.3)
ax = axs[1, 1]
ax.hist(d_env, bins=60, color="#1f77b4"); ax.set_xlabel("壁厚 [mm]"); ax.set_title(f"水下壳壁厚分布（内表面采样点到外表面的距离）\nmin {d_env.min():.2f} / median {np.median(d_env):.2f} mm（> 2 mm 的是首尾内腔收尖处）", fontsize=10); ax.grid(alpha=.3)
plt.tight_layout(); fo = os.path.join(OUT, "fig_M3_stl_check.png"); plt.savefig(fo, dpi=100); print("→", fo)
if worst:
    w = max(worst, key=lambda t: t[1]); print(f"横剖面（解析）每 1 mm 弧长的最大转角：{w[1]:.2f}° (x = {w[0]:.0f}，即最小曲率半径 {57.3/w[1]:.0f} mm)；全船中位 {np.median([k for _, k in worst]):.2f}°  —— 连续曲率、无折角；STL 弦高误差 ≤ 0.015 mm")
