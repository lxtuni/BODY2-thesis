#!/usr/bin/env python3
"""最优脚蹼几何 → 平面形轮廓 CSV + NACA 0021 截面 CSV + 三维 STL（可直接打印 / 进 CFD）。
平面形：lunate（前后缘同时后掠），面积 S 固定，AR、taper、sweep 可调；截面 NACA 0021（t_max 21% @ 0.30c，
Kelly 2023 的 Re 10⁴ 最优 = 海豚尾鳍截面 0.218c @ 0.285c），尾缘闭合。
用法：python3 fluke_design.py [--AR 4] [--taper 0.5] [--sweep 35] [--S 18.5e-4] [--out fluke_AR4]
"""
import argparse, os, math, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--AR", type=float, default=4.0)
ap.add_argument("--taper", type=float, default=0.5)
ap.add_argument("--sweep", type=float, default=35.0)
ap.add_argument("--S", type=float, default=18.5e-4)
ap.add_argument("--tc", type=float, default=0.21)
ap.add_argument("--pivot", type=float, default=0.30, help="俯仰轴（支杆）位置，展中弦的比例")
ap.add_argument("--out", default=None)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
out = a.out or f"fluke_AR{a.AR:g}_t{a.taper:g}_s{a.sweep:g}"


def naca00xx(t, n=60, closed=True):
    """NACA 00xx 半厚度分布（余弦分布的 x），闭合尾缘。"""
    beta = np.linspace(0, math.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    a4 = -0.1036 if closed else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 + a4 * x**4)
    return x, yt


def planform(AR, taper, sweep, S, n=41):
    b = math.sqrt(AR * S)
    y = np.linspace(-b / 2, b / 2, n)
    s = np.abs(2 * y / b)
    shape = 1 - (1 - taper) * s
    c_mid = S / np.trapezoid(shape, y)
    c = c_mid * shape
    x_le = np.abs(y) * math.tan(math.radians(sweep))
    return b, c_mid, y, c, x_le


b, c_mid, y, c, x_le = planform(a.AR, a.taper, a.sweep, a.S)
x_te = x_le + c
# 平面形轮廓（俯视，x 向后为正、y 展向）
outline = np.vstack([np.column_stack([x_le, y]), np.column_stack([x_te[::-1], y[::-1]])])
np.savetxt(os.path.join(HERE, out + "_planform_mm.csv"), outline * 1e3, delimiter=",", header="x_mm,y_mm", comments="", fmt="%.3f")

# 截面
xs, yt = naca00xx(a.tc)
sec = np.vstack([np.column_stack([xs[::-1], yt[::-1]]), np.column_stack([xs[1:], -yt[1:]])])
np.savetxt(os.path.join(HERE, out + "_section_NACA0021.csv"), sec, delimiter=",", header="x_over_c,y_over_c", comments="", fmt="%.5f")

# 三维 STL：逐站放样（每站 = 局部弦缩放的截面，绕俯仰轴不转，只平移），两端翼尖封闭
piv_x = a.pivot * c_mid           # 俯仰轴 x（展中），作为坐标原点
rings = []
for yi, ci, xi in zip(y, c, x_le):
    ring = np.column_stack([xi + sec[:, 0] * ci - piv_x, np.full(len(sec), yi), sec[:, 1] * ci])
    rings.append(ring)
tris = []
for r0, r1 in zip(rings[:-1], rings[1:]):
    n = len(r0)
    for i in range(n - 1):
        tris.append([r0[i], r1[i], r1[i + 1]]); tris.append([r0[i], r1[i + 1], r0[i + 1]])
for ring, flip in ((rings[0], True), (rings[-1], False)):        # 翼尖封盖（扇形三角）
    cen = ring.mean(0)
    for i in range(len(ring) - 1):
        t = [cen, ring[i], ring[i + 1]]
        tris.append(t[::-1] if flip else t)
tri = np.array(tris)
nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
rec = np.zeros(len(tri), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
rec["n"] = nrm; rec["v"] = tri
open(os.path.join(HERE, out + ".stl"), "wb").write(b"\0" * 80 + struct.pack("<I", len(tri)) + rec.tobytes())

vol = float(np.sum(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2]))) / 6)
print(f"平面形 lunate：S = {a.S*1e4:.1f} cm²  AR = {a.AR:g}  展 b = {b*1e3:.1f} mm  展中弦 {c_mid*1e3:.1f} mm  翼尖弦 {c[0]*1e3:.1f} mm  后掠 {a.sweep:g}°")
print(f"截面 NACA 00{int(a.tc*100)}：t_max {a.tc*100:.0f}% @ 0.30c → 展中最大厚度 {a.tc*c_mid*1e3:.1f} mm，前缘半径 {1.1019*a.tc**2*c_mid*1e3:.2f} mm")
print(f"俯仰轴 / 支杆：距展中弦前缘 {a.pivot:g} c = {piv_x*1e3:.1f} mm（STL 原点）")
print(f"STL：{len(tri)} 三角形，体积 {abs(vol)*1e6:.2f} cm³（PLA ≈ {abs(vol)*1e6*1.24:.1f} g）→ {out}.stl / _planform_mm.csv / _section_NACA0021.csv")
