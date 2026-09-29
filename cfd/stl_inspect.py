#!/usr/bin/env python3
"""
STL 体检 + 算例尺寸建议（给 snappyHexMesh 用）
==============================================
用法:
    python3 stl_inspect.py body.stl                       # 一个文件
    python3 stl_inspect.py meshes/*.STL                   # 多个 link 一起看
    python3 stl_inspect.py body.stl --draft 0.03          # 给定吃水(m)，顺便给水线
    python3 stl_inspect.py body.stl --scale 0.001         # 先把 mm 换成 m 再分析

只用 numpy。检查项：单位、包围盒、是否封闭(watertight)、法向是否一致、退化三角形；
然后按 DTC 教程的比例给出计算域、加密盒、locationInMesh 的建议数值。
"""
import sys, os, struct
import numpy as np


def read_stl(path):
    with open(path, "rb") as fh:
        head = fh.read(5)
    if head.lower().startswith(b"solid"):
        # 有些二进制 STL 也以 solid 开头 —— 用文件大小验证
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            fh.seek(80)
            n = struct.unpack("<I", fh.read(4))[0]
        if size == 84 + 50 * n:
            return _read_binary(path)
        return _read_ascii(path)
    return _read_binary(path)


def _read_binary(path):
    with open(path, "rb") as fh:
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        dt = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        rec = np.frombuffer(fh.read(50 * n), dtype=dt, count=n)
    return rec["v"].astype(np.float64)


def _read_ascii(path):
    import re
    txt = open(path, "r", errors="replace").read()
    v = np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt), dtype=float)
    return v.reshape(-1, 3, 3)


def topology(tri):
    """返回 (是否封闭, 是否法向一致, 非流形边数, 边界边数, 退化三角形数)"""
    v = tri.reshape(-1, 3)
    scale = np.abs(v).max() or 1.0
    key = np.round(v / scale * 1e7).astype(np.int64)          # 相对精度 1e-7 合并重合点
    _, idx = np.unique(key, axis=0, return_inverse=True)
    f = idx.reshape(-1, 3)
    bad = (f[:, 0] == f[:, 1]) | (f[:, 1] == f[:, 2]) | (f[:, 0] == f[:, 2])
    degen = int(bad.sum())
    f = f[~bad]
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])          # 有向边
    und = np.sort(e, axis=1)
    _, cnt = np.unique(und, axis=0, return_counts=True)
    boundary = int(np.sum(cnt == 1))
    nonmanifold = int(np.sum(cnt > 2))
    watertight = boundary == 0 and nonmanifold == 0
    # 法向一致：每条有向边只出现一次（相邻三角形沿相反方向遍历共享边）
    _, dcnt = np.unique(e, axis=0, return_counts=True)
    oriented = bool(np.all(dcnt == 1))
    return watertight, oriented, nonmanifold, boundary, degen


def fmt(v):
    return "(" + " ".join(f"{x:.4g}" for x in v) + ")"


def main():
    argv = sys.argv[1:]
    opts, args = {}, []
    i = 0
    while i < len(argv):
        if argv[i].startswith("--"):
            opts[argv[i]] = float(argv[i + 1])
            i += 2
        else:
            args.append(argv[i])
            i += 1
    draft = opts.get("--draft")
    scale = opts.get("--scale", 1.0)
    wl = opts.get("--waterline")
    if not args:
        print(__doc__)
        sys.exit(1)

    allv = []
    print(f"{'文件':<28}{'三角形':>9}  {'封闭':^5} {'法向':^5}  尺寸 (x × y × z) [m]")
    print("-" * 92)
    for p in args:
        tri = read_stl(p) * scale
        wt, ori, nm, bd, dg = topology(tri)
        mn, mx = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
        ext = mx - mn
        flag = lambda b: "OK" if b else "XX"
        print(f"{os.path.basename(p)[:27]:<28}{len(tri):>9}  {flag(wt):^5} {flag(ori):^5}  "
              f"{ext[0]:.4g} × {ext[1]:.4g} × {ext[2]:.4g}")
        if not wt:
            print(f"{'':<28}  !! 不封闭：{bd} 条边界边，{nm} 条非流形边 —— snappyHexMesh 会分不清内外")
        if not ori:
            print(f"{'':<28}  !! 法向不一致 —— 用  surfaceOrient {os.path.basename(p)} '(1e6 1e6 1e6)' fixed.stl  修")
        if dg:
            print(f"{'':<28}  !! {dg} 个退化三角形 —— surfaceClean 可以清掉")
        allv.append(tri.reshape(-1, 3))

    v = np.concatenate(allv)
    mn, mx = v.min(0), v.max(0)
    ext = mx - mn
    L = ext.max()
    print("-" * 92)
    print(f"整体包围盒   min {fmt(mn)}   max {fmt(mx)}")
    print(f"整体尺寸     {ext[0]:.4g} × {ext[1]:.4g} × {ext[2]:.4g} m   （最长边 L = {L:.4g} m）")

    # ---- 单位判断 ----
    if L > 50:
        print(f"\n!! 最长边 {L:.0f} —— 几乎肯定是毫米。加 --scale 0.001 重新跑，"
              f"或在 OpenFOAM 里：surfaceTransformPoints -scale 0.001 in.stl out.stl")
        return
    elif L > 20:
        print(f"\n?? 最长边 {L:.1f} m —— 如果这不是真实尺寸，可能单位是 mm/cm，检查一下")

    # ---- 水线 ----
    if wl is None and draft is not None:
        wl = mn[2] + draft
    if wl is not None:
        sub = (v[:, 2] < wl).mean() * 100
        print(f"\n水线 z = {wl:.4g} m   （约 {sub:.0f}% 的表面点在水下）")

    # ---- 计算域建议（比例取自 DTC 教程：上游 1.6L / 下游 4.2L / 侧 3.1L / 下 2.6L / 上 0.65L）----
    cx, cy = (mn[0] + mx[0]) / 2, (mn[1] + mx[1]) / 2
    z0 = wl if wl is not None else (mn[2] + mx[2]) / 2
    Lx = ext[0] if ext[0] > 0 else L
    dom_min = (mn[0] - 4.2 * Lx, cy - 3.1 * Lx, z0 - 2.6 * Lx)
    dom_max = (mx[0] + 1.6 * Lx, cy + 3.1 * Lx, z0 + 0.65 * Lx)
    print("\n== 建议：blockMeshDict 计算域（假设物体沿 x 放置、来流从 +x 吹向 -x，和 DTC 一样）==")
    print(f"   min {fmt(dom_min)}\n   max {fmt(dom_max)}")
    print("   对称的物体可以只算 y ≤ 中面 的一半（y_max 设为中面），力 ×2。")

    print("\n== 建议：topoSet/refineMesh 三级嵌套加密盒（外→内）==")
    for k, m in enumerate([1.5, 0.6, 0.15], 1):
        bmin = (mn[0] - m * Lx * 2.0, cy - m * Lx, z0 - m * Lx * 0.6)
        bmax = (mx[0] + m * Lx, cy + m * Lx, z0 + m * Lx * 0.6)
        print(f"   c{k}: box {fmt(bmin)} {fmt(bmax)}")
    if wl is not None:
        print(f"   自由液面薄层（最重要）: box {fmt((dom_min[0], dom_min[1], wl - 0.06 * Lx))} "
              f"{fmt((dom_max[0], dom_max[1], wl + 0.04 * Lx))}")

    lim = (mx[0] + 0.37 * Lx, cy + 0.0137 * Lx, z0 - 0.31 * Lx)
    print("\n== 建议：snappyHexMeshDict locationInMesh（在流体里、物体外、避开整数网格线）==")
    print(f"   locationInMesh {fmt(lim)};")

    print("\n== 建议：setFieldsDict 水面 ==")
    if wl is not None:
        print(f"   box (-999 -999 -999) (999 999 {wl:.4g});   // alpha.water = 1")
    else:
        print("   加 --draft <吃水 m> 或 --waterline <z 坐标 m> 我才能给出数值")

    print("\n提醒：所有物体 STL 放进 constant/triSurface/（或 constant/geometry/，看你版本的 snappyHexMeshDict 怎么写），"
          "\n      在 snappyHexMeshDict 的 geometry{} 里每个文件一个条目，refinementSurfaces 里给每个一个 level。")


if __name__ == "__main__":
    main()
