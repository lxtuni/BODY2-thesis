#!/usr/bin/env python3
"""由 cad_M3_full/_upper.brep + _shell_lower_ring.brep 逐体细分（各自封闭），再用 manifold3d 做网格并集 → 单壳封闭的打印 STL
   输出 cad_M3_full/BODY2_M3_hull_print.stl（单壳，并集成功时）/ BODY2_M3_hull_print_2shells.stl（两壳重叠，备用）"""
import os, sys, time, numpy as np, cadquery as cq, trimesh
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "cad_M3_full") if os.path.isdir(os.path.join(HERE, "cad_M3_full")) else (os.path.dirname(HERE) if os.path.basename(HERE) == "scripts" else HERE); import tempfile; TMP = os.path.join(tempfile.gettempdir(), "m3_stl"); os.makedirs(TMP, exist_ok=True)
tol, atol = 0.015, 0.06
def mesh_of(name):
    s = cq.Shape.importBrep(os.path.join(OUT, name + ".brep")); fn = os.path.join(TMP, name + ".stl")
    cq.exporters.export(s, fn, tolerance=tol, angularTolerance=atol)
    m = trimesh.load(fn, force="mesh"); parts = sorted(m.split(only_watertight=False), key=lambda p: -len(p.faces)); big = parts[0]
    big.fix_normals(); print(f"  {name:18s} faces {len(big.faces):6d} watertight={big.is_watertight} vol {big.volume/1e3:.2f} cm³ (B-rep {s.Volume()/1e3:.2f}); 丢弃碎片 {len(parts)-1} 个/{sum(len(p.faces) for p in parts[1:])} 面", flush=True)
    return big
t0 = time.time()
up = mesh_of("_upper"); low = mesh_of("_shell_lower_ring")
two = trimesh.util.concatenate([up, low]); two.export(os.path.join(OUT, "BODY2_M3_hull_print_2shells.stl"))
try:
    u = trimesh.boolean.union([up, low], engine="manifold")
    parts = u.split(only_watertight=False)
    print(f"  并集：faces {len(u.faces)} watertight={u.is_watertight} parts={len(parts)} vol {u.volume/1e3:.2f} cm³ (两体和 {(up.volume+low.volume)/1e3:.2f}，重叠 ≈ 搭接环)  {time.time()-t0:.0f}s", flush=True)
    if u.is_watertight and len(parts) == 1:
        u.export(os.path.join(OUT, "BODY2_M3_hull_print.stl")); print("  → BODY2_M3_hull_print.stl（单壳封闭）")
    else: print("  并集不封闭，保留两壳版本")
except Exception as e: print("  并集失败：", e)
