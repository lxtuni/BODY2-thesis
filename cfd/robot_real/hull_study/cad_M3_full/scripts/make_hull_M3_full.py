#!/usr/bin/env python3
"""
BODY2 平滑 M3 全船体：水线以下 = 新的光滑流线型（hull_design.py），水线以上 = 密封身体0826.STEP 原上部（含腿槽、法兰、安装孔、2 mm 壁）
    python3 make_hull_M3_full.py [--out cad_M3_full] [--wall 2.0] [--nst 45] [--npt 33]
输出（零件坐标系同 0826.STEP：x 船长、y 向上、z 横向；单位 mm）：
  BODY2_M3_hull_print.step / .stl   可 3D 打印的整体壳（开口向上，壁厚 --wall），SolidWorks 直接打开
  BODY2_M3_underwater_envelope.step 水线以下实体外形（CFD / 参考用，顶面 = 水线 y = −10）
  BODY2_M3_full_envelope.stl        整船外形（水下新型线 + 原上部外表面）校核用
  fig_M3_full_lines.png             三视图线图（全部光滑曲线）；hull_M3_full_spec.json 主要尺寸
"""
import os, sys, json, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=(os.path.dirname(HERE) if os.path.basename(HERE) == "scripts" else os.path.join(HERE, "cad_M3_full"))); ap.add_argument("--wall", type=float, default=2.0)
ap.add_argument("--nst", type=int, default=45, help="水下放样站位数"); ap.add_argument("--npt", type=int, default=33, help="每半剖面点数"); ap.add_argument("--step", default=os.path.join(HERE, "body0826.step"))
ap.add_argument("--no-fuse", action="store_true", help="不做布尔合并，上下两体分别输出"); ap.add_argument("--loft-tol", type=float, default=2e-3, help="放样曲面逼近公差 mm（越小 STEP 越大；2e-3 已远高于打印精度）"); a = ap.parse_args()
OUT = a.out; os.makedirs(OUT, exist_ok=True); t0 = time.time()
from hull_design import Design, WL, DECK
import cadquery as cq
from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid, BRepBuilderAPI_MakeWire
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.gp import gp_Pnt, gp_Pln, gp_Dir
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCP.TopoDS import TopoDS
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE, TopAbs_SHELL
from OCP.ShapeFix import ShapeFix_Solid
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GCPnts import GCPnts_UniformAbscissa
from scipy.interpolate import UnivariateSpline

D = Design(); S = D.summary()
print(f"设计：LWL {S['LWL']:.1f}  BWL {S['BWL']:.1f}  T3 {S['T3']:.2f}  skew a {S['a']:.3f}  最深 y {S['deepest_y']:.2f} @ x {S['deepest_x']:.1f}  V {S['V_cm3']:.1f} cm³  LCB {S['LCB']:.2f}  Cb {S['Cb']:.3f}")

# ---------------- 原船体：上部（y ≥ −10）与内壁水线轮廓 ----------------
orig = cq.importers.importStep(a.step).solids().vals()[0]
faces = cq.Shape(orig.wrapped).Faces()
bsp = sorted([f for f in faces if f.geomType() == "BSPLINE"], key=lambda f: -f.Area())
outer_surf = BRep_Tool.Surface_s(bsp[0].wrapped); inner_surf = BRep_Tool.Surface_s(bsp[-1].wrapped)
def outline_of(surf, y, ds=0.5):
    face = BRepBuilderAPI_MakeFace(surf, 1e-6).Face(); s = BRepAlgoAPI_Section(face, gp_Pln(gp_Pnt(0, y, 0), gp_Dir(0, 1, 0))); s.Build()
    ex = TopExp_Explorer(s.Shape(), TopAbs_EDGE); pts = []
    while ex.More():
        c = BRepAdaptor_Curve(TopoDS.Edge_s(ex.Current())); ua = GCPnts_UniformAbscissa(c, ds)
        if ua.IsDone(): pts += [(c.Value(ua.Parameter(i)).X(), c.Value(ua.Parameter(i)).Z()) for i in range(1, ua.NbPoints() + 1)]
        ex.Next()
    P = np.array(pts); P = P[P[:, 1] > 1e-6]; P = P[np.argsort(P[:, 0])]; xs, idx = np.unique(np.round(P[:, 0], 4), return_index=True); return xs, P[idx, 1]
xin, bin_ = outline_of(inner_surf, WL); f_in = UnivariateSpline(xin, bin_, k=3, s=0); XIN0, XIN1 = float(xin.min()), float(xin.max())
print(f"原船水线：外 x[{D.xmin:.2f},{D.xmax:.2f}] 半宽 {D.bw(D.xc + 29):.2f}；内 x[{XIN0:.2f},{XIN1:.2f}]（壁厚 ≈ {D.xmin - XIN0:.2f}/{XIN1 - D.xmax:.2f} mm）")

# ---------------- 放样工具 ----------------
def loft_solid(stations, sec_fn, x_end0, x_end1, nm, collar=None, ruled=False):
    """开口剖面（水线→龙骨→水线）放样 + 两端顶点收尖 → 顶面 = 平面 → 缝合成实体
       collar=(hc, dz)：在剖面两端各加一段直线到 y = WL+hc、横向偏移 dz（负 = 向内），顶面在 WL+hc（用于与原上部搭接）"""
    ytop = WL if collar is None else WL + collar[0]
    gen = BRepOffsetAPI_ThruSections(False, ruled, a.loft_tol); gen.CheckCompatibility(True)
    gen.AddVertex(BRepBuilderAPI_MakeVertex(gp_Pnt(x_end0, ytop, 0.0)).Vertex())
    for x in stations:
        z, y = sec_fn(float(x)); x = float(x)
        e = cq.Edge.makeSpline([cq.Vector(x, float(yi), float(zi)) for zi, yi in zip(z, y)])
        if collar is None: w = cq.Wire.assembleEdges([e])
        else:
            hc, dz = collar; sgn = 1.0 if z[-1] > 0 else -1.0
            l0 = cq.Edge.makeLine(cq.Vector(x, ytop, float(z[0]) - sgn * dz), cq.Vector(x, float(y[0]), float(z[0])))      # dz<0：领口向内收
            l1 = cq.Edge.makeLine(cq.Vector(x, float(y[-1]), float(z[-1])), cq.Vector(x, ytop, float(z[-1]) + sgn * dz))
            w = cq.Wire.assembleEdges([l0, e, l1])
        gen.AddWire(w.wrapped)
    gen.AddVertex(BRepBuilderAPI_MakeVertex(gp_Pnt(x_end1, ytop, 0.0)).Vertex()); gen.Build(); shell = gen.Shape()
    fb = ShapeAnalysis_FreeBounds(shell, 1e-4, False, False); closed = fb.GetClosedWires()
    ex = TopExp_Explorer(closed, TopAbs_EDGE); mw = BRepBuilderAPI_MakeWire()
    while ex.More(): mw.Add(TopoDS.Edge_s(ex.Current())); ex.Next()
    top = BRepBuilderAPI_MakeFace(mw.Wire(), True).Face()
    sew = BRepBuilderAPI_Sewing(1e-3); sew.Add(shell); sew.Add(top); sew.Perform(); sewed = sew.SewedShape()
    ex = TopExp_Explorer(sewed, TopAbs_SHELL); sh = TopoDS.Shell_s(ex.Current())
    solid = BRepBuilderAPI_MakeSolid(sh).Solid(); fix = ShapeFix_Solid(solid); fix.Perform(); solid = fix.Solid()
    Sd = cq.Solid(solid)
    if Sd.Volume() < 0: Sd = cq.Solid(TopoDS.Solid_s(solid.Reversed()))
    ok = BRepCheck_Analyzer(Sd.wrapped).IsValid(); bb = Sd.BoundingBox()
    print(f"  {nm}: valid={ok} 体积 {Sd.Volume()/1e3:.1f} cm³  x[{bb.xmin:.1f},{bb.xmax:.1f}] y[{bb.ymin:.2f},{bb.ymax:.2f}] z ±{bb.zmax:.2f}  ({time.time()-t0:.0f} s)")
    return Sd

# ---------------- 水下外形（envelope）与内腔 ----------------
th0 = np.arccos(1 - 1.5 / D.Lh)                                   # 首末站位离端点 1.5 mm
xs = D.xc - D.Lh * np.cos(np.linspace(th0, np.pi - th0, a.nst))
xs = np.sort(np.concatenate([xs, [D.xmin + 0.4, D.xmin + 0.8, D.xmax - 0.8, D.xmax - 0.4]]))   # 两端再各加 2 站：抑制样条放样在尖端顶边的过冲（实测 0.7 mm → 0）
env = loft_solid(xs, lambda x: D.section(x, a.npt), D.xmin, D.xmax, "水下外形")
W = a.wall
def inner_section(x, check=False):
    """外剖面向内偏置成内腔剖面：偏置量 δ = W / |n_yz| （n 为外表面的三维法向，|n_yz| 为其在剖面平面内的分量）
       → 三维法向壁厚恒为 W，即使在首尾龙骨/水线沿 x 很陡的地方也不会变薄（纯 2D 偏置在那里只剩 W·cosα）。δ 上限 2.5 W。
       之后在水线处截断并按弧长重采样；check=True 时检查偏置曲线是否打折（燕尾），打折则抛出 ValueError（该站位弃用）"""
    n4 = 4 * a.npt; h = 0.25
    z, y = D.section(x, n4); zp, yp = D.section(min(x + h, D.xmax - 1e-3), n4); zm, ym = D.section(max(x - h, D.xmin + 1e-3), n4)
    dzdx, dydx = (zp - zm) / (2 * h), (yp - ym) / (2 * h)                 # 同一 t 参数处沿 x 的变化率
    dz, dy = np.gradient(z), np.gradient(y); L = np.hypot(dz, dy); nz, ny = -dy / L, dz / L
    mid = len(z) // 2
    if ny[mid] < 0: nz, ny = -nz, -ny                                    # 指向内（龙骨处向上）
    nx = (dz * dydx - dy * dzdx) / L                                      # 三维法向 n ∝ (n_z, n_y, n_x) = (dy, −dz, dz·dydx − dy·dzdx)/L
    delta = W * np.minimum(np.sqrt(1 + nx ** 2), 2.5)
    zi, yi = z + delta * nz, y + delta * ny
    if check and np.any(np.diff(zi) * np.diff(z) + np.diff(yi) * np.diff(y) <= 0): raise ValueError("offset curve folds")
    def cross(i0, i1):
        f = (WL - yi[i0]) / (yi[i1] - yi[i0]); return zi[i0] + f * (zi[i1] - zi[i0])
    below = yi <= WL - 1e-6; idx = np.where(below)[0]
    if len(idx) < 8: raise ValueError("section too small")
    i_a, i_b = int(idx[0]), int(idx[-1])
    zL = cross(i_a - 1, i_a) if i_a > 0 else zi[0]; zR = cross(i_b + 1, i_b) if i_b < len(zi) - 1 else zi[-1]
    pz = np.concatenate([[zL], zi[i_a:i_b + 1], [zR]]); py = np.concatenate([[WL], yi[i_a:i_b + 1], [WL]])
    if check and (pz.max() - pz.min() < 4.0 or WL - py.min() < 1.5): raise ValueError("section too small")
    sarc = np.concatenate([[0], np.cumsum(np.hypot(np.diff(pz), np.diff(py)))]); su = np.linspace(0, sarc[-1], 2 * a.npt - 1)
    return np.interp(su, sarc, pz), np.interp(su, sarc, py)
def _ok(x):
    try: inner_section(x, check=True); return True
    except ValueError: return False
ok_in = [x for x in xs if _ok(float(x))]
XIN0, XIN1 = D.xmin + 1.5 * W, D.xmax - 1.5 * W                            # 内腔两端顶点：离外形端点 1.5 W（端壁 ≥ W）
print(f"  内腔站位 {len(ok_in)}/{len(xs)}：x[{ok_in[0]:.1f},{ok_in[-1]:.1f}]，顶点 x = {XIN0:.1f}/{XIN1:.1f}", flush=True)
cavity = loft_solid(ok_in, inner_section, XIN0, XIN1, "内腔", ruled=True)   # 内腔用直纹放样：无样条过冲（内表面小平面对实验无影响，壁厚偏差 < 0.05 mm）
shell_lower = env.cut(cavity)
print(f"  水下壳（顶面 = 水线）：valid={shell_lower.isValid()} 材料体积 {shell_lower.Volume()/1e3:.1f} cm³  ({time.time()-t0:.0f} s)", flush=True)
# 搭接环：由 env / cavity 的顶面（y = −10 平面）直接向上拉伸 HC，藏在原壁内（外面：原壁外飘；里面：比原内壁多出 ≤ 0.35 mm 的小台阶）
HC = 1.5
env_top = [f for f in env.Faces() if f.geomType() == "PLANE" and abs(f.Center().y - WL) < 1e-3][0]
cav_top = [f for f in cavity.Faces() if f.geomType() == "PLANE" and abs(f.Center().y - WL) < 1e-3][0]
ring = cq.Solid.extrudeLinear(env_top, cq.Vector(0, HC, 0)).cut(cq.Solid.extrudeLinear(cav_top, cq.Vector(0, HC, 0)))
print(f"  搭接环：valid={ring.isValid()} 体积 {ring.Volume()/1e3:.2f} cm³  ({time.time()-t0:.0f} s)", flush=True)
shell_lower_ring = None
try:
    lw = shell_lower.fuse(ring, glue=True)
    print(f"  水下壳+搭接环：valid={lw.isValid()} solids={len(lw.Solids())} 体积 {lw.Volume()/1e3:.1f}（和 {(shell_lower.Volume()+ring.Volume())/1e3:.1f}）  ({time.time()-t0:.0f} s)", flush=True)
    if lw.isValid() and len(lw.Solids()) == 1: shell_lower_ring = lw.Solids()[0]
except Exception as e: print("  搭接环合并失败：", e, flush=True)

# ---------------- 原上部（y ≥ −10）+ 输出 ----------------
keep = cq.Solid.makeBox(600, 100, 400, cq.Vector(-300, WL, -200))
upper = cq.Shape(orig.wrapped).intersect(keep).Solids()[0]
print(f"  原上部：体积 {upper.Volume()/1e3:.1f} cm³（原整体 {orig.Volume()/1e3:.1f}）  ({time.time()-t0:.0f} s)")
shell_lower = shell_lower.Solids()[0] if len(shell_lower.Solids()) else shell_lower
shell_lower.exportBrep(os.path.join(OUT, "_shell_lower.brep")); upper.exportBrep(os.path.join(OUT, "_upper.brep")); env.exportBrep(os.path.join(OUT, "_env.brep")); ring.exportBrep(os.path.join(OUT, "_ring.brep"))
if shell_lower_ring is not None: shell_lower_ring.exportBrep(os.path.join(OUT, "_shell_lower_ring.brep"))
low_body = shell_lower_ring if shell_lower_ring is not None else shell_lower
two = cq.Compound.makeCompound([upper, low_body])                          # STEP：两体（原上部 + 新水下壳含 1.5 mm 搭接环，重叠在原壁内）→ SolidWorks 组合
three = cq.Compound.makeCompound([upper, shell_lower, ring])               # STL 初版：三壳重叠；最终单壳封闭 STL 由 make_print_stl.py 用网格并集生成
cq.exporters.export(two, os.path.join(OUT, "BODY2_M3_hull_print.step"))
cq.exporters.export(three, os.path.join(OUT, "BODY2_M3_hull_print_3shells.stl"), tolerance=0.02, angularTolerance=0.08)
cq.exporters.export(env, os.path.join(OUT, "BODY2_M3_underwater_envelope.step"))
cq.exporters.export(cq.Compound.makeCompound([env, upper]), os.path.join(OUT, "BODY2_M3_full_envelope.stl"), tolerance=0.03, angularTolerance=0.1)
mass = dict(shell_cm3=shell_lower.Volume() / 1e3, upper_cm3=upper.Volume() / 1e3, total_cm3=(shell_lower.Volume() + upper.Volume()) / 1e3)
mass["PLA_g"] = mass["total_cm3"] * 1.24; mass["PETG_g"] = mass["total_cm3"] * 1.27
mass["ring_cm3"] = ring.Volume() / 1e3
spec = dict(design=S, wall_mm=W, collar_mm=HC, mass=mass, frame="零件坐标同 0826.STEP：x 船长（尖端 −x）、y 向上（水线 −10、顶缘 +15）、z 横向；x_CFD = −x, z_CFD = y",
            original_body=dict(volume_cm3=orig.Volume() / 1e3, displacement_below_WL_cm3=385.1, LCB=11.34, bbox="262 × 116 × 50 (y −35…+15)"))
json.dump(spec, open(os.path.join(OUT, "hull_M3_full_spec.json"), "w"), indent=1, ensure_ascii=False)
print(json.dumps(mass, indent=1)); print("两体 STEP/STL 已输出；单体合并另跑 fuse_bodies.py。done", f"{time.time()-t0:.0f} s")
