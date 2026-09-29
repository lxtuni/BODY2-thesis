#!/usr/bin/env python3
"""月牙尾鳍脚：直接替换原“扇形足”的三件套（cadquery）
 零件坐标系 = 原扇形足零件坐标系：x 沿安装片/支杆（装机后 ψ = 90° 时竖直向下），y 弦向（向后），z 展向（横向）
 件 1 strut  ：原安装片（两个 Ø2.7 孔，间距 15，原样复制）+ 流线支杆 + 铰座块（含 ±40° 限位爪）      PETG/PLA
 件 2 fluke  ：月牙尾鳍 AR 4、S 18.5 cm²、NACA 0021 + 根部夹块                                          PETG/PLA
 件 3 flexure：TPU 95A 柔性铰片（两端 2 mm 厚夹持片 + 中间薄片）                                       TPU
 装配：flexure 前片插入 strut 铰座槽、后片插入 fluke 根部槽，各 2 颗 M2×6 自攻螺钉；A/B 两个镜像版本（头朝 +y / −y）
"""
import os, math, json, numpy as np
import cadquery as cq

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(HERE, "fan_foot_orig.step")

# ---------------- 参数（mm）
X_E = -12.21                 # 原零件 E 孔（平行四边形末端）x 坐标
TAB_END = 8.0                # 保留原安装片到 x = 8
Z_C = 1.5                    # 原安装片厚度中心 z（0…3）
STRUT_A, STRUT_B = 6.0, 3.0  # 支杆椭圆半轴：y（顺流）6，z（横向）3
L_HINGE = 115.0              # E → 铰轴 距离（沿 x）
X_H = X_E + L_HINGE          # 铰轴 x（= 柔性片中面）
BLK_H = 8.0                  # 铰座块半高（x 方向）
BLK_W = 11.0                 # 铰座块半宽（z 方向）
GAP = 6.0                    # 柔性片自由段长度（y）
TAB_L = 8.0                  # 柔性片夹持片长度（y）
T_FLEX = 1.0                 # 柔性片自由段厚度（x）——目标 k ≈ 8 mN·m/rad（E_TPU 25 MPa）；另出 1.2 版
T_TAB = 2.0                  # 夹持片厚度
FLEX_W = 20.0                # 柔性片宽度（z）
STOP_DEG = 40.0              # 限位角
ROOT_H, ROOT_L = 5.0, 10.0   # 尾鳍根部夹块半高、长度（y，从前缘起）
Y_LE = TAB_L + GAP           # 尾鳍展中前缘 y（相对铰座块后表面 y = TAB_L... 铰座块 y ∈ [-6, TAB_L]）
# 尾鳍
AR, TAPER, SWEEP, S = 4.0, 0.5, 35.0, 1850.0   # mm²
TC = 0.21
SCREW_D = 1.6                # M2 自攻底孔
E_TPU = 25.0                 # TPU 95A 弯曲模量 [MPa]（15–40 之间，实测校正）
K_FLEX = E_TPU * FLEX_W * T_FLEX**3 / (12 * GAP) * 1e-3   # N·mm/rad → ×1e-3 … 见 report

def naca00xx(t, n=50):
    beta = np.linspace(0, math.pi, n); x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    return x, yt

def planform(n=13):
    b = math.sqrt(AR * S); zz = np.linspace(-b / 2, b / 2, n); s = np.abs(2 * zz / b)
    shape = 1 - (1 - TAPER) * s; c_mid = S / np.trapezoid(shape, zz); c = c_mid * shape
    y_le = np.abs(zz) * math.tan(math.radians(SWEEP))
    return b, c_mid, zz, c, y_le

def section_wire(zk, ck, yle):
    xc, yt = naca00xx(TC)
    up = [(X_H + ck * t, Y_LE + yle + ck * x, zk) for x, t in zip(xc, yt)]
    lo = [(X_H - ck * t, Y_LE + yle + ck * x, zk) for x, t in zip(xc[::-1], yt[::-1])]
    pts = up + lo[1:-1]
    e = cq.Edge.makeSpline([cq.Vector(*p) for p in pts], periodic=True)
    return cq.Wire.assembleEdges([e])

def build_fluke():
    b, c_mid, zz, c, y_le = planform()
    wires = [section_wire(Z_C + zk, ck, yl) for zk, ck, yl in zip(zz, c, y_le)]
    foil = cq.Solid.makeLoft(wires, ruled=False)
    foil = cq.Workplane().add(foil)
    # 根部夹块
    root = (cq.Workplane("XY").box(2 * ROOT_H, ROOT_L, 2 * BLK_W, centered=(True, False, True))
            .translate((X_H, Y_LE, Z_C)).edges("|Z and <Y").fillet(2.0).edges("|Y").chamfer(0.8))
    part = foil.union(root)
    # TPU 后片槽（从前缘面开口向后）
    slot = cq.Workplane("XY").box(T_TAB + 0.25, TAB_L + 0.3, FLEX_W + 0.3, centered=(True, False, True)).translate((X_H, Y_LE - 0.3, Z_C))
    part = part.cut(slot)
    # 自攻螺钉孔（沿 x 贯穿夹块）
    for dz in (-6.0, 6.0):
        part = part.cut(cq.Workplane("XY").cylinder(2 * ROOT_H + 2, SCREW_D / 2, direct=(1, 0, 0)).translate((X_H, Y_LE + TAB_L / 2 + 0.5, Z_C + dz)))
    return part, dict(b=b, c_mid=c_mid, c_tip=float(c[0]))

def build_strut():
    orig = cq.importers.importStep(ORIG)
    tab = orig.intersect(cq.Workplane("XY").box(TAB_END + 20, 20, 20, centered=(False, True, True)).translate((-20, 0, Z_C)))
    # 支杆：椭圆截面，从 x = 5 到铰座块
    strut = (cq.Workplane("YZ").ellipse(STRUT_A, STRUT_B).extrude(X_H - BLK_H + 1 - 5.0)
             .translate((5.0, 0, Z_C)))
    # 铰座块 y ∈ [-6, TAB_L]
    blk = cq.Workplane("XY").box(2 * BLK_H, TAB_L + 6.0, 2 * BLK_W, centered=(True, False, True)).translate((X_H, -6.0, Z_C))
    blk = blk.edges("|Z and <Y").fillet(BLK_H - 0.3).edges("|Y").chamfer(0.8)   # 迎流面整圆、展向棱倒角
    part = tab.union(strut).union(blk)
    # 限位爪：内面在 x_h ± s，向后伸到 y = Y_LE - 1
    s_in = ROOT_H * math.cos(math.radians(STOP_DEG)) + (Y_LE - (TAB_L + GAP / 2)) * math.sin(math.radians(STOP_DEG))   # 夹块前角在 40° 时的高度
    for sgn in (+1, -1):
        prong = cq.Workplane("XY").box(BLK_H - s_in, Y_LE - 1 - TAB_L, 2 * BLK_W, centered=(False, False, True))
        prong = prong.translate((X_H + (s_in if sgn > 0 else -BLK_H), TAB_L, Z_C))
        part = part.union(prong)
    # TPU 前片槽（从后表面 y = TAB_L 开口向前）
    slot = cq.Workplane("XY").box(T_TAB + 0.25, TAB_L + 0.3, FLEX_W + 0.3, centered=(True, False, True)).translate((X_H, 0.0, Z_C))
    part = part.cut(slot)
    for dz in (-6.0, 6.0):
        part = part.cut(cq.Workplane("XY").cylinder(2 * BLK_H + 2, SCREW_D / 2, direct=(1, 0, 0)).translate((X_H, TAB_L / 2, Z_C + dz)))
    return part, s_in

def build_flexure(t=None):
    t = T_FLEX if t is None else t
    f = cq.Workplane("XY").box(T_TAB, TAB_L, FLEX_W, centered=(True, False, True)).translate((X_H, 0.0, Z_C))
    f = f.union(cq.Workplane("XY").box(t, GAP + 0.4, FLEX_W, centered=(True, False, True)).translate((X_H, TAB_L - 0.2, Z_C)))
    f = f.union(cq.Workplane("XY").box(T_TAB, TAB_L, FLEX_W, centered=(True, False, True)).translate((X_H, TAB_L + GAP, Z_C)))
    for y in (TAB_L / 2, TAB_L + GAP + TAB_L / 2):
        for dz in (-6.0, 6.0):
            f = f.cut(cq.Workplane("XY").cylinder(6, 1.1, direct=(1, 0, 0)).translate((X_H, y, Z_C + dz)))
    return f

if __name__ == "__main__":
    fluke, info = build_fluke(); strut, s_in = build_strut(); flex = build_flexure()
    parts = dict(fluke=fluke, strut=strut, **{"flexure_t1.0": flex, "flexure_t1.2": build_flexure(1.2)})
    for name, p in parts.items():
        for ver, mir in (("A", False), ("B", True)):
            if name.startswith("flexure") and mir: continue   # 柔性片对称，不分 A/B
            q = p.mirror("XZ") if mir else p          # B：头朝 −y
            tag = f"{name}_{ver}" if not name.startswith("flexure") else name
            cq.exporters.export(q, os.path.join(HERE, f"{tag}.stl"), tolerance=0.02, angularTolerance=0.1)
            cq.exporters.export(q, os.path.join(HERE, f"{tag}.step"))
        v = p.val().Volume(); bb = p.val().BoundingBox()
        print(f"{name}: 体积 {v/1e3:.2f} cm³  bbox x[{bb.xmin:.1f},{bb.xmax:.1f}] y[{bb.ymin:.1f},{bb.ymax:.1f}] z[{bb.zmin:.1f},{bb.zmax:.1f}]")
    asm = cq.Assembly()
    asm.add(strut, name="strut", color=cq.Color(0.6, 0.6, 0.62)); asm.add(fluke, name="fluke", color=cq.Color(0.76, 0.1, 0.36)); asm.add(flex, name="flexure", color=cq.Color(0.83, 0.63, 0.09))
    asm.save(os.path.join(HERE, "fluke_foot_assembly_A.step"))
    k_flex = E_TPU * FLEX_W * T_FLEX**3 / (12 * GAP)      # N·mm/rad (MPa·mm³/mm)
    rep = dict(k_flex_t12_mNm_per_rad=E_TPU * FLEX_W * 1.2**3 / (12 * GAP), X_E=X_E, X_H=X_H, L_HINGE=L_HINGE, Y_LE=Y_LE, GAP=GAP, T_FLEX=T_FLEX, T_TAB=T_TAB, FLEX_W=FLEX_W, STOP_DEG=STOP_DEG, s_in=s_in,
               k_flex_mNm_per_rad=k_flex, E_TPU=E_TPU, strain_at_25deg_pct=100 * math.radians(25) / GAP * T_FLEX / 2, fluke=info,
               x_p_mm=-(GAP / 2))
    json.dump(rep, open(os.path.join(HERE, "fluke_foot_cad.json"), "w"), indent=1)
    print(json.dumps(rep, indent=1, ensure_ascii=False))
