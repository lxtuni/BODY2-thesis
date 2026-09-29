#!/usr/bin/env python3
"""A3 = A2 的尾鳍整体外移：翼根（销）在 z = 0，翼面从 z = −15 伸到 +45（原 −30…+30）。
   目的：尾鳍内侧尖离船体中线 ≥ 55 mm，脱开船底舭部，让杆可以绕水平位置摆 → 推力纯向前。
   其余（安装片、支杆、横梁、叶簧、叉耳、销、指节、T 臂）与 A2 完全相同。"""
import os, math, json, numpy as np, cadquery as cq
import fluke_straight_cad as base
HERE = os.path.dirname(os.path.abspath(__file__))
Z_IN, Z_OUT, S, TAPER, SWEEP = -15.0, 45.0, 1150.0, 0.5, 30.0
def planform(n=13):
    zz = np.linspace(Z_IN, Z_OUT, n)
    shape = 1 - (1 - TAPER) * np.abs(zz) / Z_OUT          # 弦长从翼根线性收到外梢 0.5；内侧短翼同斜率
    c_root = S / np.trapezoid(shape, zz); c = c_root * shape
    return Z_OUT - Z_IN, c_root, zz, c, np.abs(zz) * math.tan(math.radians(SWEEP))   # 后掠顶点在翼根 z=0
base.planform = planform
if __name__ == "__main__":
    part, strut, fluke, info = base.build()
    cq.exporters.export(part, os.path.join(HERE, "fluke_straight_A3.stl"), tolerance=0.02, angularTolerance=0.1)
    cq.exporters.export(part, os.path.join(HERE, "fluke_straight_A3.step"))
    pin_solid = cq.Workplane("XY").cylinder(base.EAR_GAP + 2 * base.EAR_T + 2.0, 1.75 / 2, direct=(0, 0, 1)).translate((base.X_PIN, 0, base.Z_C))
    asm = cq.Assembly(); asm.add(part, name="fluke_straight_A3", color=cq.Color(0.7, 0.7, 0.72)); asm.add(pin_solid, name="pin_filament_1.75", color=cq.Color(0.9, 0.2, 0.2))
    asm.save(os.path.join(HERE, "fluke_straight_A3_with_pin.step"))
    bb = part.val().BoundingBox()
    rep = dict(z_in=Z_IN, z_out=Z_OUT, span=Z_OUT - Z_IN, c_root=info["c_mid"], c_tip=info["c_tip"], S_mm2=S, volume_cm3=part.val().Volume() / 1e3, mass_g_PETG=part.val().Volume() / 1e3 * 1.27, bbox=[bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax])
    json.dump(rep, open(os.path.join(HERE, "fluke_straight_A3.json"), "w"), indent=1); print(json.dumps(rep, indent=1))
