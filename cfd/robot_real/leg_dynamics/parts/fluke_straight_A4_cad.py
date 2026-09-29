#!/usr/bin/env python3
"""A4 = A2 的连杆整体向外（+z）平移 DZ：安装片不动，支杆在 x 8…12 处做一个横向曲拐，
   之后的支杆、横梁、叶簧、叉耳、销、指节、尾鳍全部原样、只是 z 方向外移 DZ。尾鳍仍以杆为中心（±30）。"""
import os, math, json, cadquery as cq
import fluke_straight_cad as b
HERE = os.path.dirname(os.path.abspath(__file__))
DZ = 15.0                                   # 外移量：尾鳍内侧尖离船体中线 40 → 55 mm
X_J0, X_J1 = 8.0, 12.0                      # 曲拐横杆的 x 范围
def build():
    Z = b.Z_C + DZ
    orig = cq.importers.importStep(b.ORIG)
    tab = orig.intersect(cq.Workplane("XY").box(b.TAB_END + 20, 20, 20, centered=(False, True, True)).translate((-20, 0, b.Z_C)))
    stub = cq.Workplane("YZ").ellipse(b.STEM_Y / 2, b.STEM_Z / 2).extrude(X_J1 - 5.0).translate((5.0, 0, b.Z_C))                  # 安装片后的短杆（原 z）
    jog = (cq.Workplane("XY").box(X_J1 - X_J0, b.STEM_Y, DZ + b.STEM_Z + 0.6, centered=(False, True, False)).translate((X_J0, 0, b.Z_C - b.STEM_Z / 2 - 0.3)))                                                                                                 # 横向曲拐 4 × 6 × (DZ+3)
    stem = cq.Workplane("YZ").ellipse(b.STEM_Y / 2, b.STEM_Z / 2).extrude(b.X_STEM_END - X_J0).translate((X_J0, 0, Z))              # 外移后的支杆
    zw = 2 * (b.WEB_ZOFF + b.WEB_B / 2) + 1.0
    bar = cq.Workplane("XY").ellipse((b.X_BAR[1] - b.X_BAR[0]) / 2, 1.2).extrude(zw).translate(((b.X_BAR[0] + b.X_BAR[1]) / 2, 0, Z - zw / 2))
    ears = None
    for zc in (Z - b.EAR_GAP / 2 - b.EAR_T / 2, Z + b.EAR_GAP / 2 + b.EAR_T / 2):
        e = cq.Workplane("XY").box(b.X_PIN + 3.0 - (b.X_STEM_END - 4.0), 2 * b.EAR_Y, b.EAR_T, centered=(False, True, True)).translate((b.X_STEM_END - 4.0, 0, zc)).edges("|Z and >X").fillet(b.EAR_Y - 0.2)
        ears = e if ears is None else ears.union(e)
    webs = None
    for zc in (Z - b.WEB_ZOFF, Z + b.WEB_ZOFF):
        w = cq.Workplane("XY").box(b.WEB_END - b.X_BAR[1] + 1.0, b.WEB_T, b.WEB_B, centered=(False, True, True)).translate((b.X_BAR[1] - 1.0, 0, zc))
        webs = w if webs is None else webs.union(w)
    fork_base = cq.Workplane("XY").box(4.0, 2 * b.EAR_Y, b.EAR_GAP + 2 * b.EAR_T, centered=(False, True, True)).translate((b.X_STEM_END - 4.0, 0, Z)).edges("|Z and <X").fillet(1.5).edges("|Y and <X").fillet(1.2)
    strut = tab.union(stub).union(jog).union(stem).union(bar).union(fork_base).union(ears).union(webs)
    bb_, c_mid, zz, c, xle = b.planform()
    foil = cq.Workplane().add(cq.Solid.makeLoft([b.section_wire(Z + zk, ck, xl) for zk, ck, xl in zip(zz, c, xle)], ruled=False))
    kn = cq.Workplane("XY").cylinder(5.6, b.KNUCKLE_R, direct=(0, 0, 1)).translate((b.X_PIN, 0, Z))
    bridge = cq.Workplane("XY").box(61.0 - b.X_PIN, 5.0, 5.6, centered=(False, True, True)).translate((b.X_PIN, 0, Z)).edges("|Z and >X").chamfer(2.2, 7.0)
    tarm = cq.Workplane("XY").box(2.0, 2 * b.T_ARM[1], 5.6, centered=(True, True, True)).translate((b.X_PIN + b.T_ARM[0], 0, Z))
    fluke = foil.union(kn).union(bridge).union(tarm)
    pin = cq.Workplane("XY").cylinder(b.EAR_GAP + 2 * b.EAR_T + 1.0, b.PIN_D / 2, direct=(0, 0, 1)).translate((b.X_PIN, 0, Z))
    strut = strut.cut(pin); fluke = fluke.cut(pin)
    return strut.union(fluke), Z
if __name__ == "__main__":
    part, Z = build()
    for name, p in (("fluke_straight_A4", part), ("fluke_straight_A4_mirror", part.mirror("XY", basePointVector=(0, 0, b.Z_C)))):
        cq.exporters.export(p, os.path.join(HERE, name + ".stl"), tolerance=0.02, angularTolerance=0.1)
        cq.exporters.export(p, os.path.join(HERE, name + ".step"))
    bb = part.val().BoundingBox()
    rep = dict(DZ=DZ, fluke_z_center=Z, fluke_z_range=[Z - 30, Z + 30], volume_cm3=part.val().Volume() / 1e3, mass_g_PETG=part.val().Volume() / 1e3 * 1.27, bbox=[bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax])
    json.dump(rep, open(os.path.join(HERE, "fluke_straight_A4.json"), "w"), indent=1); print(json.dumps(rep, indent=1))
