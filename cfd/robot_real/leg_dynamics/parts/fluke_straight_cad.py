#!/usr/bin/env python3
"""直排一体件（A2，摇橹式，缩小版）：原安装片 → 支杆 → 耗材销铰 + 两片 PETG 叶簧 → 月牙尾鳍，全 PETG 单件
 零件坐标 = 原扇形足：x 沿杆（装机时水平朝前），y 尾鳍厚度方向（竖直），z 展向
 铰轴：x = X_PIN，Ø1.75 耗材销穿 z 向（叉耳 — 尾鳍根部指节 — 叉耳）
 弹簧：两片叶簧（0.5 × 8 × 46）在 z = z_c ± 10.5，从支杆横梁伸到尾鳍前缘；尾鳍绕销转 θ 时叶簧末端在 y 向弯 9·θ
 限位：尾鳍指节上的 T 形挡臂碰支杆端面，±40°
"""
import os, math, json, numpy as np
import cadquery as cq

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(HERE, "fan_foot_orig.step")
X_E = -12.21; Z_C = 1.5
TAB_END = 8.0
STEM_Y, STEM_Z = 6.0, 3.0            # 支杆椭圆截面（y 竖直 × z 横向）
X_PIN = 50.0                          # 铰轴（距 E 62 mm）
X_STEM_END = 44.0                     # 支杆端面 = 限位面（销前 6 mm）
EAR_T, EAR_GAP = 2.5, 6.4             # 叉耳厚 2.5，耳间距 6.4（指节 5.6 + 每侧 0.4 间隙）
EAR_Y = 3.75
PIN_D = 1.9                           # Ø1.75 耗材，打印孔收缩留量
KNUCKLE_R = 3.25
T_ARM = (-3.2, 5.5)                   # 限位 T 臂末端相对销 (dx, ±dy) → 碰 x = X_STEM_END 面时 ±40°
X_BAR = (14.0, 18.0)                  # 横梁 x 范围（叶簧根部）
WEB_T, WEB_B, WEB_ZOFF = 0.50, 6.0, 10.5  # 叶簧厚、宽（z）、中心离 z_c 的偏距（窄一点减少上下拍水面积，厚一点补刚度）
WEB_END = 64.0                        # 叶簧伸到尾鳍前缘里
X_LE = 53.0                           # 尾鳍展中前缘（销后 3 mm）
AR, TAPER, SWEEP, S, TC = 3600.0 / 1150.0, 0.5, 30.0, 1150.0, 0.21   # A2：展 60、11.5 cm²（A1 是展 86、18.5 cm²）
E_PETG = 2000.0

def naca00xx(t, n=50):
    beta = np.linspace(0, math.pi, n); x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    return x, yt

def planform(n=13):
    b = math.sqrt(AR * S); zz = np.linspace(-b / 2, b / 2, n); s = np.abs(2 * zz / b)
    shape = 1 - (1 - TAPER) * s; c_mid = S / np.trapezoid(shape, zz); c = c_mid * shape
    return b, c_mid, zz, c, np.abs(zz) * math.tan(math.radians(SWEEP))

def section_wire(zk, ck, xle):
    xc, yt = naca00xx(TC)
    up = [(X_LE + xle + ck * x, ck * t, zk) for x, t in zip(xc, yt)]
    lo = [(X_LE + xle + ck * x, -ck * t, zk) for x, t in zip(xc[::-1], yt[::-1])]
    e = cq.Edge.makeSpline([cq.Vector(*p) for p in up + lo[1:-1]], periodic=True)
    return cq.Wire.assembleEdges([e])

def build():
    orig = cq.importers.importStep(ORIG)
    tab = orig.intersect(cq.Workplane("XY").box(TAB_END + 20, 20, 20, centered=(False, True, True)).translate((-20, 0, Z_C)))
    # 支杆 x 5 → 44，截面 8×4 圆角
    stem = cq.Workplane("YZ").ellipse(STEM_Y / 2, STEM_Z / 2).extrude(X_STEM_END - 5.0).translate((5.0, 0, Z_C))
    # 横梁（叶簧根）
    # 横梁：x–y 面内椭圆截面（弦 6 沿 x 顺流，厚 2.4 沿 y），沿 z 拉伸，像一片小翼
    zw = 2 * (WEB_ZOFF + WEB_B / 2) + 1.0
    bar = cq.Workplane("XY").ellipse((X_BAR[1] - X_BAR[0]) / 2, 1.2).extrude(zw).translate(((X_BAR[0] + X_BAR[1]) / 2, 0, Z_C - zw / 2))
    # 叉耳（两片），在支杆末端向后伸到销后 3 mm
    ears = None
    for zc in (Z_C - EAR_GAP / 2 - EAR_T / 2, Z_C + EAR_GAP / 2 + EAR_T / 2):
        e = cq.Workplane("XY").box(X_PIN + 3.0 - (X_STEM_END - 4.0), 2 * EAR_Y, EAR_T, centered=(False, True, True)).translate((X_STEM_END - 4.0, 0, zc))
        e = e.edges("|Z and >X").fillet(EAR_Y - 0.2)
        ears = e if ears is None else ears.union(e)
    # 叶簧 ×2
    webs = None
    for zc in (Z_C - WEB_ZOFF, Z_C + WEB_ZOFF):
        w = cq.Workplane("XY").box(WEB_END - X_BAR[1] + 1.0, WEB_T, WEB_B, centered=(False, True, True)).translate((X_BAR[1] - 1.0, 0, zc))
        webs = w if webs is None else webs.union(w)
    # 叉座：把两只耳和支杆末端连成一体（x 40–44，跨整个 z 宽度）；其后表面 x = 44 就是 T 臂的限位面
    fork_base = cq.Workplane("XY").box(4.0, 2 * EAR_Y, EAR_GAP + 2 * EAR_T, centered=(False, True, True)).translate((X_STEM_END - 4.0, 0, Z_C)).edges("|Z and <X").fillet(1.5).edges("|Y and <X").fillet(1.2)
    strut = tab.union(stem).union(bar).union(fork_base).union(ears).union(webs)
    # 尾鳍 + 指节 + 桥 + 限位 T 臂
    b, c_mid, zz, c, xle = planform()
    foil = cq.Workplane().add(cq.Solid.makeLoft([section_wire(Z_C + zk, ck, xl) for zk, ck, xl in zip(zz, c, xle)], ruled=False))
    kn = cq.Workplane("XY").cylinder(5.6, KNUCKLE_R, direct=(0, 0, 1)).translate((X_PIN, 0, Z_C))
    # 桥：从指节到尾鳍前缘并埋入翼内（到 x = 61），上下后缘做长楔形倒角，让它顺着翼型表面没入，不露台阶
    bridge = (cq.Workplane("XY").box(61.0 - X_PIN, 5.0, 5.6, centered=(False, True, True)).translate((X_PIN, 0, Z_C))
              .edges("|Z and >X").chamfer(2.2, 7.0))
    tarm = cq.Workplane("XY").box(2.0, 2 * T_ARM[1], 5.6, centered=(True, True, True)).translate((X_PIN + T_ARM[0], 0, Z_C))
    fluke = foil.union(kn).union(bridge).union(tarm)
    # 销孔贯穿耳与指节
    pin = cq.Workplane("XY").cylinder(EAR_GAP + 2 * EAR_T + 1.0, PIN_D / 2, direct=(0, 0, 1)).translate((X_PIN, 0, Z_C))   # 只穿耳和指节，不碰叶簧
    strut = strut.cut(pin); fluke = fluke.cut(pin)
    part = strut.union(fluke)          # 一体 STL（叶簧把两者连成一件；销孔为通孔）
    return part, strut, fluke, dict(b=b, c_mid=c_mid, c_tip=float(c[0]))

if __name__ == "__main__":
    part, strut, fluke, info = build()
    cq.exporters.export(part, os.path.join(HERE, "fluke_straight_A2.stl"), tolerance=0.02, angularTolerance=0.1)
    cq.exporters.export(part, os.path.join(HERE, "fluke_straight_A2.step"))
    pin_solid = cq.Workplane("XY").cylinder(EAR_GAP + 2 * EAR_T + 2.0, 1.75 / 2, direct=(0, 0, 1)).translate((X_PIN, 0, Z_C))
    cq.exporters.export(pin_solid, os.path.join(HERE, "pin_1p75.stl"))
    asm = cq.Assembly(); asm.add(part, name="fluke_straight_A2", color=cq.Color(0.7, 0.7, 0.72)); asm.add(pin_solid, name="pin_filament_1.75", color=cq.Color(0.9, 0.2, 0.2))
    asm.save(os.path.join(HERE, "fluke_straight_A2_with_pin.step"))
    L = WEB_END - X_BAR[1]; arm = X_LE + WEB_ZOFF * math.tan(math.radians(SWEEP)) - X_PIN   # 叶簧末端离销的 x 距离
    EI = E_PETG * WEB_B * WEB_T**3 / 12; k_web = 12 * EI / L**3; k_rot = 2 * k_web * arm**2
    th_stop = math.degrees(math.pi - math.atan2(T_ARM[1], T_ARM[0]) - math.acos((X_PIN - X_STEM_END) / math.hypot(*T_ARM)))
    rep = dict(X_E=X_E, X_PIN=X_PIN, pin_from_E=X_PIN - X_E, X_LE=X_LE, x_p_mm=X_PIN - X_LE, web=dict(t=WEB_T, b=WEB_B, L=L, arm=arm, k_web_N_per_mm=k_web),
               k_rot_mNm_per_rad=k_rot, strain_pct_at_35deg=100 * 6 * WEB_T * arm * math.radians(35) / L**2, strain_pct_at_27deg=100 * 6 * WEB_T * arm * math.radians(27) / L**2,
               stop_deg=th_stop, fluke=info, volume_cm3=part.val().Volume() / 1e3, mass_g_PETG=part.val().Volume() / 1e3 * 1.27)
    bb = part.val().BoundingBox(); rep["bbox"] = [bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax]
    json.dump(rep, open(os.path.join(HERE, "fluke_straight_A2.json"), "w"), indent=1)
    print(json.dumps(rep, indent=1, ensure_ascii=False))
