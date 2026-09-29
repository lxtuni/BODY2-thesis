#!/usr/bin/env python3
"""E 点（平行四边形末端 = 桨铰轴）沿竖直线沉浮的两自由度逆解（ψ 摇臂 + φ 平行四边形）
 E = O2 + L_RC·rot(ψ) + L_PAR·rot(ψ+φ)，φ ∈ (12°,168°)
 机体坐标：x 向前、z 向上，O2 = (63.4, 4) mm，水线 z = −10 mm
 杆长：CAD 值 (14.7, 39) 与照片估计值 (40, 60) 都算；真实值请量。
"""
import math, json, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
O2 = np.array([63.4, 4.0]); WL = -10.0
PHI_MIN, PHI_MAX = 12.0, 168.0

def ik(E, L_rc, L_par, branch=+1):
    """返回 (ψ, φ) [deg]；branch=+1: φ 取正解（E 在 C 的逆时针侧）"""
    d = E - O2; r = np.hypot(*d)
    cphi = (r**2 - L_rc**2 - L_par**2) / (2 * L_rc * L_par)
    if abs(cphi) > 1: return None
    phi = branch * math.degrees(math.acos(cphi))
    # ψ：E − O2 = L_rc rot(ψ) + L_par rot(ψ+φ) → 角 ψ = atan2(d) − atan2(L_par sinφ, L_rc + L_par cosφ)
    psi = math.degrees(math.atan2(d[1], d[0]) - math.atan2(L_par * math.sin(math.radians(phi)), L_rc + L_par * math.cos(math.radians(phi))))
    return (psi + 180) % 360 - 180, phi

def fk(psi, phi, L_rc, L_par):
    p, q = math.radians(psi), math.radians(psi + phi)
    C = O2 + L_rc * np.array([math.cos(p), math.sin(p)])
    E = C + L_par * np.array([math.cos(q), math.sin(q)])
    return C, E

def traj(x_E, z_c, h0, L_rc, L_par, f, n=360):
    t = np.linspace(0, 1 / f, n, endpoint=False)
    z = z_c + h0 * np.sin(2 * math.pi * f * t)
    out = []
    for zi in z:
        s = None
        for br in (+1, -1):
            r = ik(np.array([x_E, zi]), L_rc, L_par, br)
            if r and PHI_MIN <= abs(r[1]) <= PHI_MAX:
                s = r; break
        if s is None: return None
        out.append(s)
    out = np.array(out); psi, phi = out[:, 0], out[:, 1]
    psi = np.degrees(np.unwrap(np.radians(psi)))
    dt = t[1] - t[0]
    return dict(t=t, z=z, psi=psi, phi=phi,
                psi_rate=float(np.abs(np.gradient(psi, dt)).max()), phi_rate=float(np.abs(np.gradient(phi, dt)).max()),
                psi_rng=(float(psi.min()), float(psi.max())), phi_rng=(float(phi.min()), float(phi.max())))

def arc_traj(psi_fix, phi_mid, dphi, L_rc, L_par, f, n=360):
    """ψ 固定、只动 φ：E 在以 C 为圆心的圆弧上（切向竖直时 ≈ 沉浮 + 二阶纵荡）"""
    t = np.linspace(0, 1 / f, n, endpoint=False)
    phi = phi_mid + dphi * np.sin(2 * math.pi * f * t)
    E = np.array([fk(psi_fix, p, L_rc, L_par)[1] for p in phi])
    dt = t[1] - t[0]
    return dict(t=t, x=E[:, 0], z=E[:, 1], psi=np.full_like(t, psi_fix), phi=phi,
                phi_rate=float(np.abs(np.gradient(phi, dt)).max()), psi_rate=0.0,
                surge=float((E[:, 0].max() - E[:, 0].min()) / 2), h0=float((E[:, 1].max() - E[:, 1].min()) / 2))

if __name__ == "__main__":
    for tag, (L_rc, L_par) in (("CAD", (14.7, 39.0)), ("照片估计", (40.0, 60.0))):
        print(f"\n=== 杆长 {tag}: L_RC {L_rc}, L_PAR {L_par}  (E 最深 z = {O2[1]-L_rc-L_par:.0f} mm) ===")
        print("竖直直线沉浮（ψ+φ 联动）：")
        for x_E in (O2[0] - 30, O2[0], O2[0] + 30, O2[0] + 60):
            for z_c in (-30, -40, -50, -60, -70):
                for h0 in (15, 20, 25):
                    r = traj(x_E, z_c, h0, L_rc, L_par, 1.5)
                    if r:
                        print(f"  x_E {x_E:5.1f} z_c {z_c:4d} h0 {h0}: ψ {r['psi_rng'][0]:6.0f}…{r['psi_rng'][1]:6.0f} (摆 {r['psi_rng'][1]-r['psi_rng'][0]:3.0f}°)  φ {r['phi_rng'][0]:5.0f}…{r['phi_rng'][1]:5.0f} (摆 {r['phi_rng'][1]-r['phi_rng'][0]:3.0f}°)  峰值角速度 1.5 Hz: ψ {r['psi_rate']:4.0f} φ {r['phi_rate']:4.0f} °/s")
        print("只动 φ（ψ 固定，圆弧沉浮）：")
        for psi_fix in (180, 225, 270):
            for phi_mid in (90, 120, 135):
                for dphi in (15, 20, 25):
                    if not (PHI_MIN <= phi_mid - dphi and phi_mid + dphi <= PHI_MAX): continue
                    r = arc_traj(psi_fix, phi_mid, dphi, L_rc, L_par, 1.5)
                    Cc, Em = fk(psi_fix, phi_mid, L_rc, L_par)
                    print(f"  ψ {psi_fix} φ {phi_mid}±{dphi}: E 中心 ({Em[0]:.0f},{Em[1]:.0f})  h0 {r['h0']:.1f} mm  纵荡 ±{r['surge']:.1f} mm  φ 峰值角速度 {r['phi_rate']:.0f} °/s")
