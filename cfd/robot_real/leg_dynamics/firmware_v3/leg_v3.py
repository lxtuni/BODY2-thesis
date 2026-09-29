"""新腿型（gait_task.c Leg_GeometryV3 / Leg_Solve2DIK_v3 逐行移植）
   O1=(-21.8,22) O2=(0,0)，曲柄 O1A=30，连杆 AC=35，摆杆 O2C=30，大腿 O2E=30，小腿 EM=45
   C = -30·dir(q_r)，E = 30·dir(gamma)，M = E - 45·dir(q_r)   （+x 前进，+y 上，原点 O2）"""
import math, numpy as np
O1 = np.array([-21.8, 22.0]); O1A, AC, O2C, O2E, EM = 30.0, 35.0, 30.0, 30.0, 45.0
def d(a): return np.array([math.cos(math.radians(a)), math.sin(math.radians(a))])
def ang(u, v): return math.degrees(math.acos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
def fourbar_from_rocker(qr, branch_crank=1.0):
    C = -O2C * d(qr); v = C - O1; dd = np.linalg.norm(v)
    ca = (O1A**2 + dd**2 - AC**2) / (2 * O1A * dd)
    if abs(ca) > 1: return None
    crank = math.atan2(v[1], v[0]) - branch_crank * math.acos(ca)
    A = O1 + O1A * np.array([math.cos(crank), math.sin(crank)])
    q1 = (math.degrees(crank) - 180 + 180) % 360 - 180
    return dict(q1=q1, A=A, C=C, ACO2=ang(A - C, -C), O1C=dd)
def ik(xM, yM, branch_rocker=1.0, branch_crank=1.0):
    r = math.hypot(xM, yM); k = (O2E**2 - EM**2 - r**2) / (2 * r * EM)
    if abs(k) > 1: return None
    qr = math.degrees(math.atan2(yM, xM) + branch_rocker * math.acos(k))
    E = np.array([xM, yM]) + EM * d(qr); gamma = (math.degrees(math.atan2(E[1], E[0])) + 180) % 360 - 180
    fb = fourbar_from_rocker(qr, branch_crank)
    if fb is None: return None
    return dict(q1=fb["q1"], gamma=gamma, qr=(qr + 180) % 360 - 180, E=E, ACO2=fb["ACO2"], O2EM=ang(-E, np.array([xM, yM]) - E), r=r)
def fk(q1, gamma):
    """q1, gamma → M（摆杆角由四杆正解：A 已知，C 在 O2 圆上且 |AC|=35，两解取传动角合规且与 IK 一致者）"""
    A = O1 - O1A * d(q1); out = []
    dd = np.linalg.norm(A); a = (O2C**2 - AC**2 + dd**2) / (2 * dd); h2 = O2C**2 - a**2
    if h2 < 0: return out
    h = math.sqrt(h2); u = A / dd; p = np.array([-u[1], u[0]])
    for s in (1, -1):
        C = a * u + s * h * p; qr = math.degrees(math.atan2(-C[1], -C[0]))
        E = O2E * d(gamma); M = E - EM * d(qr)
        out.append(dict(qr=qr, M=M, E=E, C=C, A=A, ACO2=ang(A - C, -C)))
    return out
def ok(s): return s and 10 <= s["ACO2"] <= 170 and s["O2EM"] >= 45 and 15 <= s["r"] <= 75
if __name__ == "__main__":
    user = [(23.8, -49.0), (20.6, -65.7), (37.4, -58.3), (51.1, -47.1), (32.5, -10.4), (29.7, -25.6), (26.6, -36.8)]
    orig = [(66*math.cos(math.radians(p)), 66*math.sin(math.radians(p))) for p in (-60, -85, -110, -135)] + [None] + [(50.8*math.cos(math.radians(p)), 50.8*math.sin(math.radians(p))) for p in (-101, -80)]
    for name, pts in (("用户输入的 7 点", user), ("注释里的原始点（|O2M|=66/50.8, φ 从 +x）", orig)):
        print("==", name)
        for i, p in enumerate(pts):
            if p is None: print(f" cp4: 由锁定圆决定"); continue
            s = ik(*p)
            print(f" cp{i}: ({p[0]:6.1f},{p[1]:6.1f}) |O2M|={math.hypot(*p):5.1f}  " + (f"q1={s['q1']:6.1f} γ={s['gamma']:7.1f} q_r={s['qr']:7.1f} 膝={s['O2EM']:5.1f}° 传动角={s['ACO2']:5.1f}° {'OK' if ok(s) else '✗ 违反判据'}" if s else "✗ 不可达"))
