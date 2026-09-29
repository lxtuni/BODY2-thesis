#!/usr/bin/env python3
"""升力型（lift-based）——用 BODY2 原来的双连杆（摇臂 + 平行四边形）实现。
桨板枢轴放在 1/4 弦（不在中间），髋摆角 psi 做上下拍动，脚踝只负责把攻角控住。
准定常受力，机理示意，非 CFD。
"""
import numpy as np

# ---- 机构几何（取自 leg_dynamics.py，单位 m）
O2       = np.array([0.0, 0.0])      # 髋（摇臂）转轴，平移到原点
L_ROCK_B = 0.0301                    # O2→B 摇臂长端
L_ROCK_C = 0.0147                    # O2→C 摇臂短端 = 平行四边形地杆
L_PAR    = 0.039                     # 从动杆（平行四边形长边）
STEM, HEAD = 0.050, 0.042            # 桨杆长、桨头长（= 弦长 chord）
SPAN, THK  = 0.044, 0.0056           # 桨头宽（= 展长 span）、厚
PHI_EXT    = np.deg2rad(105.0)       # 平行四边形张角，划水位姿

# ---- 升力型运动参数
RHO   = 1000.0
U     = 0.25                         # 前进速度 [m/s]
T     = 1.0; F = 1/T; OM = 2*np.pi*F
PSI_M = np.deg2rad(214.0)      # 摆动中心：让桨板枢轴的圆弧最接近竖直（前伸水平翼）
PSI_A = np.deg2rad(30.0)             # 髋摆幅 ±30°
PIVOT = 0.25   # 枢轴距前缘的弦比例（0.5 = 中间）。前缘是"外侧"端（迎流），
               # 所以枢轴沿桨轴的位置 = STEM + (1-PIVOT)*HEAD
AL0   = np.deg2rad(13.0)             # 目标攻角幅值
S     = HEAD * SPAN                  # 迎流面积
AR    = SPAN / HEAD
CLA   = 2*np.pi*AR/(2 + np.sqrt(AR**2 + 4))
AL_ST = np.deg2rad(16.0)

def rot(a): return np.array([np.cos(a), np.sin(a)])
def wrap(a): return (a + np.pi) % (2*np.pi) - np.pi

def psi_of(t):  return PSI_M + PSI_A*np.sin(OM*t)

def linkage(psi, phi=PHI_EXT):
    """摇臂角 psi（O2→C）+ 平行四边形张角 phi → 各关节点与桨轴方向 u"""
    C = O2 + L_ROCK_C*rot(psi)
    D = O2 + L_PAR*rot(psi + phi)
    E = C + (D - O2)                      # 平行四边形闭合
    B = O2 + L_ROCK_B*rot(psi + np.pi)
    u = rot(psi + np.pi)                  # 桨轴（E → 桨尖），与摇臂平行
    return dict(B=B, C=C, D=D, E=E, u=u)

def piv_r(pivot=None):
    """枢轴沿桨轴、距 E 的距离"""
    return STEM + (1 - (PIVOT if pivot is None else pivot))*HEAD

def pivot_pos(t, pivot=None):
    """桨板俯仰枢轴的位置（前缘在外侧，枢轴在距前缘 PIVOT 弦处）"""
    lk = linkage(psi_of(t))
    return lk["E"] + piv_r(pivot)*lk["u"], lk

def state(t, u_inf=U, pivot=PIVOT, al0=AL0, dt=2e-4):
    P,  lk = pivot_pos(t, pivot)
    P1, _  = pivot_pos(t + dt, pivot)
    P0, _  = pivot_pos(t - dt, pivot)
    V = (P1 - P0)/(2*dt)                          # 枢轴速度（含髋摆带来的前后分量）
    W = np.array([-u_inf, 0.0]) - V               # 板感受到的相对来流
    Wm = np.hypot(*W); ŵ = W/Wm
    gW = np.arctan2(W[1], W[0])
    # 脚踝的任务：把攻角控到 al_target（随扫动速度变号，端点自然归零）
    al_t = al0 * np.clip(V[1]/(PSI_A*OM*piv_r(pivot)), -1, 1)
    g_c  = gW - al_t - np.pi                      # 弦线方向（前缘→后缘）
    psi  = psi_of(t)
    th   = wrap(g_c - psi - np.pi)                # 脚踝俯仰角（相对桨轴）
    al   = al_t
    ale  = np.clip(al, -AL_ST, AL_ST)
    CL = CLA*np.sin(ale)
    CD = 0.02 + CL**2/(np.pi*AR*0.8) + 1.17*np.sin(al)**2
    q  = 0.5*RHO*Wm**2*S
    n̂  = np.array([-ŵ[1], ŵ[0]])
    Fv = q*CL*n̂ + q*CD*ŵ
    # 压心位置（弦比例）：小攻角 1/4 弦，接近失速时后移
    xcp = 0.25 + 0.20*min(abs(al)/AL_ST, 1.0)
    M_piv = q*CL*HEAD*(xcp - pivot)               # 绕枢轴的水动力俯仰力矩
    d̂ = rot(g_c)
    return dict(P=P, lk=lk, V=V, W=W, Wm=Wm, ŵ=ŵ, n̂=n̂, al=al, th=th, psi=psi,
                L=q*CL, D=q*CD, F=Fv, Fx=Fv[0], Fz=Fv[1], M=M_piv, d̂=d̂, xcp=xcp)

if __name__ == "__main__":
    tt = np.linspace(0, T, 721)
    S_ = [state(t) for t in tt]
    Pz = np.array([s["P"][1] for s in S_]); Px = np.array([s["P"][0] for s in S_])
    th = np.degrees([s["th"] for s in S_]); al = np.degrees([s["al"] for s in S_])
    Fx = np.array([s["Fx"] for s in S_]);   M  = np.array([s["M"] for s in S_])
    Mmid = np.array([state(t, pivot=0.5)["M"] for t in tt])
    thr  = np.degrees(np.gradient(np.radians(th), tt))
    print(f"枢轴：距前缘 {PIVOT:.2f} 弦 = {PIVOT*HEAD*1e3:.1f} mm，沿桨轴距 E {piv_r()*1e3:.1f} mm（中点会是 {(STEM+0.5*HEAD)*1e3:.1f} mm），板 {HEAD*1e3:.0f}×{SPAN*1e3:.0f} mm, AR {AR:.2f}")
    print(f"髋摆角 psi = {np.degrees(PSI_M):.0f}° ± {np.degrees(PSI_A):.0f}°")
    print(f"桨板枢轴 峰峰：垂直 {(Pz.max()-Pz.min())*1e3:.1f} mm，水平 {(Px.max()-Px.min())*1e3:.1f} mm")
    print(f"St = {F*(Pz.max()-Pz.min())/U:.3f}")
    print(f"脚踝俯仰角 θ ∈ [{th.min():+.1f}, {th.max():+.1f}]°，峰值角速度 {np.abs(thr).max():.0f} °/s")
    print(f"攻角 α 幅值 {np.abs(al).max():.1f}°")
    print(f"平均推力 {Fx.mean()*1e3:+.1f} mN（准定常，非 CFD）")
    print(f"脚踝力矩峰值：1/4 弦 {np.abs(M).max()*1e3:.1f} mN·m   中点 {np.abs(Mmid).max()*1e3:.1f} mN·m"
          f"  → 差 {np.abs(Mmid).max()/max(np.abs(M).max(),1e-9):.1f} 倍")
