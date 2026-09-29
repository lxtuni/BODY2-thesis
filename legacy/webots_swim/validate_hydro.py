#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_hydro.py  --  量化 Webots 内置流体模型与"真实水动力"的差距
Quantify the gap between the Webots built-in fluid model and real hydrodynamics,
for the robot0526 swim demo. Pure Python, no Webots needed.

参照三层 / three reference tiers:
  (1) 解析叶素-阻力模型 (blade-element drag, Li 2019 / Qu 2025)
  (2) 无量纲相似性 (Re, Fr, St)
  (3) 缺失物理项估计 (added mass, unsteady, lift, free surface, ...)
"""
import math

# ---------- physical constants ----------
rho = 1000.0      # water density [kg/m^3]
mu  = 1.3e-3      # dynamic viscosity [Pa.s]
g   = 9.81

# ---------- robot / gait (matches swim_demo.wbt + paddle_swim.py) ----------
A_paddle = 0.05 * 0.07        # paddle broad-face area [m^2]
w_paddle, L_paddle = 0.05, 0.07
r_hip2pad = 0.11              # hip pivot -> paddle center [m]
body_w, body_h_sub = 0.10, 0.025   # body width, submerged height (~50%)
f    = 0.7                    # paddling frequency [Hz]
hip_amp = 0.50               # hip swing amplitude [rad]
Cd_flat = 1.10               # flat-plate drag coeff (Li 2019: 1.17)
Cd_body = 1.10

v_webots = 0.015             # MEASURED steady speed in Webots [m/s]

# ---------- (A) paddle kinematics ----------
Omega_pk = hip_amp * 2*math.pi*f          # peak angular velocity [rad/s]
alpha_pk = hip_amp * (2*math.pi*f)**2     # peak angular acceleration [rad/s^2]
v_pad_pk = Omega_pk * r_hip2pad           # peak paddle speed [m/s]
a_pad_pk = alpha_pk * r_hip2pad           # peak paddle accel [m/s^2]

# ---------- (B) blade-element quasi-steady drag thrust ----------
F_drag_peak = 0.5*rho*Cd_flat*A_paddle*v_pad_pk**2        # peak [N]
F_drag_avg  = 0.5*F_drag_peak                              # <sin^2>=0.5 over power stroke
thrust_net  = 2*F_drag_avg*0.5   # ~2 legs in power, ~50% duty -> rough net [N]

# ---------- (C) ADDED MASS (the term Webots omits entirely) ----------
# flat plate added mass ~ rho * pi * (w/2)^2 * length   (normal to broad face)
m_added = rho*math.pi*(w_paddle/2)**2*L_paddle
F_added_pk = m_added*a_pad_pk                              # peak added-mass force [N]

# ---------- (D) terminal speed from analytical balance: thrust = body drag ----------
A_body = body_w*body_h_sub
V_term = math.sqrt(thrust_net/(0.5*rho*Cd_body*A_body))    # [m/s]

# ---------- (E) dimensionless numbers ----------
Re_pad = rho*v_pad_pk*w_paddle/mu                          # paddle Reynolds
Fr     = v_webots/math.sqrt(g*0.22)                        # body Froude (L=0.22)
stroke_excursion = 2*hip_amp*r_hip2pad                     # peak-peak paddle travel [m]
St     = f*stroke_excursion/max(v_webots,1e-6)            # Strouhal

def line(): print("-"*64)
print("="*64); print(" robot0526 水动力保真度核对 / hydro-fidelity check"); print("="*64)
print(f"步态: f={f} Hz, hip_amp={hip_amp} rad  ->  桨峰值速度 {v_pad_pk:.3f} m/s, 峰值加速度 {a_pad_pk:.2f} m/s^2")
line()
print(" (1) 解析叶素-阻力参照 / blade-element reference")
print(f"   单桨峰值推力 F=1/2 rho Cd A v^2  = {F_drag_peak*1000:6.1f} mN")
print(f"   估计净推力                       ~ {thrust_net*1000:6.1f} mN")
print(f"   解析终端速度 V_term              = {V_term:6.3f} m/s")
print(f"   Webots 实测速度                  = {v_webots:6.3f} m/s")
print(f"   >> 速度比 (解析/仿真)            = {V_term/v_webots:5.1f} x")
line()
print(" (2) 无量纲数 / dimensionless numbers")
print(f"   Re(桨)  = {Re_pad:8.0f}   (~1e4, 平板 Cd~1.1 适用; 与同尺度真实机器人同量级 -> 阻力定律可迁移)")
print(f"   Fr(体)  = {Fr:8.3f}   (<<1 -> 兴波阻力可忽略, 但平面水面仍无法溅水/造波)")
print(f"   St      = {St:8.2f}   (高效游动应 ~0.2-0.4; 这里远超 -> 桨在'打滑'/效率低)")
line()
print(" (3) Webots 完全没建模的项 / physics omitted by the built-in fluid")
print(f"   * 附加质量 added mass: m_add~{m_added*1000:.0f} g = {m_added/0.030:.1f}x 桨自重;")
print(f"       峰值附加质量力 F_add = m_add*a = {F_added_pk*1000:.0f} mN")
print(f"       >> 与阻力推力同量级甚至更大 ({F_added_pk/ max(F_drag_peak,1e-9):.1f}x)! 这是最大缺口")
print("   * 非定常/涡脱 (Wagner, 涡环推力): 准稳态模型完全没有")
print("   * 升力/环量 lift: 只有沿轴阻力, 无环量升力")
print("   * 自由液面兴波/溅水: 水面是静止平面")
print("   * Cd 取常数: 不随 Re 变化; 桨用长方体投影面积, 非真实蹼形")
print("   * 腿间流场干扰: 各连杆独立, 无非线性叠加 (Qu 2025 指出真实是非线性)")
print("   * 人工阻尼污染: defaultDamping(linear0.3/ang0.6) 是为稳定加的非物理阻力")
print("="*64)
