#!/usr/bin/env python3
"""摇橹式游泳轨迹 → 电控用的足端 J 轨迹 + 关节角 q1/q2 + 各腿舵机角
 几何 = 固件 Leg_GeometryV2（AB 30, BC 35, CE 52, EH 16.2, HI 52, IJ 51, AH 50），腿平面坐标 y 向前、z 向上，A = 髋轴。
 运动：长杆 H–I 保持竖直，杆（A–H 与脚 I–J 同向，角度 q2）绕 A 摆 θ(t) = θ_mean + θ0·sin(2πft)（θ_mean = −18°：顶端水平、底端 −36°）
       → J = (101·cosθ, −52 + 101·sinθ)：以 (0,−52) 为圆心、半径 101 的圆弧，中点切线竖直。
 输出：scull_traj_Nppc.csv（每周期 N 点）、scull_gait_patch.c（可直接替换 swim_leg_pose 的代码）
"""
import math, csv, json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
AB, BC, CE, EH, HI, IJ, AH = 30.0, 35.0, 52.0, 16.2, 52.0, 51.0, 50.0
TH0_DEG = 18.0            # 杆摆幅（→ J 沉浮 ±31 mm）
TH_MEAN_DEG = -18.0       # 杆平均角：冲程顶端杆水平（0°），底端 −36°；使尾鳍在冲程顶端仍在船底（−35）以下 ≥ 10 mm
F_HZ = 1.5
N = 100
ZC = -(HI)                # 圆心高度（长杆竖直时 I 在 H 正下方 52）
R = AH + IJ               # 101
MAPS = dict(LF=(145.0, -1.0, 53.0, -1.0), RF=(32.0, 1.0, 142.0, 1.0), LH=(120.0, -1.0, 48.0, -1.0), RH=(59.0, 1.0, 134.0, 1.0))   # offset1,sign1,offset2,sign2（固件当前值）
PHASE_OFF = dict(LF=0.0, RF=0.5, LH=0.5, RH=0.0)

def ik_v2(yJ, zJ, bq2=1.0, bq1=1.0):
    """固件 Leg_Solve2DIK_v2 的逐行移植"""
    rAJ = math.hypot(yJ, zJ); AHpIJ = AH + IJ
    c = (rAJ**2 + AHpIJ**2 - HI**2) / (2 * AHpIJ * rAJ)
    if abs(c) > 1: return None
    q2 = math.atan2(zJ, yJ) + bq2 * math.acos(c)
    c2, s2 = math.cos(q2), math.sin(q2)
    dy, dz = yJ - AHpIJ * c2, zJ - AHpIJ * s2; rHJ = math.hypot(dy, dz)
    yE, zE = AH * c2 - EH * dy / rHJ, AH * s2 - EH * dz / rHJ
    zE_AB = zE - AB; rBE = math.hypot(yE, zE_AB)
    c1 = (BC**2 + rBE**2 - CE**2) / (2 * BC * rBE)
    if abs(c1) > 1: return None
    q1 = math.atan2(zE_AB, yE) + bq1 * math.acos(c1)
    return math.degrees(q1), math.degrees(q2), (yE, zE)

def fk(q1, q2):
    """LinkageSim 正解（固件尺寸）：返回 H, E, I, J"""
    e2 = np.array([math.cos(math.radians(q2)), math.sin(math.radians(q2))]); e1 = np.array([math.cos(math.radians(q1)), math.sin(math.radians(q1))])
    B = np.array([0.0, AB]); H = AH * e2; C = B + BC * e1
    d = H - C; dist = np.linalg.norm(d); a = (CE**2 - EH**2 + dist**2) / (2 * dist); h = math.sqrt(max(CE**2 - a**2, 0))
    mid = C + a * d / dist; perp = np.array([-d[1], d[0]]) / dist
    E = mid + h * perp                      # side +1（与 Python LinkageSim 一致；用 IK 的 E 校验）
    u = (H - E) / np.linalg.norm(H - E); I = H + HI * u; J = I + IJ * e2
    return H, E, I, J

def traj(phase):
    th = math.radians(TH_MEAN_DEG) + math.radians(TH0_DEG) * math.sin(2 * math.pi * phase)
    return R * math.cos(th), ZC + R * math.sin(th), math.degrees(th)

if __name__ == "__main__":
    rows = []; ok = True
    for i in range(N):
        ph = i / N; yJ, zJ, th = traj(ph)
        s = ik_v2(yJ, zJ)
        if s is None: ok = False; print("IK fail at phase", ph); continue
        q1, q2, (yE, zE) = s
        H, E, I, J = fk(q1, q2)
        row = dict(phase=ph, t_s=ph / F_HZ, yJ=yJ, zJ=zJ, stem_deg=th, q1=q1, q2=q2, yE=yE, zE=zE, yI=I[0], zI=I[1], fk_err=float(np.hypot(J[0] - yJ, J[1] - zJ)))
        for leg, (o1, s1, o2, s2) in MAPS.items():
            row[f"theta1_{leg}"] = o1 + s1 * q1; row[f"theta2_{leg}"] = o2 + s2 * q2
        rows.append(row)
    q1s = np.array([r["q1"] for r in rows]); q2s = np.array([r["q2"] for r in rows]); dt = 1 / (F_HZ * N)
    rate1 = np.abs(np.gradient(q1s, dt)).max(); rate2 = np.abs(np.gradient(q2s, dt)).max()
    print(f"θ0 {TH0_DEG}° → J: y {min(r['yJ'] for r in rows):.1f}…{max(r['yJ'] for r in rows):.1f}, z {min(r['zJ'] for r in rows):.1f}…{max(r['zJ'] for r in rows):.1f}")
    print(f"q2 {q2s.min():.1f}…{q2s.max():.1f}°（= 杆角）  q1 {q1s.min():.1f}…{q1s.max():.1f}°  峰值角速度 @{F_HZ} Hz: q1 {rate1:.0f} °/s, q2 {rate2:.0f} °/s；@2 Hz ×{2/F_HZ:.2f}")
    print(f"I（安装片前孔）z 范围 {min(r['zI'] for r in rows):.1f}…{max(r['zI'] for r in rows):.1f}；FK 复核最大误差 {max(r['fk_err'] for r in rows):.3f} mm")
    for leg in MAPS:
        t1 = [r[f"theta1_{leg}"] for r in rows]; t2 = [r[f"theta2_{leg}"] for r in rows]
        flag = "" if (min(t1) >= 0 and max(t1) <= 180 and min(t2) >= 0 and max(t2) <= 180) else "  ← 超出 0–180，需重装舵机臂/改 offset"
        print(f"  {leg}: theta1 {min(t1):.0f}…{max(t1):.0f}   theta2 {min(t2):.0f}…{max(t2):.0f}{flag}")
    fn = os.path.join(HERE, f"scull_traj_{N}ppc.csv")
    with open(fn, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}) for r in rows]
    print("→", fn)
    # C 补丁：替换 swim_leg_pose 为正弦圆弧
    patch = f"""/* ===== 摇橹式（尾鳍脚）游泳轨迹：替换 gait_task.c 里的 swim_leg_pose() =====
 * 长杆 H–I 保持竖直，杆 A–H（= 脚 I–J 方向，q2）绕髋轴 A 正弦摆动 scull_theta_mean_deg ± scull_theta0_deg（−18° ± 18°：顶端杆水平，底端 −36°），
 * 这样冲程顶端尾鳍仍在船底以下（船底 y = −35，尾鳍顶端 ≈ −48）。
 * 足端 J 走以 (0, −HI) 为圆心、半径 AH+IJ 的圆弧，中点切线竖直 → 尾鳍铰轴上下沉浮 ±{R*math.sin(math.radians(TH0_DEG)):.0f} mm。
 * 频率仍由 swim_freq_hz / CH5 控制（建议 1.5 Hz，最高 2 Hz），相位偏移沿用 swim_phase_off_*（对角同相）。
 * 船的前进方向 = 尾鳍拖行的反方向（脚朝前时船向后游）。 */
volatile float scull_theta0_deg = {TH0_DEG:.1f}f;     /* 杆摆幅；18° → 沉浮 ±31 mm，22° → ±38 mm（2 Hz 时舵机 q1 峰值 ~{rate1*2/F_HZ*22/18:.0f} °/s，注意上限） */
volatile float scull_theta_mean_deg = {TH_MEAN_DEG:.1f}f;  /* 杆平均角（负 = 向下）；改成 0 就回到对称摆动，但尾鳍会顶到船底 */
volatile float scull_center_z   = -{HI:.1f}f;         /* 圆心高度 = −HI（长杆竖直） */
volatile float scull_radius     = {R:.1f}f;           /* AH + IJ */

static void scull_leg_pose(float phase, float *yJ, float *zJ, int *in_power)
{{
    float p  = phase - floorf(phase);
    float th = deg_to_rad(scull_theta_mean_deg) + deg_to_rad(scull_theta0_deg) * sinf(2.0f * LEG_PI * p);
    *yJ = scull_radius * cosf(th);
    *zJ = scull_center_z + scull_radius * sinf(th);
    if (in_power) *in_power = (th > deg_to_rad(scull_theta_mean_deg));           /* 仅用于调试灯，升力型没有“划水相” */
}}
/* 在 MODE 5 的 walk_trot_select == 2 分支里，把四个 swim_leg_pose(...) 改成 scull_leg_pose(...) 即可。
 * IK、舵机映射、相位推进都不用改。 */
"""
    open(os.path.join(HERE, "scull_gait_patch.c"), "w").write(patch)
    json.dump(dict(theta0_deg=TH0_DEG, theta_mean_deg=TH_MEAN_DEG, f_hz=F_HZ, radius=R, center_z=ZC, q1_range=[float(q1s.min()), float(q1s.max())], q2_range=[float(q2s.min()), float(q2s.max())],
                   rate_q1_deg_s=float(rate1), rate_q2_deg_s=float(rate2)), open(os.path.join(HERE, "scull_traj_summary.json"), "w"), indent=1)
    print("→ scull_gait_patch.c")
