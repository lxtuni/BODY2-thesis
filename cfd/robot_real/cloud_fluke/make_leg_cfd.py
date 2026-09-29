#!/usr/bin/env python3
"""
BODY2 单腿非定常 CFD 算例生成器（第 2 级）—— 重叠网格 overset + overPimpleDyMFoam
==================================================================================
桨板按 leg_dynamics.py 的运动学（髋摆动 ψ(t) + 平行四边形 φ(t)）做刚体运动，运动表由本脚本导出为
tabulated6DoFMotion 的 6DoF.dat；脚蹼两态用两次计算包络：--fan open（全程张开）/ --fan closed（全程收拢），
后处理时划水相取 open、回收相取 closed 合成。单相水、层流（桨头 Re≈10⁴ 过渡区）。

    python3 make_leg_cfd.py --case ~/run/leg_open   --fan open   --U 0
    python3 make_leg_cfd.py --case ~/run/leg_closed --fan closed --U 0
    步态参数与 leg_dynamics.py 相同，用环境变量：DUTY=0.3 FAN_CLOSE=0.15 SWEEP=90 PHI_RET=60 python3 make_leg_cfd.py ...
    其他：--cycles 3  --bg 0.004（背景网格最细尺寸）--comp 0.0025（部件网格基准尺寸）--turb laminar|kOmegaSST --np 8

坐标：机体坐标系，x 前进、z 向上、y 左；O2（下舵机轴）在 (0.0634, 0, 0.004)。来流从 +x 吹向 −x（U>0 时）。
"""
import argparse, os, sys, json, struct, shutil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True)
ap.add_argument("--fan", choices=["open", "closed", "plate"], required=True, help="open/closed = 第一套折扇两态；plate = 第二套刚性矩形板（配 --feather 做脚踝羽化）")
ap.add_argument("--plate-w", type=float, default=0.044, help="刚性板宽 [m]（默认 0.044 → 与 120° 扇等面积 18.5 cm²；0.073 = 扇尖同宽）")
ap.add_argument("--plate-l", type=float, default=None, help="刚性板长 [m]（默认 = HEAD 42 mm）")
ap.add_argument("--feather", default=None, help="脚踝羽化 TC,TO,TF[,lin]：τ_c 开始绕桨轴转 90°（边缘迎流），τ_o 开始转回；过渡 TF 秒。第 4 项 lin = 匀速（限速舵机，如 7465W 0.07 s/60° → 90° 0.105 s），缺省 cos = 余弦斜坡。例 0.351,-0.066,0.11,lin")
ap.add_argument("--U", type=float, default=0.0)
ap.add_argument("--cycles", type=float, default=3.0)
ap.add_argument("--bg", type=float, default=0.004, help="背景网格在桨扫掠区的单元尺寸 [m]")
ap.add_argument("--comp", type=float, default=0.0025, help="部件网格基准单元尺寸 [m]（表面再细化 1 级）")
ap.add_argument("--turb", default="laminar", choices=["laminar", "kOmegaSST"])
ap.add_argument("--np", type=int, default=8)
ap.add_argument("--nsub", type=int, default=800, help="运动表每周期采样点数")
ap.add_argument("--maxCo", type=float, default=0.8)
ap.add_argument("--writes", type=int, default=20, help="每周期写出帧数")
ap.add_argument("--vox", type=float, default=None, help="洞切体素尺寸 [m]（默认 min(桨厚/4.5, comp/2) ≈ 1.2 mm；不要大于 1.5 mm）")
ap.add_argument("--leg", action="store_true", help="把划水相里没入水中的曲柄2、从动杆段（按零件实测截面）作为刚体附在桨上一起算（只对 --fan open 有意义）")
ap.add_argument("--wl", type=float, default=-0.010, help="水线 z [m]，用来决定 --leg 时杆件的湿长（船体中部 −0.010）")
ap.add_argument("--traj", default=None, help="直接给桨头中心轨迹 CSV（results/papers/optimal_trajectory.py 的输出）：覆盖 leg_dynamics 的步态；列 t_s,x_abs_m,z_abs_m,theta_deg,fan_open,phase")
a = ap.parse_args()

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leg_dynamics as ld          # 运动学 + 参数（读环境变量）

CASE = os.path.abspath(os.path.expanduser(a.case))
T = ld.T_PER
Y_LEG = 0.0                         # 单腿算例：把腿放在 y=0 的对称面上（机体坐标下 FL 腿在 y=+0.066，这里只算单腿）
RHO, NU = ld.RHO, ld.NU

# ------------------------------------------------------------------ 运动学采样
TRAJ = None
if a.traj:
    # 轨迹模式：桨头中心 P(t)、桨轴角 θ（0 = 桨尖朝下）→ 桨根 E = P − (STEM + HEAD/2)·u，ψ = 90° − θ（u = rot(ψ+π)）；φ 只用于 --leg，取 PHI_EXT
    D = np.genfromtxt(a.traj, delimiter=",", names=True, dtype=None, encoding="utf-8")
    tc = D["t_s"]; T = float(tc[-1] + (tc[1] - tc[0]))
    TRAJ = dict(t=tc, x=D["x_abs_m"], z=D["z_abs_m"], th=D["theta_deg"], fan=D["fan_open"], power=(D["phase"] == "power").astype(float))
    print(f"--traj {a.traj}：{len(tc)} 行，周期 T = {T:.4f} s（覆盖 leg_dynamics 的 T/DUTY/SWEEP）")
nT = int(a.nsub * a.cycles) + 1
t = np.linspace(0, a.cycles * T, nT)
if TRAJ is None:
    psi, phi, power = ld.gait(t)
    P = [ld.pose(p, q) for p, q in zip(psi, phi)]
    E = np.array([p["E"] for p in P]); u = np.array([p["u"] for p in P])
else:
    tm = t % T; ip = lambda y: np.interp(tm, TRAJ["t"], y, period=T)
    th = np.radians(ip(TRAJ["th"])); u = np.stack([-np.sin(th), -np.cos(th)], 1)
    Pc = np.stack([ip(TRAJ["x"]), ip(TRAJ["z"])], 1); E = Pc - (ld.STEM + ld.HEAD / 2) * u
    psi = np.pi / 2 - th; phi = np.full_like(psi, np.deg2rad(ld.PHI_EXT)); power = ip(TRAJ["power"]) > 0.5
    P = [dict(E=E[k], u=u[k], n=np.array([-u[k, 1], u[k, 0]]), th2=psi[k] + phi[k]) for k in range(nT)]
E0, psi0 = E[0], psi[0]
u0 = u[0]; n0 = np.array([-u0[1], u0[0]])
tip = E + (ld.STEM + ld.HEAD) * u
ry = -np.degrees(psi - psi0)                     # 绕 +y 右手转角；x→z 的逆时针（x-z 图）= 绕 +y 的负转角
# 脚踝羽化 β(t)：绕桨自身轴（过 E0、方向 u0）转 90° → 桨板与 x–z 面共面，对面内任何回程运动都是边缘迎流
beta = np.zeros(nT); FEATHER = None
if a.feather:
    fparts = a.feather.split(","); fc, fo, tf = [float(v) for v in fparts[:3]]; fshape = fparts[3] if len(fparts) > 3 else "cos"
    FEATHER = dict(TAU_CLOSE=fc, TAU_OPEN=fo, TF=tf, SHAPE=fshape)
    tau_ = (t / T) % 1.0; tfr = tf / T
    def ramp(x):                                                     # 0→1，x 为相位差/过渡长；lin = 限速舵机匀速，cos = 余弦斜坡（峰值速率 1.57×）
        return np.clip(x, 0, 1) if fshape == "lin" else 0.5 - 0.5 * np.cos(np.pi * np.clip(x, 0, 1))
    dc = (tau_ - fc) % 1.0; do = (tau_ - fo) % 1.0; Lw = (fc - fo) % 1.0
    beta = np.where(do < Lw, 90.0 * (1 - ramp(do / tfr)), 90.0 * ramp(dc / tfr))   # 张开窗 [fo, fc)：从 90 转回 0；其余：从 0 转到 90
    if tf > 0.9 * T * min(Lw, 1 - Lw): print(f"!! 过渡 {tf*1e3:.0f} ms 比开/闭窗口还长，β 达不到 0 或 90°")
    from scipy.spatial.transform import Rotation as Rot
    uw = np.array([u0[0], 0.0, u0[1]])
    Mt = [Rot.from_euler("y", ry[k], degrees=True) * Rot.from_rotvec(uw * np.radians(beta[k])) for k in range(nT)]
    eul = np.array([m.as_euler("XYZ", degrees=True) for m in Mt])   # OpenFOAM tabulated6DoFMotion: q = qx·qy·qz → R = Rx·Ry·Rz = scipy 内旋 XYZ
    eul = np.degrees(np.unwrap(np.radians(eul), axis=0))
    rx, ry, rz = eul[:, 0], eul[:, 1], eul[:, 2]
    print(f"--feather：τ_c {fc}，τ_o {fo}，过渡 {tf*1e3:.0f} ms（{fshape}，均速 {90/tf:.0f} °/s）；欧拉角范围 rx[{rx.min():.1f},{rx.max():.1f}] ry[{ry.min():.1f},{ry.max():.1f}] rz[{rz.min():.1f},{rz.max():.1f}]")
else:
    rx = np.zeros(nT); rz = np.zeros(nT)

# ------------------------------------------------------------------ 桨板 STL（世界坐标，t=0 位姿）
def paddle_outline(fan_open):
    """桨面内轮廓 (lateral y, s 沿桨轴自 E 起)"""
    S, H, SW, HW = ld.STEM, ld.HEAD, ld.STEM_W, ld.HEAD_W
    if a.fan == "plate":
        Wp_, Lp_ = a.plate_w, (a.plate_l or ld.HEAD)
        return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S), (Wp_ / 2, S), (Wp_ / 2, S + Lp_), (-Wp_ / 2, S + Lp_), (-SW / 2, S)]
    if fan_open:
        ang = np.deg2rad(ld.FAN_DEG); th = np.linspace(ang / 2, -ang / 2, 48)   # 从右角到左角，避免多边形自交
        fan = [(ld.FAN_R * np.sin(x), S + ld.FAN_R * np.cos(x)) for x in th]
        return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + fan + [(-SW / 2, S)]
    k = np.linspace(0, 1, 24)
    sh = [(HW / 2 * np.sin(np.pi / 2 * q), S + 0.025 * (1 - np.cos(np.pi / 2 * q))) for q in k]
    right = [(w, z) for w, z in sh if w > SW / 2]
    return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + right + [(HW / 2, S + H), (-HW / 2, S + H)] + [(-w, z) for w, z in right[::-1]] + [(-SW / 2, S)]

def extrude_world(outline, thk):
    """挤出并放到世界坐标：局部 (lat, s, nrm) → E0 + s·u0 + nrm·n0 + lat·ŷ；逐面按已知外法向定向"""
    Pl = np.array(outline, dtype=float)
    # 多边形取逆时针（在 (lat, s) 平面内）
    area2 = np.sum(Pl[:, 0] * np.roll(Pl[:, 1], -1) - np.roll(Pl[:, 0], -1) * Pl[:, 1])
    if area2 < 0: Pl = Pl[::-1]
    c = Pl.mean(0); N = len(Pl)
    u0w = np.array([u0[0], 0.0, u0[1]]); n0w = np.array([n0[0], 0.0, n0[1]]); yw = np.array([0.0, 1.0, 0.0])
    def Wp(lat, s, nrm):
        return E0[0] + s * u0[0] + nrm * n0[0], lat + Y_LEG, E0[1] + s * u0[1] + nrm * n0[1]
    tris, want = [], []
    for i in range(N):
        (la, sa), (lb, sb) = Pl[i], Pl[(i + 1) % N]
        tris.append([Wp(c[0], c[1], thk / 2), Wp(la, sa, thk / 2), Wp(lb, sb, thk / 2)]); want.append(n0w)
        tris.append([Wp(c[0], c[1], -thk / 2), Wp(la, sa, -thk / 2), Wp(lb, sb, -thk / 2)]); want.append(-n0w)
        dlat, ds = lb - la, sb - sa
        out = ds * yw - dlat * u0w                       # 逆时针多边形的边外法向 (ds, -dlat) 映射到世界
        p1, p2 = Wp(la, sa, thk / 2), Wp(lb, sb, thk / 2); q1, q2 = Wp(la, sa, -thk / 2), Wp(lb, sb, -thk / 2)
        tris.append([p1, q1, q2]); want.append(out)
        tris.append([p1, q2, p2]); want.append(out)
    tri = np.array(tris); want = np.array(want)
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    flip = np.einsum("ij,ij->i", nrm, want) < 0
    tri[flip] = tri[flip][:, ::-1]
    return tri

def write_stl(tri, fn):
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
    rec = np.zeros(len(tri), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])); rec["n"] = nrm; rec["v"] = tri
    open(fn, "wb").write(b"\0" * 80 + struct.pack("<I", len(tri)) + rec.tobytes())

tri = extrude_world(paddle_outline(a.fan == "open"), ld.THK)
pmin, pmax = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
half_w = (pmax[1] - pmin[1]) / 2

# ------------------------------------------------------------------ 连杆段（--leg）：曲柄2 从 D 指向 O2、从动杆从 E 指向 C，方向都是 rot(θ2+π)
# 零件实测（parts/*.stl）：曲柄2 横向厚 7.8 × 面内宽 9.0 mm，从动杆 3.0 × 5.0 mm；横向都在桨的 −y 侧，与桨杆略有重叠以便合成一个封闭体
def box_world(origin, d, length, width, y0, y1):
    """面内矩形（沿 d 长 length、沿 rot(d,90°) 宽 width，居中）× 横向 [y0,y1] 的六面体，返回带外法向的 12 个三角形"""
    m = np.array([-d[1], d[0]])
    def Wp(s, q, y): p = origin + s * d + q * m; return (p[0], y + Y_LEG, p[1])
    c = [Wp(s, q, y) for s in (0.0, length) for q in (-width / 2, width / 2) for y in (y0, y1)]   # index = 4*is + 2*iq + iy
    faces = [((0, 1, 3, 2), -np.array([d[0], 0, d[1]])), ((4, 6, 7, 5), np.array([d[0], 0, d[1]])),
             ((0, 4, 5, 1), -np.array([m[0], 0, m[1]])), ((2, 3, 7, 6), np.array([m[0], 0, m[1]])),
             ((0, 2, 6, 4), np.array([0, -1.0, 0])), ((1, 5, 7, 3), np.array([0, 1.0, 0]))]
    tris, want = [], []
    for (i, j, k, l), nrm in faces:
        tris += [[c[i], c[j], c[k]], [c[i], c[k], c[l]]]; want += [nrm, nrm]
    tri_ = np.array(tris); want = np.array(want)
    nn = np.cross(tri_[:, 1] - tri_[:, 0], tri_[:, 2] - tri_[:, 0]); flip = np.einsum("ij,ij->i", nn, want) < 0
    tri_[flip] = tri_[flip][:, ::-1]
    return tri_

tri_links = None; leg_info = {}
if a.leg and a.fan == "closed":
    print("!! --leg 只对 --fan open 有意义：closed 算例用于回收相，那时曲柄2/从动杆在水面以上，附上去只会多出虚假阻力 → 忽略 --leg"); a.leg = False
if a.leg:
    th2_0 = P[0]["th2"]; d2 = ld.rot(th2_0 + np.pi)
    D0 = ld.O2 + ld.L_PAR * ld.rot(th2_0)
    # 湿长：在划水相（φ 不变，杆与桨刚性同转）内、水线 a.wl 以下的最大长度；零件比销距长的部分各伸一半
    def wet_max(origin_t, dir_t, pins, part, mask):
        ext = (part - pins) / 2; s = np.linspace(-ext, pins + ext, 80); best = 0.0
        for k in np.where(mask)[0][::4]:
            z = origin_t[k][1] + s * dir_t[k][1]; wet = (z < a.wl).sum() * (s[1] - s[0]); best = max(best, wet)
        return best
    if TRAJ is not None: sys.exit("--leg 与 --traj 不能同时用（轨迹模式没有连杆几何）")
    t1 = np.linspace(0, T, 801); psi1, phi1, power1 = ld.gait(t1)           # 单独取一个完整周期算湿长（与 --cycles 无关）
    P1 = [ld.pose(p, q) for p, q in zip(psi1, phi1)]
    th2_t = np.array([p["th2"] for p in P1]); dir_t = np.stack([np.cos(th2_t + np.pi), np.sin(th2_t + np.pi)], 1)
    D_t = ld.O2 + ld.L_PAR * np.stack([np.cos(th2_t), np.sin(th2_t)], 1)
    E_t = np.array([p["E"] for p in P1])
    Lc2 = wet_max(D_t, dir_t, ld.L_PAR, 0.047, power1); Lf = wet_max(E_t, dir_t, ld.L_PAR, 0.045, power1)
    if Lc2 < 0.003 and Lf < 0.003:
        print("!! --leg：划水相里曲柄2 / 从动杆都不入水（水线 %.3f），不加连杆" % a.wl); a.leg = False
    else:
        Lc2 += 0.002; Lf += 0.002
        parts = []
        if Lc2 > 0.003: parts.append(box_world(D0 - 0.004 * d2, d2, Lc2 + 0.004, 0.009, -0.0090, 0.0000)); leg_info["crank2_len"] = Lc2
        if Lf > 0.003:  parts.append(box_world(E0 - 0.004 * d2, d2, Lf + 0.004, 0.005, -0.0045, 0.0000)); leg_info["follower_len"] = Lf
        tri_links = np.concatenate(parts)
        lmin, lmax = tri_links.reshape(-1, 3).min(0), tri_links.reshape(-1, 3).max(0)
        print(f"--leg：曲柄2 段 {Lc2*1000:.1f} mm、从动杆段 {Lf*1000:.1f} mm 附在桨上（水线 {a.wl:+.3f}，划水相最大湿长 + 2 mm）；连杆 STL 包围盒 x[{lmin[0]:.3f},{lmax[0]:.3f}] z[{lmin[2]:.3f},{lmax[2]:.3f}]")

# ------------------------------------------------------------------ 部件网格盒（与桨轴对齐的旋转六面体）
m_ax, m_n, m_lat = 0.016, 0.016, 0.012
s_lo, s_hi = -m_ax, ld.STEM + ld.HEAD + m_ax
n_lo, n_hi = -(ld.THK / 2 + m_n), ld.THK / 2 + m_n
w_lo, w_hi = -(half_w + m_lat), half_w + m_lat
if tri_links is not None:                      # 盒子扩到把连杆段也包进去（局部 (s, n, lat) 坐标）
    Pl = tri_links.reshape(-1, 3); rel = Pl[:, [0, 2]] - E0
    s_l = rel @ u0; n_l = rel @ n0; y_l = Pl[:, 1] - Y_LEG
    s_lo, s_hi = min(s_lo, s_l.min() - m_ax), max(s_hi, s_l.max() + m_ax)
    n_lo, n_hi = min(n_lo, n_l.min() - m_n), max(n_hi, n_l.max() + m_n)
    w_lo, w_hi = min(w_lo, y_l.min() - m_lat), max(w_hi, y_l.max() + m_lat)
def box_pt(s, nn, lat):
    p = E0 + s * u0 + nn * n0; return (p[0], lat + Y_LEG, p[1])
# blockMesh 顶点顺序：底面 (0..3) 逆时针，顶面 (4..7)；这里“底/顶”取法向 -n0 / +n0
V = [box_pt(s_lo, n_lo, w_lo), box_pt(s_hi, n_lo, w_lo), box_pt(s_hi, n_lo, w_hi), box_pt(s_lo, n_lo, w_hi),
     box_pt(s_lo, n_hi, w_lo), box_pt(s_hi, n_hi, w_lo), box_pt(s_hi, n_hi, w_hi), box_pt(s_lo, n_hi, w_hi)]
n_s = int(round((s_hi - s_lo) / a.comp)); n_w = int(round((w_hi - w_lo) / a.comp)); n_n = int(round((n_hi - n_lo) / a.comp))
# 检查 hex 手性（右手）：(v1-v0)×(v3-v0)·(v4-v0) > 0
v0, v1, v3, v4 = map(np.array, (V[0], V[1], V[3], V[4]))
if np.dot(np.cross(v1 - v0, v3 - v0), v4 - v0) < 0:
    V = V[4:] + V[:4]

# 部件盒在整个运动中的包络（用于确定背景域 / 加密盒）
corners_local = np.array([(s, nn, lat) for s in (s_lo, s_hi) for nn in (n_lo, n_hi) for lat in (w_lo, w_hi)])
env_min = np.full(3, np.inf); env_max = np.full(3, -np.inf)
for k in range(0, nT, max(1, nT // 400)):
    uk = u[k]; nk = np.array([-uk[1], uk[0]]); cb, sb_ = np.cos(np.radians(beta[k])), np.sin(np.radians(beta[k]))
    for s, nn, lat in corners_local:
        nn2, lat2 = nn * cb - lat * sb_, nn * sb_ + lat * cb        # 羽化：部件盒绕桨轴随板转
        p = E[k] + s * uk + nn2 * nk
        w = np.array([p[0], lat2 + Y_LEG, p[1]])
        env_min = np.minimum(env_min, w); env_max = np.maximum(env_max, w)

# ------------------------------------------------------------------ 背景域
pad_in = 0.03                                   # 加密盒对部件包络的余量
ref2 = (env_min - pad_in, env_max + pad_in)     # level-2 加密盒（a.bg）
ref1 = (env_min - 3 * pad_in, env_max + 3 * pad_in)
L = 0.262
dom_min = np.array([env_min[0] - 1.4 * L, min(env_min[1] - 0.9 * L, -0.25), env_min[2] - 1.0 * L])
dom_max = np.array([env_max[0] + 1.0 * L, max(env_max[1] + 0.9 * L, 0.25), max(env_max[2] + 0.12, 0.06)])
# 洞切搜索盒：运动包络 + 3 cm 余量。体素必须明显小于桨厚（5.6 mm）：沙箱实测 1.2 mm 正常、1.85 mm 及以上全域变洞。
sb_min = env_min - 0.03; sb_max = env_max + 0.03
vox = a.vox if a.vox else min(ld.THK / 4.5, a.comp / 2)
sbd = [int(min(240, max(40, np.ceil((sb_max[i] - sb_min[i]) / vox)))) for i in range(3)]
base = a.bg * 4                                  # 背景 blockMesh 基准（两级八叉树后到 a.bg）
nx, ny, nz = [max(8, int(round((dom_max[i] - dom_min[i]) / base))) for i in range(3)]

# ------------------------------------------------------------------ 6DoF 运动表
trans = E - E0
# 时间步、写出
tip_speed = np.max(np.linalg.norm(np.gradient(tip, t, axis=0), axis=1))
dt0 = 0.5 * a.comp / 2 / max(tip_speed, 0.05) * a.maxCo    # 表面单元 comp/2
dt0 = float(f"{dt0:.1e}")
endTime = a.cycles * T
writeInterval = T / a.writes

# ------------------------------------------------------------------ 写文件
def W(path, txt):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w", encoding="utf-8").write(txt)
def hdr(cls, obj, loc="system"):
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
"""
if os.path.exists(CASE):
    print(f"!! {CASE} 已存在，先删除"); shutil.rmtree(CASE)
os.makedirs(CASE)
fmt = lambda v: "(" + " ".join(f"{x:.6g}" for x in v) + ")"

# ---- 部件网格子算例 paddleMesh
PM = os.path.join(CASE, "paddleMesh")
os.makedirs(os.path.join(PM, "constant", "triSurface"))
write_stl(tri, os.path.join(PM, "constant", "triSurface", "paddle.stl"))
if tri_links is not None: write_stl(tri_links, os.path.join(PM, "constant", "triSurface", "links.stl"))
LINKS_GEOM = "    links.stl { type triSurfaceMesh; name links; }\n" if tri_links is not None else ""
LINKS_FEAT = ' { file "links.eMesh"; level 1; }' if tri_links is not None else ""
LINKS_SURF = "        links { level (1 1); patchInfo { type wall; inGroups (wall); } }\n" if tri_links is not None else ""
WALLS = '"(paddle|links)"' if tri_links is not None else "paddle"
LINKS_FO = f"""    forcesLinks               // 只有连杆段的力（桨 = forces − forcesLinks）
    {{
        type            forces;
        libs            ("libforces.so");
        patches         (links);
        rho             rhoInf;
        rhoInf          {RHO};
        CofR            ({ld.O2[0]:.4f} 0 {ld.O2[1]:.4f});
        writeControl    timeStep;
        writeInterval   1;
        log             false;
    }}
""" if tri_links is not None else ""
W(os.path.join(PM, "system", "blockMeshDict"), hdr("dictionary", "blockMeshDict") + f"""
scale 1;
vertices
(
{chr(10).join('    ' + fmt(v) for v in V)}
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({n_s} {n_w} {n_n}) simpleGrading (1 1 1)
);
edges ();
boundary
(
    oversetPatch
    {{
        type overset;
        faces
        (
            (0 4 7 3)
            (1 2 6 5)
            (0 1 5 4)
            (3 7 6 2)
            (0 3 2 1)
            (4 5 6 7)
        );
    }}
);
mergePatchPairs ();
""")
W(os.path.join(PM, "system", "surfaceFeatureExtractDict"), hdr("dictionary", "surfaceFeatureExtractDict") + "".join(f"""
{stl}
{{
    extractionMethod    extractFromSurface;
    includedAngle       150;
    subsetFeatures {{ nonManifoldEdges no; openEdges yes; }}
    writeObj            no;
}}
""" for stl in (["paddle.stl", "links.stl"] if tri_links is not None else ["paddle.stl"])))
lim = box_pt(s_hi - 0.004, n_hi - 0.004, w_hi - 0.004)     # 部件盒内、桨外的点
W(os.path.join(PM, "system", "snappyHexMeshDict"), hdr("dictionary", "snappyHexMeshDict") + f"""
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
    paddle.stl {{ type triSurfaceMesh; name paddle; }}
{LINKS_GEOM}}}

castellatedMeshControls
{{
    maxLocalCells       2000000;
    maxGlobalCells      4000000;
    minRefinementCells  10;
    maxLoadUnbalance    0.1;
    nCellsBetweenLevels 3;
    features ( {{ file "paddle.eMesh"; level 1; }}{LINKS_FEAT} );
    refinementSurfaces
    {{
        paddle {{ level (1 1); patchInfo {{ type wall; inGroups (wall); }} }}
{LINKS_SURF}    }}
    resolveFeatureAngle 30;
    refinementRegions {{ }}
    locationInMesh {fmt(lim)};
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch    3;
    tolerance       2.0;
    nSolveIter      50;
    nRelaxIter      5;
    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    relativeSizes true; layers {{ }} expansionRatio 1.2; finalLayerThickness 0.4; minThickness 0.1; nGrow 0;
    featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
    maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-13; minTetQuality 1e-15;
    minArea -1; minTwist 0.02; minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1;
    nSmoothScale 4; errorReduction 0.75;
}}

writeFlags ( scalarLevels );
mergeTolerance 1e-6;
""")
W(os.path.join(PM, "system", "controlDict"), hdr("dictionary", "controlDict") + """
application     snappyHexMesh;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1;
deltaT          1;
writeControl    timeStep;
writeInterval   1;
writeFormat     binary;
writePrecision  8;
""")
W(os.path.join(PM, "system", "fvSchemes"), hdr("dictionary", "fvSchemes") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
W(os.path.join(PM, "system", "fvSolution"), hdr("dictionary", "fvSolution") + "\n")
W(os.path.join(PM, "system", "meshQualityDict"), hdr("dictionary", "meshQualityDict") + "\n")

# ---- 背景
W(os.path.join(CASE, "system", "blockMeshDict"), hdr("dictionary", "blockMeshDict") + f"""
scale 1;
vertices
(
    ({dom_min[0]:.4f} {dom_min[1]:.4f} {dom_min[2]:.4f})
    ({dom_max[0]:.4f} {dom_min[1]:.4f} {dom_min[2]:.4f})
    ({dom_max[0]:.4f} {dom_max[1]:.4f} {dom_min[2]:.4f})
    ({dom_min[0]:.4f} {dom_max[1]:.4f} {dom_min[2]:.4f})
    ({dom_min[0]:.4f} {dom_min[1]:.4f} {dom_max[2]:.4f})
    ({dom_max[0]:.4f} {dom_min[1]:.4f} {dom_max[2]:.4f})
    ({dom_max[0]:.4f} {dom_max[1]:.4f} {dom_max[2]:.4f})
    ({dom_min[0]:.4f} {dom_max[1]:.4f} {dom_max[2]:.4f})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
edges ();
boundary
(
    inlet  {{ type patch; faces ((1 2 6 5)); }}    // +x 端：来流从这里吹向 -x
    outlet {{ type patch; faces ((0 4 7 3)); }}    // -x 端
    sides  {{ type patch; faces ((0 1 5 4) (3 7 6 2)); }}
    bottom {{ type patch; faces ((0 3 2 1)); }}
    top    {{ type patch; faces ((4 5 6 7)); }}    // 刚盖近似自由液面（单相）
);
mergePatchPairs ();
""")
W(os.path.join(CASE, "system", "snappyHexMeshDict"), hdr("dictionary", "snappyHexMeshDict") + f"""
// 背景网格：只做区域加密（无表面），把桨扫掠区细化到 {a.bg*1000:.1f} mm
castellatedMesh true;
snap            false;
addLayers       false;
geometry
{{
    ref1 {{ type searchableBox; min {fmt(ref1[0])}; max {fmt(ref1[1])}; }}
    ref2 {{ type searchableBox; min {fmt(ref2[0])}; max {fmt(ref2[1])}; }}
}}
castellatedMeshControls
{{
    maxLocalCells 3000000; maxGlobalCells 6000000; minRefinementCells 10; maxLoadUnbalance 0.1;
    nCellsBetweenLevels 4;
    features ();
    refinementSurfaces {{ }}
    resolveFeatureAngle 30;
    refinementRegions
    {{
        ref1 {{ mode inside; levels ((1e15 1)); }}
        ref2 {{ mode inside; levels ((1e15 2)); }}
    }}
    locationInMesh ({dom_min[0]+0.05:.4f} {dom_min[1]+0.05:.4f} {dom_min[2]+0.05:.4f});
    allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; }}
addLayersControls {{ relativeSizes true; layers {{ }} expansionRatio 1.0; finalLayerThickness 0.3; minThickness 0.1; nGrow 0; featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50; }}
meshQualityControls {{ maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist 0.02; minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75; }}
writeFlags ();
mergeTolerance 1e-6;
""")
W(os.path.join(CASE, "system", "topoSetDict"), hdr("dictionary", "topoSetDict") + f"""
actions
(
    {{ name c0; type cellSet; action new; source regionToCell; insidePoints (({dom_min[0]+0.05:.4f} {dom_min[1]+0.05:.4f} {dom_min[2]+0.05:.4f})); }}
    {{ name c1; type cellSet; action new; source cellToCell; set c0; }}
    {{ name c1; type cellSet; action invert; }}
    {{ name c0; type cellZoneSet; action new; source setToCellZone; set c0; }}
    {{ name c1; type cellZoneSet; action new; source setToCellZone; set c1; }}
);
""")
W(os.path.join(CASE, "system", "setFieldsDict"), hdr("dictionary", "setFieldsDict") + """
defaultFieldValues ( volScalarFieldValue zoneID 123 );
regions
(
    cellToCell { set c0; fieldValues ( volScalarFieldValue zoneID 0 ); }
    cellToCell { set c1; fieldValues ( volScalarFieldValue zoneID 1 ); }
);
""")
W(os.path.join(CASE, "system", "controlDict"), hdr("dictionary", "controlDict") + f"""
application     overPimpleDyMFoam;
libs            ("liboverset.so" "libfvMotionSolvers.so");

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {endTime:.4f};
deltaT          {dt0};
writeControl    adjustableRunTime;
writeInterval   {writeInterval:.5f};
purgeWrite      0;
writeFormat     binary;
writePrecision  8;
writeCompression on;
timeFormat      general;
timePrecision   7;
runTimeModifiable true;
adjustTimeStep  yes;
maxCo           {a.maxCo};
maxDeltaT       {4*dt0};

functions
{{
    forces                    // 整条腿（桨 + 连杆段）的总力
    {{
        type            forces;
        libs            ("libforces.so");
        patches         ({"paddle links" if tri_links is not None else "paddle"});
        rho             rhoInf;
        rhoInf          {RHO};
        CofR            ({ld.O2[0]:.4f} 0 {ld.O2[1]:.4f});     // 力矩参考点 = O2 髋轴
        writeControl    timeStep;
        writeInterval   1;
        log             false;
    }}
    vorticity1
    {{
        type            vorticity;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
    }}
{LINKS_FO}}}
""")
divU = "Gauss limitedLinearV 1"
turbDiv = "" if a.turb == "laminar" else """
    div(phi,k)      Gauss limitedLinear 1;
    div(phi,omega)  Gauss limitedLinear 1;"""
W(os.path.join(CASE, "system", "fvSchemes"), hdr("dictionary", "fvSchemes") + f"""
ddtSchemes {{ default Euler; }}
gradSchemes {{ default Gauss linear; }}
divSchemes
{{
    default         none;
    div(phi,U)      {divU};{turbDiv}
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default corrected; }}
// 重叠网格洞切 (hole cutting)：默认 30³ 体素铺满整个背景域，分辩不了 5.6 mm 厚的桨板 → 泄漏成全洞。
// 把搜索盒限制在桨的运动包络内、体素 ~{vox*1000:.1f} mm。v2606 若报 searchBoxDivisions 类型错误，改成每个 zone 一组：((nx ny nz) (nx ny nz))
oversetInterpolation
{{
    method              inverseDistance;
    searchBox           {fmt(sb_min)} {fmt(sb_max)};
    searchBoxDivisions  ({sbd[0]} {sbd[1]} {sbd[2]});
}}
fluxRequired {{ default no; p; }}
wallDist {{ method meshWave; }}
""")
W(os.path.join(CASE, "system", "fvSolution"), hdr("dictionary", "fvSolution") + """
solvers
{
    cellDisplacement { solver PCG; preconditioner DIC; tolerance 1e-06; relTol 0; maxIter 100; }
    p       { solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.01; }
    pFinal  { $p; relTol 0; }
    "(U|k|omega)"       { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1; }
    "(U|k|omega)Final"  { $U; relTol 0; }
}
PIMPLE
{
    momentumPredictor   false;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 0;
    ddtCorr             true;
    correctPhi          false;
    oversetAdjustPhi    false;
    checkMeshCourantNo  yes;
}
relaxationFactors { equations { ".*" 1; } }
""")
def hier_n(n):
    """把 np 拆成 (nx ny nz)，x 方向优先（域沿 x 最长）"""
    best = (n, 1, 1)
    for nx in range(n, 0, -1):                       # 同样均匀时优先切 x（域沿 x 最长）
        for ny in range(n, 0, -1):
            if n % (nx * ny) == 0:
                nz = n // (nx * ny)
                if (max(nx, ny, nz), -nx) < (max(best), -best[0]): best = (nx, ny, nz)
    return best
hn = hier_n(a.np)
W(os.path.join(CASE, "system", "decomposeParDict"), hdr("dictionary", "decomposeParDict") + f"""
// 重叠网格并行：scotch 在某些分区下会触发 "Attempt to cast type patch to type lduInterface"（v2606，closed 桨算例复现），
// 改用 hierarchical 规则分区，沙箱与实机均正常。
numberOfSubdomains {a.np};
method          hierarchical;
hierarchicalCoeffs {{ n ({hn[0]} {hn[1]} {hn[2]}); order xyz; }}
""")

# ---- constant
W(os.path.join(CASE, "constant", "dynamicMeshDict"), hdr("dictionary", "dynamicMeshDict", "constant") + f"""
dynamicFvMesh       dynamicOversetFvMesh;
dynamicOversetFvMeshCoeffs {{ }}

solver          solidBody;
solidBodyCoeffs
{{
    cellZone        c1;
    solidBodyMotionFunction tabulated6DoFMotion;
    tabulated6DoFMotionCoeffs
    {{
        CofG            ({E0[0]:.6f} {Y_LEG:.6f} {E0[1]:.6f});
        timeDataFileName "<constant>/6DoF.dat";
    }}
}}
""")
with open(os.path.join(CASE, "constant", "6DoF.dat"), "w") as f:
    f.write(f"{nT}\n(\n")
    for k in range(nT):
        f.write(f"({t[k]:.6f} (({trans[k,0]:.6f} 0 {trans[k,1]:.6f}) ({rx[k]:.5f} {ry[k]:.5f} {rz[k]:.5f})))\n")
    f.write(")\n")
W(os.path.join(CASE, "constant", "transportProperties"), hdr("dictionary", "transportProperties", "constant") + f"""
transportModel  Newtonian;
nu              {NU};
""")
if a.turb == "laminar":
    W(os.path.join(CASE, "constant", "turbulenceProperties"), hdr("dictionary", "turbulenceProperties", "constant") + "\nsimulationType laminar;\n")
else:
    W(os.path.join(CASE, "constant", "turbulenceProperties"), hdr("dictionary", "turbulenceProperties", "constant") + """
simulationType RAS;
RAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }
""")

# ---- 0.orig
Uin = f"(-{a.U:.4f} 0 0)" if a.U > 0 else "(0 0 0)"
W(os.path.join(CASE, "0.orig", "U"), hdr("volVectorField", "U", "0") + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform {Uin};
boundaryField
{{
    oversetPatch {{ type overset; }}
    {WALLS} {{ type movingWallVelocity; value uniform (0 0 0); }}
    inlet        {{ type fixedValue; value uniform {Uin}; }}
    outlet       {{ type inletOutlet; inletValue uniform (0 0 0); value uniform {Uin}; }}
    "(sides|bottom|top)" {{ type slip; }}
}}
""")
W(os.path.join(CASE, "0.orig", "p"), hdr("volScalarField", "p", "0") + """
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    oversetPatch { type overset; }
    """ + WALLS + """ { type zeroGradient; }
    inlet        { type zeroGradient; }
    outlet       { type fixedValue; value uniform 0; }
    "(sides|bottom|top)" { type zeroGradient; }
}
""")
W(os.path.join(CASE, "0.orig", "zoneID"), hdr("volScalarField", "zoneID", "0") + """
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    oversetPatch { type overset; value uniform 0; }
    ".*"         { type zeroGradient; }
}
""")
if a.turb != "laminar":
    k_in = 1.5 * (0.02 * max(a.U, 0.1))**2; om_in = np.sqrt(k_in) / (0.09**0.25 * 0.07 * 0.03)
    W(os.path.join(CASE, "0.orig", "k"), hdr("volScalarField", "k", "0") + f"""
dimensions [0 2 -2 0 0 0 0];
internalField uniform {k_in:.3e};
boundaryField
{{
    oversetPatch {{ type overset; }}
    {WALLS} {{ type kqRWallFunction; value uniform {k_in:.3e}; }}
    inlet {{ type fixedValue; value uniform {k_in:.3e}; }}
    outlet {{ type inletOutlet; inletValue uniform {k_in:.3e}; value uniform {k_in:.3e}; }}
    "(sides|bottom|top)" {{ type slip; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "omega"), hdr("volScalarField", "omega", "0") + f"""
dimensions [0 0 -1 0 0 0 0];
internalField uniform {om_in:.3e};
boundaryField
{{
    oversetPatch {{ type overset; }}
    {WALLS} {{ type omegaWallFunction; value uniform {om_in:.3e}; }}
    inlet {{ type fixedValue; value uniform {om_in:.3e}; }}
    outlet {{ type inletOutlet; inletValue uniform {om_in:.3e}; value uniform {om_in:.3e}; }}
    "(sides|bottom|top)" {{ type slip; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "nut"), hdr("volScalarField", "nut", "0") + """
dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{
    oversetPatch { type overset; }
    """ + WALLS + """ { type nutkWallFunction; value uniform 0; }
    ".*" { type calculated; value uniform 0; }
}
""")

# ---- Allrun
W(os.path.join(CASE, "Allrun"), f"""#!/bin/sh
cd "${{0%/*}}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 1 部件网格（桨板，贴体）
( cd paddleMesh || exit 1
  runApplication blockMesh
  runApplication surfaceFeatureExtract
  runApplication snappyHexMesh -overwrite
  runApplication checkMesh )

# 2 背景网格（只做区域加密）
runApplication blockMesh
runApplication snappyHexMesh -overwrite

# 3 合并 → 分区 → zoneID
runApplication mergeMeshes . paddleMesh -overwrite
runApplication topoSet
restore0Dir
runApplication setFields
runApplication -s merged checkMesh

# 4 求解
runApplication decomposePar
runParallel $(getApplication)
runParallel -s reconstruct redistributePar -reconstruct
touch case.foam
""")
os.chmod(os.path.join(CASE, "Allrun"), 0o755)
W(os.path.join(CASE, "Allclean"), """#!/bin/sh
cd "${0%/*}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/CleanFunctions"
cleanCase
rm -rf paddleMesh/constant/polyMesh paddleMesh/constant/extendedFeatureEdgeMesh paddleMesh/constant/triSurface/*.eMesh paddleMesh/log.* postProcessing 0
""")
os.chmod(os.path.join(CASE, "Allclean"), 0o755)

# ---- 记录 & 校核数据
gj_extra = {}
if TRAJ is not None:   # 两态合成的切换相位：开度过 0.5 的时刻（τ_open 为负 = 上一周期末段就开始张开）
    f1 = TRAJ["fan"]; tt1 = TRAJ["t"]
    down = np.where((f1[:-1] >= 0.5) & (f1[1:] < 0.5))[0]; up = np.where((f1[:-1] < 0.5) & (f1[1:] >= 0.5))[0]
    tau_c = float(tt1[down[0] + 1] / T) if len(down) else 1.0
    tau_o = (float(tt1[up[-1] + 1] / T) if len(up) else 0.0); tau_o = tau_o - 1.0 if tau_o > 0.5 else tau_o
    gj_extra = dict(traj=os.path.abspath(a.traj), TAU_CLOSE=tau_c, TAU_OPEN=tau_o)
if FEATHER is not None: gj_extra.update(FEATHER=FEATHER, PLATE_W=a.plate_w, PLATE_L=(a.plate_l or ld.HEAD))
elif a.fan == "plate": gj_extra.update(PLATE_W=a.plate_w, PLATE_L=(a.plate_l or ld.HEAD))
json.dump(dict(fan=a.fan, U=a.U, cycles=a.cycles, T=T, DUTY=(gj_extra["TAU_CLOSE"] if TRAJ is not None else ld.DUTY), SWEEP=ld.SWEEP, PHI_EXT=ld.PHI_EXT, PHI_RET=ld.PHI_RET, **gj_extra,
               FAN_DEG=ld.FAN_DEG, FAN_CLOSE=ld.FAN_CLOSE, CLOSED_W=ld.CLOSED_W, SWITCH_FRAC=getattr(ld, "SWITCH_FRAC", 0.3), VEL_LAW=getattr(ld, "VEL_LAW", "cos"), RAMP=getattr(ld, "RAMP", 0.15), turb=a.turb, bg=a.bg, comp=a.comp,
               E0=E0.tolist(), psi0_deg=float(np.degrees(psi0)), dt0=dt0, tip_speed_max=float(tip_speed), leg=bool(tri_links is not None), WL=a.wl, **leg_info),
          open(os.path.join(CASE, "gait.json"), "w"), indent=1, ensure_ascii=False)
np.savetxt(os.path.join(CASE, "kinematics_check.csv"), np.column_stack([t, E, tip, np.degrees(psi), np.degrees(phi), power, beta, rx, ry, rz]),
           delimiter=",", header="t,Ex,Ez,tipx,tipz,psi_deg,phi_deg,power,beta_deg,rx,ry,rz", comments="", fmt="%.6f")

print(f"算例：{CASE}   脚蹼 {a.fan}   U = {a.U} m/s   {a.cycles} 个周期 (T = {T} s) → endTime {endTime:.3f} s")
print(f"桨板 STL：{len(tri)} 三角形，世界包围盒 x[{pmin[0]:.3f},{pmax[0]:.3f}] y[{pmin[1]:.3f},{pmax[1]:.3f}] z[{pmin[2]:.3f},{pmax[2]:.3f}]")
print(f"部件网格盒 {n_s}×{n_w}×{n_n} = {n_s*n_w*n_n} 单元（{a.comp*1000:.2f} mm，表面再细化 1 级）" + ("  含连杆段 links.stl（patch links）" if tri_links is not None else ""))
print(f"部件运动包络 x[{env_min[0]:.3f},{env_max[0]:.3f}] y[{env_min[1]:.3f},{env_max[1]:.3f}] z[{env_min[2]:.3f},{env_max[2]:.3f}]")
print(f"背景域 x[{dom_min[0]:.3f},{dom_max[0]:.3f}] y[{dom_min[1]:.3f},{dom_max[1]:.3f}] z[{dom_min[2]:.3f},{dom_max[2]:.3f}]  blockMesh {nx}×{ny}×{nz}，加密到 {a.bg*1000:.1f} mm")
print(f"桨尖最大速度 {tip_speed:.2f} m/s → 初始 Δt = {dt0} s，maxCo {a.maxCo}；每周期写 {a.writes} 帧")
print(f"洞切搜索盒 {fmt(sb_min)}–{fmt(sb_max)}，体素 {sbd[0]}×{sbd[1]}×{sbd[2]}（≈{vox*1000:.1f} mm）")
print(f"运动表 6DoF.dat：{nT} 行，平动 = E(t)−E0，绕 y 转角 = −(ψ−ψ0)  [CofG = E0 = ({E0[0]:.4f}, 0, {E0[1]:.4f})]")
print("跑：cd", CASE, "&& ./Allrun")
