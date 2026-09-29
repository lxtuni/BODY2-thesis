#!/usr/bin/env python3
"""
整机（船体 + 四条腿对角步态 + 自由液面）overset 算例生成器 —— Stage 3b
    python3 make_full_cfd.py --case ~/run/full_open_U01  --fan open   --U 0.1
    python3 make_full_cfd.py --case ~/run/full_closed_U01 --fan closed --U 0.1
    cd ~/run/full_open_U01 && ./Allrun

内容：
  · 船体 hull_body_only.stl（密封身体 + 盖板）固定在背景网格里，水线 z = --wl（真机船体中部 −0.010），k-ω SST
  · 四块脚蹼各是一个 overset 部件（paddleMesh_FL/FR/RL/RR），运动由 leg_dynamics.py 的连杆运动学给出，
    对角步态：FL+RR 相位 0，FR+RL 相位 0.5（--gait 可改）；后腿 = 前腿机构沿 x 平移 --dx-rear
  · 两态脚蹼仍用 open / closed 两个算例按每条腿自己的相位拼接（full_forces.py）
  · 自由液面：overInterDyMFoam（VOF），来流 U 从 +x 吹向 −x（机器人固定，拖曳水池等效）
  · 受力：forces_hull / forces_FL / forces_FR / forces_RL / forces_RR 分开输出（moment 参考点 = 船体中点）
选项：
  --no-hull      不放船体（Stage 3a：四腿干扰）        --no-fs   单相 overPimpleDyMFoam（只能和 --no-hull 一起用）
  --legs FL,RR   只放部分腿                             --cycles 2  周期数（第 1 周期为启动瞬态）
  --bg 0.004 / --comp 0.0025 / --base 0.016             背景加密尺寸 / 部件基准尺寸 / 背景 blockMesh 基准
坐标：机体坐标系，x 前进（艏 +x）、z 向上、y 左。FL 腿 O2 = (0.0634, +0.0688, 0.004)。
"""
import argparse, os, sys, json, struct, shutil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True)
ap.add_argument("--fan", choices=["open", "closed"], required=True)
ap.add_argument("--U", type=float, default=0.0)
ap.add_argument("--cycles", type=float, default=2.0)
ap.add_argument("--bg", type=float, default=0.004, help="背景网格在桨扫掠区 / 自由面 / 船体表面的单元尺寸 [m]")
ap.add_argument("--comp", type=float, default=0.0025, help="部件网格基准单元尺寸 [m]（桨面再细化 1 级）")
ap.add_argument("--base", type=float, default=None, help="背景 blockMesh 基准单元 [m]，默认 4×bg")
ap.add_argument("--turb", default="kOmegaSST", choices=["laminar", "kOmegaSST"])
ap.add_argument("--np", type=int, default=8)
ap.add_argument("--nsub", type=int, default=800)
ap.add_argument("--maxCo", type=float, default=1.0)
ap.add_argument("--maxAlphaCo", type=float, default=1.0, help="界面 Courant 上限（alpha 有 2 个子循环，1.0 等效每子循环 0.5）")
ap.add_argument("--limitU", type=float, default=1.2, help="fvOptions limitVelocity 的速度上限 [m/s]，压住空气区/薄片单元里的伪速度免得 Δt 被一个单元卡死；0 = 不加")
ap.add_argument("--overset", default="inverseDistance", choices=["inverseDistance", "trackingInverseDistance"], help="洞切/供体方法；tracking 版（v2006+）增量更新供体，mesh.update 快数倍")
ap.add_argument("--writes", type=int, default=20, help="每周期写出帧数")
ap.add_argument("--vox", type=float, default=None, help="洞切体素 [m]（默认 min(桨厚/4.5, comp/2) ≈ 1.2 mm）")
ap.add_argument("--wl", type=float, default=-0.010, help="水线 z [m]")
ap.add_argument("--hull", default=None, help="船体 STL（默认 ../hull_body_only.stl）")
ap.add_argument("--no-hull", action="store_true")
ap.add_argument("--no-fs", action="store_true", help="单相（无自由面），只能与 --no-hull 同用")
ap.add_argument("--legs", default="FL,FR,RL,RR")
ap.add_argument("--gait", default="diag", choices=["diag", "sync", "bound", "walk"],
                help="diag: FL+RR=0, FR+RL=0.5 | sync: 全 0 | bound: 前 0 后 0.5 | walk: FL 0, RR .25, FR .5, RL .75")
ap.add_argument("--dx-rear", type=float, default=-0.1435, help="后腿机构相对前腿的 x 平移 [m]（STEP 量得）")
ap.add_argument("--y-front", type=float, default=0.075, help="前腿桨中心 |y| [m]。CAD 是 0.0688，但那样收拢桨离船体只有 7.9 mm、张开扇会碰船底棱，overset 供体链会断；默认外移到 0.075")
ap.add_argument("--y-rear", type=float, default=0.0700)
ap.add_argument("--min-gap", type=float, default=None, help="桨与船体的最小允许间隙 [m]（默认 2.5×comp + 2×bg ≈ 14 mm）；张开扇内侧超出的部分自动裁掉")
ap.add_argument("--no-clip", action="store_true", help="不自动裁扇：间隙不够就直接退出")
ap.add_argument("--hull-level", type=int, default=3, help="船体表面 snappy 加密级（3 = bg/2 = 2 mm）")
ap.add_argument("--paddle", default="cad", choices=["cad", "model"], help="收拢态桨：cad = STEP 零件 扇形足 的精确平面轮廓与板厚（paddle_closed_cad.json）；model = 理想化轮廓")
ap.add_argument("--thk-open", type=float, default=None, help="张开扇板厚 [m]（默认 = CAD 桨头板厚 2.7 mm；model 模式 5.6 mm）")
ap.add_argument("--paddle-level", type=int, default=1, help="部件网格里桨面的 snappy 加密级（1 = comp/2 = 1.25 mm；板厚 2.7 mm 时 2 层单元，若 checkMesh 不过改 2）")
a = ap.parse_args()
if a.no_fs and not a.no_hull: sys.exit("--no-fs（单相）只能和 --no-hull 一起用：船体一半在水上，单相算不了")

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import leg_dynamics as ld
CASE = os.path.abspath(os.path.expanduser(a.case)); T = ld.T_PER; RHO, NU = ld.RHO, ld.NU
WL = a.wl; FS = not a.no_fs
HULL = a.hull or os.path.join(HERE, "..", "hull_body_only.stl")
if not a.no_hull and not os.path.exists(HULL): sys.exit(f"找不到船体 STL：{HULL}（robot_real/hull_body_only.stl）")
PHASE = {"diag": dict(FL=0.0, RR=0.0, FR=0.5, RL=0.5), "sync": dict(FL=0, FR=0, RL=0, RR=0),
         "bound": dict(FL=0, FR=0, RL=0.5, RR=0.5), "walk": dict(FL=0.0, RR=0.25, FR=0.5, RL=0.75)}[a.gait]
LEGS = [s.strip().upper() for s in a.legs.split(",") if s.strip()]
for L in LEGS:
    if L not in ("FL", "FR", "RL", "RR"): sys.exit(f"腿代号只能是 FL/FR/RL/RR：{L}")
def leg_geom(L):
    dx = 0.0 if L[0] == "F" else a.dx_rear
    y = (a.y_front if L[0] == "F" else a.y_rear) * (1 if L[1] == "L" else -1)
    return dx, y

# ------------------------------------------------------------------ 运动学（每条腿：相位偏移 + x 平移）
nT = int(a.nsub * a.cycles) + 1
t = np.linspace(0, a.cycles * T, nT)
legs = {}
for L in LEGS:
    dx, y = leg_geom(L); ph = PHASE[L]
    psi, phi, power = ld.gait(t + ph * T)
    P = [ld.pose(p, q) for p, q in zip(psi, phi)]
    E = np.array([p["E"] for p in P]) + np.array([dx, 0.0]); u = np.array([p["u"] for p in P])
    n = np.stack([-u[:, 1], u[:, 0]], 1)
    legs[L] = dict(dx=dx, y=y, ph=ph, psi=psi, phi=phi, E=E, u=u, n=n, E0=E[0], u0=u[0], n0=n[0], psi0=psi[0],
                   tip=E + (ld.STEM + ld.HEAD) * u, O2=ld.O2 + np.array([dx, 0.0]))
tip_speed = max(np.max(np.linalg.norm(np.gradient(l["tip"], t, axis=0), axis=1)) for l in legs.values())

# ------------------------------------------------------------------ 桨板几何（同 make_leg_cfd.py）
def paddle_outline(fan_open):
    S, H, SW, HW = ld.STEM, ld.HEAD, ld.STEM_W, ld.HEAD_W
    if fan_open:
        ang = np.deg2rad(ld.FAN_DEG); th = np.linspace(ang / 2, -ang / 2, 48)
        fan = [(ld.FAN_R * np.sin(x), S + ld.FAN_R * np.cos(x)) for x in th]
        return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + fan + [(-SW / 2, S)]
    k = np.linspace(0, 1, 24)
    sh = [(HW / 2 * np.sin(np.pi / 2 * q), S + 0.025 * (1 - np.cos(np.pi / 2 * q))) for q in k]
    right = [(w, z) for w, z in sh if w > SW / 2]
    return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + right + [(HW / 2, S + H), (-HW / 2, S + H)] + [(-w, z) for w, z in right[::-1]] + [(-SW / 2, S)]

def extrude_world(outline, thk, E0, u0, n0, Y, nrm_c=0.0):
    """把 (lat, s) 轮廓沿桨法向挤出 thk（中心在 nrm_c）并放到世界坐标"""
    Pl = np.array(outline, dtype=float)
    area2 = np.sum(Pl[:, 0] * np.roll(Pl[:, 1], -1) - np.roll(Pl[:, 0], -1) * Pl[:, 1])
    if area2 < 0: Pl = Pl[::-1]
    c = Pl.mean(0); N = len(Pl)
    u0w = np.array([u0[0], 0.0, u0[1]]); n0w = np.array([n0[0], 0.0, n0[1]]); yw = np.array([0.0, 1.0, 0.0])
    def Wp(lat, s, nrm): nrm = nrm + nrm_c; return E0[0] + s * u0[0] + nrm * n0[0], lat + Y, E0[1] + s * u0[1] + nrm * n0[1]
    tris, want = [], []
    for i in range(N):
        (la, sa), (lb, sb) = Pl[i], Pl[(i + 1) % N]
        tris.append([Wp(c[0], c[1], thk / 2), Wp(la, sa, thk / 2), Wp(lb, sb, thk / 2)]); want.append(n0w)
        tris.append([Wp(c[0], c[1], -thk / 2), Wp(la, sa, -thk / 2), Wp(lb, sb, -thk / 2)]); want.append(-n0w)
        dlat, ds = lb - la, sb - sa; out = ds * yw - dlat * u0w
        p1, p2 = Wp(la, sa, thk / 2), Wp(lb, sb, thk / 2); q1, q2 = Wp(la, sa, -thk / 2), Wp(lb, sb, -thk / 2)
        tris.append([p1, q1, q2]); want.append(out); tris.append([p1, q2, p2]); want.append(out)
    tri = np.array(tris); want = np.array(want)
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); flip = np.einsum("ij,ij->i", nrm, want) < 0
    tri[flip] = tri[flip][:, ::-1]
    return tri

def write_stl(tri, fn):
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
    rec = np.zeros(len(tri), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])); rec["n"] = nrm; rec["v"] = tri
    open(fn, "wb").write(b"\0" * 80 + struct.pack("<I", len(tri)) + rec.tobytes())

# overset 供体链要求的最小桨–船间隙：船体在部件里切出的洞层(≈体素+部件单元) + 部件前沿单元 + 其供体背景单元不能是桨洞(体素+背景单元)
vox0 = a.vox if a.vox else min(ld.THK / 4.5, a.comp / 2)          # ≈1.24 mm；洞判定按单元包围盒与体素重叠，与板厚无关
MIN_GAP = a.min_gap if a.min_gap else max(0.010, 2 * vox0 + 2 * (a.comp / 2) + (a.bg / 2) + 0.002)
fmt = lambda v: "(" + " ".join(f"{x:.6g}" for x in v) + ")"
def paddle_outline_clipped(fan_open, lat_min, lat_max):
    """桨面轮廓，横向限制在 [lat_min, lat_max]（内侧靠船体的一边可能被裁掉）"""
    pts = paddle_outline(fan_open)
    if fan_open:
        S = ld.STEM; out = []
        for (w, z) in pts:
            if w < lat_min or w > lat_max:
                # 扇边界上的点被裁：投到裁切线上（与扇弧的交点由相邻点插值近似）
                wc = min(max(w, lat_min), lat_max)
                zc = S + np.sqrt(max(ld.FAN_R**2 - wc**2, 0.0))
                out.append((wc, zc))
            else: out.append((w, z))
        # 去掉裁切线上重复的点，并补上裁切线与扇形侧边（顶点→弧端的半径线）的交点，使裁掉的是一整条竖直边
        res = []
        for p in out:
            if res and abs(p[0] - res[-1][0]) < 1e-9 and abs(p[1] - res[-1][1]) < 1e-9: continue
            res.append(p)
        ang = np.deg2rad(ld.FAN_DEG); R = ld.FAN_R; hw0 = R * np.sin(ang / 2)
        def side_pt(latc):                      # 半径线 (0,S)→(±hw0, S+R cos(ang/2)) 上 lat = latc 的点
            k = abs(latc) / hw0; return (latc, S + k * R * np.cos(ang / 2))
        fixed = []
        for i, p in enumerate(res):
            if lat_max < hw0 - 1e-9 and abs(p[0] - lat_max) < 1e-9 and i > 0 and res[i - 1][1] <= S + 1e-9:
                fixed.append(side_pt(lat_max))          # 外侧裁切：先从半径线交点上到弧
            fixed.append(p)
            if lat_min > -hw0 + 1e-9 and abs(p[0] - lat_min) < 1e-9 and i + 1 < len(res) and res[i + 1][1] <= S + 1e-9:
                fixed.append(side_pt(lat_min))          # 内侧裁切：从弧下到半径线交点再回顶点
        return fixed
    return [(min(max(w, lat_min), lat_max), z) for (w, z) in pts]

# ---- 船体表面点云（间隙检查用）
hull_tree = None
if not a.no_hull:
    b_ = open(HULL, "rb").read()
    if b_[:5] == b"solid" and b"facet" in b_[:400]:
        TT = np.array([[float(x) for x in ln.split()[1:4]] for ln in b_.decode(errors="ignore").splitlines() if ln.strip().startswith("vertex")]).reshape(-1, 3, 3)
    else:
        nh_ = struct.unpack("<I", b_[80:84])[0]
        TT = np.frombuffer(b_[84:84 + nh_ * 50], dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]))["v"].reshape(-1, 3, 3).astype(float)
    uu, vv = np.meshgrid(np.linspace(0, 1, 4), np.linspace(0, 1, 4)); uu = uu.ravel(); vv = vv.ravel(); k_ = uu + vv <= 1
    uu, vv = uu[k_], vv[k_]
    pc = (TT[:, None, 0] + uu[None, :, None] * (TT[:, None, 1] - TT[:, None, 0]) + vv[None, :, None] * (TT[:, None, 2] - TT[:, None, 0])).reshape(-1, 3)
    try:
        from scipy.spatial import cKDTree
        hull_tree = cKDTree(pc)
    except Exception:
        class _Brute:                      # 没有 scipy 时的退路：分块 numpy 最近点（船体点云抽稀到 ~2 万点）
            def __init__(self, P):
                self.P = P[:: max(1, len(P) // 20000)]
            def query(self, S):
                S = np.asarray(S); d = np.empty(len(S))
                for i in range(0, len(S), 256):
                    blk = S[i:i + 256]; dd = np.sqrt(((blk[:, None, :] - self.P[None, :, :]) ** 2).sum(-1)); d[i:i + 256] = dd.min(1)
                return d, None
        hull_tree = _Brute(pc); print("（没有 scipy，间隙检查用 numpy 逐块计算，慢一点但结果一样）")

# ---- 收拢态桨的精确几何（STEP 零件 扇形足：平面轮廓 + 桨头板厚/偏置），文件由 robot_real/parts/扇形足_1.stl 抽出
CAD = None
if a.paddle == "cad":
    cad_fn = os.path.join(HERE, "paddle_closed_cad.json")
    if os.path.exists(cad_fn):
        CAD = json.load(open(cad_fn))
        CAD["outline"] = [tuple(p) for p in CAD["outline"]]
    else:
        print("!! 找不到 paddle_closed_cad.json，收拢态退回理想化轮廓（--paddle model）"); a.paddle = "model"
if CAD is not None:
    HEAD_NRM = tuple(CAD["head_nrm"])                       # 桨头板 nrm 范围 [-3.04, -0.35] mm（板厚 2.7，偏在 -n 侧）
    STEM_NRM = tuple(CAD["stem_nrm"])                       # 杆 nrm 范围 [-3.0, +1.5] mm
    THK_HEAD = HEAD_NRM[1] - HEAD_NRM[0]; NRM_C = 0.5 * (HEAD_NRM[0] + HEAD_NRM[1])
else:
    HEAD_NRM = (-ld.THK / 2, ld.THK / 2); STEM_NRM = HEAD_NRM; THK_HEAD = ld.THK; NRM_C = 0.0
THK_OPEN = a.thk_open if a.thk_open else THK_HEAD           # 张开扇板厚：默认与 CAD 桨头板一样
def closed_outline():
    return list(CAD["outline"]) if CAD is not None else paddle_outline(False)
def lat_range_of(poly, s_):
    """多边形 (lat, s) 与直线 s = s_ 的交点范围"""
    xs = []
    n_ = len(poly)
    for i in range(n_):
        (l1, s1), (l2, s2) = poly[i], poly[(i + 1) % n_]
        if (s1 - s_) * (s2 - s_) <= 0 and s1 != s2: xs.append(l1 + (l2 - l1) * (s_ - s1) / (s2 - s1))
    return (min(xs), max(xs)) if xs else None
def sample_plate(poly, nrm_lo, nrm_hi, n_s=30, n_lat=7):
    """桨板两个表面的采样点（局部 (s, lat, nrm)）"""
    ss = [p[1] for p in poly]; pts = []
    for s_ in np.linspace(min(ss) + 1e-4, max(ss) - 1e-4, n_s):
        r = lat_range_of(poly, s_)
        if r is None: continue
        for la in np.linspace(r[0], r[1], n_lat):
            for off in (nrm_lo, nrm_hi): pts.append((s_, la, off))
    return np.array(pts)

def clearance(l, poly):
    """一条腿整个周期内桨板（轮廓 poly，含板厚）与船体的最小距离；返回 (距离, tau, 点)。poly 的 lat 已是该腿的局部横向坐标"""
    if hull_tree is None: return (1.0, 0.0, None)
    best = (1e9, 0.0, None)
    if "E1" not in l:                                          # 单独采一个完整周期（与 --cycles 无关）
        t1 = np.linspace(0, T, 400); ps1, ph1, _ = ld.gait(t1 + l["ph"] * T)
        P1 = [ld.pose(p, q) for p, q in zip(ps1, ph1)]
        l["E1"] = np.array([p["E"] for p in P1]) + np.array([l["dx"], 0.0]); l["u1"] = np.array([p["u"] for p in P1]); l["t1"] = t1
    loc = sample_plate(poly, min(STEM_NRM[0], HEAD_NRM[0]), max(STEM_NRM[1], HEAD_NRM[1]))
    for k in range(len(l["t1"])):
        E, u = l["E1"][k], l["u1"][k]; nn_ = np.array([-u[1], u[0]])
        S = np.column_stack([E[0] + loc[:, 0] * u[0] + loc[:, 2] * nn_[0], l["y"] + loc[:, 1], E[1] + loc[:, 0] * u[1] + loc[:, 2] * nn_[1]])
        d, _ = hull_tree.query(S); i = int(d.argmin())
        if d[i] < best[0]: best = (float(d[i]), float(l["t1"][k] / T), S[i])
    return best

half_w_full = max(abs(p[0]) for p in (paddle_outline(True) if a.fan == "open" else closed_outline()))
clip_info = {}
for L, l in legs.items():
    sgn_in = -1.0 if l["y"] > 0 else 1.0
    hw = half_w_full
    if a.fan == "closed":
        l["outline"] = closed_outline(); g_full = clearance(l, l["outline"]); hw_in = hw
        if hull_tree is not None and g_full[0] < MIN_GAP:
            sys.exit(f"!! {L} 腿 closed 态桨与船体最小间隙只有 {g_full[0]*1000:.1f} mm（tau={g_full[1]:.2f}，点 {np.round(g_full[2],3)}），"
                     f"小于 overset 需要的 {MIN_GAP*1000:.0f} mm。办法：--y-front/--y-rear 把腿向外移；--comp/--bg 更细；或减小 SWEEP/PHI_RET")
    else:
        def opn(hw_in_):
            lat_min, lat_max = (-hw_in_, hw) if sgn_in < 0 else (-hw, hw_in_)     # 局部 lat 坐标 = y - y_leg；内侧 = 朝船体
            return paddle_outline_clipped(True, lat_min, lat_max)
        g_full = clearance(l, opn(hw)); hw_in = hw
        if hull_tree is not None and g_full[0] < MIN_GAP:
            if a.no_clip:
                sys.exit(f"!! {L} 腿 open 态扇与船体最小间隙只有 {g_full[0]*1000:.1f} mm（tau={g_full[1]:.2f}），小于 {MIN_GAP*1000:.0f} mm（--no-clip 已禁止裁扇）")
            for hw_try in np.arange(hw, 0.004, -0.002):                  # 逐步裁内侧半宽直到满足间隙
                if clearance(l, opn(hw_try))[0] >= MIN_GAP: hw_in = hw_try; break
            else:
                sys.exit(f"!! {L} 腿：即使把张开扇内侧裁到 4 mm 也离船体不够 {MIN_GAP*1000:.0f} mm（桨杆本身太近），需要移腿或改步态")
            clip_info[L] = dict(hw_in=float(hw_in), hw_out=float(hw), gap_full=float(g_full[0]), tau=float(g_full[1]))
            print(f"!! {L} 腿：张开扇整周期与船体最小间隙 {g_full[0]*1000:.1f} mm < {MIN_GAP*1000:.0f} mm（tau={g_full[1]:.2f}），"
                  f"内侧半宽裁到 {hw_in*1000:.0f} mm（外侧仍 {hw*1000:.0f} mm）→ 面积 ×{(hw_in+hw)/(2*hw):.2f}")
        l["outline"] = opn(hw_in)
    l["gap"] = clearance(l, l["outline"])
    print(f"{L}: {a.fan} 态桨与船体最小间隙 {l['gap'][0]*1000:.1f} mm（tau={l['gap'][1]:.2f}，含板厚）" if hull_tree is not None else f"{L}: 无船体")
outline = None
half_w = half_w_full
# 部件盒对桨的余量：必须 ≳ 2.5 个背景单元——overset 的 walkFront 只给洞周围一层单元找供体（donor），
# 那层单元的中心若落在部件盒外就没有供体，会被判成洞并连锁把整个背景吃掉（粗网格测试里复现过）
m_ax = m_n = max(0.016, 3 * a.bg + 0.004); m_lat = max(0.012, 3 * a.bg + 0.004)
s_max_pad = max(max(p[1] for p in l["outline"]) for l in legs.values())
s_lo, s_hi = -m_ax, s_max_pad + m_ax
nrm_lo = min(STEM_NRM[0], HEAD_NRM[0], NRM_C - THK_OPEN / 2); nrm_hi = max(STEM_NRM[1], HEAD_NRM[1], NRM_C + THK_OPEN / 2)
n_lo_box, n_hi = nrm_lo - m_n, nrm_hi + m_n; w_hi = half_w + m_lat
env_min_all = np.full(3, np.inf); env_max_all = np.full(3, -np.inf)
for L, l in legs.items():
    E0, u0, n0, Y = l["E0"], l["u0"], l["n0"], l["y"]
    if a.fan == "closed":
        # 桨头板：CAD 平面轮廓 × 板厚（偏置 NRM_C）；杆：3 mm 宽 × STEM_NRM 的矩形杆，从 E 后 2.5 mm 到桨头根部（与板重叠成一体）
        head = extrude_world(l["outline"], THK_HEAD, E0, u0, n0, Y, nrm_c=NRM_C)
        s_head0 = CAD["s_head0"] if CAD is not None else ld.STEM
        stem = extrude_world([(-0.0015, -0.0025), (0.0015, -0.0025), (0.0015, s_head0 + 0.004), (-0.0015, s_head0 + 0.004)],
                             STEM_NRM[1] - STEM_NRM[0], E0, u0, n0, Y, nrm_c=0.5 * (STEM_NRM[0] + STEM_NRM[1]))
        l["tri"] = np.concatenate([head, stem])
    else:
        plate = extrude_world(l["outline"], THK_OPEN, E0, u0, n0, Y, nrm_c=NRM_C)
        stem = extrude_world([(-0.0015, -0.0025), (0.0015, -0.0025), (0.0015, ld.STEM + 0.004), (-0.0015, ld.STEM + 0.004)],
                             STEM_NRM[1] - STEM_NRM[0], E0, u0, n0, Y, nrm_c=0.5 * (STEM_NRM[0] + STEM_NRM[1]))
        l["tri"] = np.concatenate([plate, stem])
    def box_pt(s, nn, lat, E0=E0, u0=u0, n0=n0, Y=Y):
        p = E0 + s * u0 + nn * n0; return (p[0], lat + Y, p[1])
    V = [box_pt(s_lo, n_lo_box, -w_hi), box_pt(s_hi, n_lo_box, -w_hi), box_pt(s_hi, n_lo_box, w_hi), box_pt(s_lo, n_lo_box, w_hi),
         box_pt(s_lo, n_hi, -w_hi), box_pt(s_hi, n_hi, -w_hi), box_pt(s_hi, n_hi, w_hi), box_pt(s_lo, n_hi, w_hi)]
    v0, v1, v3, v4 = map(np.array, (V[0], V[1], V[3], V[4]))
    if np.dot(np.cross(v1 - v0, v3 - v0), v4 - v0) < 0: V = V[4:] + V[:4]
    l["V"] = V; l["lim"] = box_pt(s_hi - 0.004, n_hi - 0.004, w_hi - 0.004)
    # 部件盒内靠船体一侧 15 mm 的条带再加密一级（1.25 mm）：桨–船间隙里的供体单元要够细
    if not a.no_hull:
        Vb = np.array(V); sgn_in = -1.0 if Y > 0 else 1.0
        y_in = Y - sgn_in * (-w_hi) if False else (Y - w_hi if sgn_in < 0 else Y + w_hi)     # 盒子内侧面的 y
        ylo, yhi = (y_in, y_in + 0.015) if sgn_in < 0 else (y_in - 0.015, y_in)
        bmin = np.array([Vb[:, 0].min(), ylo, Vb[:, 2].min()]); bmax = np.array([Vb[:, 0].max(), yhi, Vb[:, 2].max()])
        l["strip_geom"] = f"strip_{L} {{ type searchableBox; min {fmt(bmin)}; max {fmt(bmax)}; }}"
        l["strip"] = f"strip_{L} {{ mode inside; levels ((1e15 1)); }}"
    else:
        l["strip_geom"] = ""; l["strip"] = ""
    # 运动包络
    corners = np.array([(s, nn, lat) for s in (s_lo, s_hi) for nn in (n_lo_box, n_hi) for lat in (-w_hi, w_hi)])
    emin = np.full(3, np.inf); emax = np.full(3, -np.inf)
    for k in range(0, nT, max(1, nT // 400)):
        uk, nk = l["u"][k], l["n"][k]
        for s, nn, lat in corners:
            p = l["E"][k] + s * uk + nn * nk; w = np.array([p[0], lat + Y, p[1]])
            emin = np.minimum(emin, w); emax = np.maximum(emax, w)
    l["env"] = (emin, emax); env_min_all = np.minimum(env_min_all, emin); env_max_all = np.maximum(env_max_all, emax)
n_s = int(round((s_hi - s_lo) / a.comp)); n_w = int(round(2 * w_hi / a.comp)); n_n = int(round((n_hi - n_lo_box) / a.comp))

# ------------------------------------------------------------------ 船体 & 背景域
def read_stl_bbox(fn):
    b = open(fn, "rb").read()
    if b[:5] == b"solid" and b"facet" in b[:400]:
        V = np.array([[float(x) for x in ln.split()[1:4]] for ln in b.decode(errors="ignore").splitlines() if ln.strip().startswith("vertex")])
    else:
        n = struct.unpack("<I", b[80:84])[0]
        V = np.frombuffer(b[84:84 + n * 50], dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]))["v"].reshape(-1, 3).astype(float)
    return V.min(0), V.max(0)
if a.no_hull:
    hmin, hmax = np.array([-0.126, -0.058, -0.035]), np.array([0.136, 0.058, 0.015])   # 只用来定域
else:
    hmin, hmax = read_stl_bbox(HULL)
Lh = hmax[0] - hmin[0]
gmin = np.minimum(hmin, env_min_all); gmax = np.maximum(hmax, env_max_all)
dom_min = np.array([gmin[0] - 2.5 * Lh, -(gmax[1] + 1.2 * Lh), gmin[2] - 1.2 * Lh])
dom_max = np.array([gmax[0] + 1.5 * Lh, gmax[1] + 1.2 * Lh, (WL + 0.5 * Lh) if FS else (gmax[2] + 0.12)])
base = a.base or a.bg * 4
nx, ny, nz = [max(8, int(round((dom_max[i] - dom_min[i]) / base))) for i in range(3)]
# 加密盒：每条腿的扫掠包络（level 2 = bg，外一圈 level 1）、自由面板、船体周围
pad_in = 0.03
refboxes = {}
for L, l in legs.items():
    emin, emax = l["env"]
    refboxes[f"ref2_{L}"] = (emin - pad_in, emax + pad_in, 2)
    refboxes[f"ref1_{L}"] = (emin - 3 * pad_in, emax + 3 * pad_in, 1)
if FS:
    refboxes["fs2"] = (np.array([gmin[0] - 0.8 * Lh, dom_min[1] * 0.6, WL - 0.03]), np.array([gmax[0] + 0.3 * Lh, dom_max[1] * 0.6, WL + 0.03]), 2)
    refboxes["fs1"] = (np.array([gmin[0] - 1.5 * Lh, dom_min[1] * 0.9, WL - 0.08]), np.array([gmax[0] + 0.5 * Lh, dom_max[1] * 0.9, WL + 0.06]), 1)
if not a.no_hull:
    refboxes["hull1"] = (hmin - np.array([0.10, 0.06, 0.06]), hmax + np.array([0.06, 0.06, 0.03]), 1)
# 洞切搜索盒：所有腿包络 + 3 cm
sb_min = env_min_all - 0.03; sb_max = env_max_all + 0.03
vox = a.vox if a.vox else min(ld.THK / 4.5, a.comp / 2)
sbd = [int(min(420, max(40, np.ceil((sb_max[i] - sb_min[i]) / vox)))) for i in range(3)]
vox_eff = max((sb_max[i] - sb_min[i]) / sbd[i] for i in range(3))
if vox_eff > 0.0016: print(f"!! 洞切体素被上限压到 {vox_eff*1000:.2f} mm（>1.5 mm 可能漏成全洞），考虑减少腿数或 --vox")

# 时间步
dt0 = 0.5 * a.comp / 2 / max(tip_speed + a.U, 0.05) * a.maxCo; dt0 = float(f"{dt0:.1e}")
endTime = a.cycles * T; writeInterval = T / a.writes

# ------------------------------------------------------------------ 写文件
def W(path, txt):
    os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w", encoding="utf-8").write(txt)
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
if os.path.exists(CASE): print(f"!! {CASE} 已存在，先删除"); shutil.rmtree(CASE)
os.makedirs(CASE)

MESHQ = "maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-13; minTetQuality 1e-15; minArea -1; minTwist 0.02; minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;"
LAYERS = "relativeSizes true; layers { } expansionRatio 1.2; finalLayerThickness 0.4; minThickness 0.1; nGrow 0; featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;"
def mesh_aux(d, app="snappyHexMesh"):
    W(os.path.join(d, "system", "controlDict"), hdr("dictionary", "controlDict") + f"\napplication {app};\nstartFrom startTime;\nstartTime 0;\nstopAt endTime;\nendTime 1;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\nwriteFormat binary;\nwritePrecision 8;\n")
    W(os.path.join(d, "system", "fvSchemes"), hdr("dictionary", "fvSchemes") + "\nddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\ndivSchemes { default none; }\nlaplacianSchemes { default Gauss linear corrected; }\ninterpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    W(os.path.join(d, "system", "fvSolution"), hdr("dictionary", "fvSolution") + "\n")
    W(os.path.join(d, "system", "meshQualityDict"), hdr("dictionary", "meshQualityDict") + "\n")

# ---- 部件网格子算例 paddleMesh_<leg>
for L, l in legs.items():
    PM = os.path.join(CASE, f"paddleMesh_{L}")
    os.makedirs(os.path.join(PM, "constant", "triSurface"))
    write_stl(l["tri"], os.path.join(PM, "constant", "triSurface", f"paddle_{L}.stl"))
    W(os.path.join(PM, "system", "blockMeshDict"), hdr("dictionary", "blockMeshDict") + f"""
scale 1;
vertices
(
{chr(10).join('    ' + fmt(v) for v in l["V"])}
);
blocks ( hex (0 1 2 3 4 5 6 7) ({n_s} {n_w} {n_n}) simpleGrading (1 1 1) );
edges ();
boundary
(
    overset_{L} {{ type overset; faces ( (0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7) ); }}
);
mergePatchPairs ();
""")
    W(os.path.join(PM, "system", "surfaceFeatureExtractDict"), hdr("dictionary", "surfaceFeatureExtractDict") + f"""
paddle_{L}.stl
{{
    extractionMethod    extractFromSurface;
    includedAngle       150;
    subsetFeatures {{ nonManifoldEdges no; openEdges yes; }}
    writeObj            no;
}}
""")
    W(os.path.join(PM, "system", "snappyHexMeshDict"), hdr("dictionary", "snappyHexMeshDict") + f"""
castellatedMesh true;
snap            true;
addLayers       false;
geometry {{ paddle_{L}.stl {{ type triSurfaceMesh; name paddle_{L}; }} {l["strip_geom"]} }}
castellatedMeshControls
{{
    maxLocalCells 2000000; maxGlobalCells 4000000; minRefinementCells 10; maxLoadUnbalance 0.1;
    nCellsBetweenLevels 3;
    features ( {{ file "paddle_{L}.eMesh"; level 1; }} );
    refinementSurfaces {{ paddle_{L} {{ level ({a.paddle_level} {a.paddle_level}); patchInfo {{ type wall; inGroups (wall); }} }} }}
    resolveFeatureAngle 30;
    refinementRegions {{ {l["strip"]} }}
    locationInMesh {fmt(l["lim"])};
    allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 50; nRelaxIter 5; nFeatureSnapIter 10; implicitFeatureSnap false; explicitFeatureSnap true; multiRegionFeatureSnap false; }}
addLayersControls {{ {LAYERS} }}
meshQualityControls {{ {MESHQ} }}
writeFlags ( scalarLevels );
mergeTolerance 1e-6;
""")
    mesh_aux(PM)

# ---- 背景网格
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
    inlet      {{ type patch; faces ((1 2 6 5)); }}    // +x 端（艏前方）：来流从这里吹向 -x
    outlet     {{ type patch; faces ((0 4 7 3)); }}    // -x 端（艉后方）
    sides      {{ type patch; faces ((0 1 5 4) (3 7 6 2)); }}
    bottom     {{ type patch; faces ((0 3 2 1)); }}
    {"atmosphere" if FS else "top"} {{ type patch; faces ((4 5 6 7)); }}
);
mergePatchPairs ();
""")
if not a.no_hull:
    os.makedirs(os.path.join(CASE, "constant", "triSurface"), exist_ok=True)
    shutil.copyfile(HULL, os.path.join(CASE, "constant", "triSurface", "hull.stl"))
    W(os.path.join(CASE, "system", "surfaceFeatureExtractDict"), hdr("dictionary", "surfaceFeatureExtractDict") + """
hull.stl
{
    extractionMethod    extractFromSurface;
    includedAngle       150;
    subsetFeatures { nonManifoldEdges no; openEdges yes; }
    writeObj            no;
}
""")
geom_bg = "".join(f"    {k} {{ type searchableBox; min {fmt(v[0])}; max {fmt(v[1])}; }}\n" for k, v in refboxes.items())
regs_bg = "".join(f"        {k} {{ mode inside; levels ((1e15 {v[2]})); }}\n" for k, v in refboxes.items())
loc_bg = (dom_min[0] + 0.05, dom_min[1] + 0.05, dom_min[2] + 0.05)
W(os.path.join(CASE, "system", "snappyHexMeshDict"), hdr("dictionary", "snappyHexMeshDict") + f"""
// 背景网格：船体贴体（level 2 = {a.bg*1000:.0f} mm）+ 区域加密（腿扫掠区 / 自由面 / 船体周围）
castellatedMesh true;
snap            {"true" if not a.no_hull else "false"};
addLayers       false;
geometry
{{
{"    hull.stl { type triSurfaceMesh; name hull; }" if not a.no_hull else ""}
{geom_bg}}}
castellatedMeshControls
{{
    maxLocalCells 4000000; maxGlobalCells 8000000; minRefinementCells 10; maxLoadUnbalance 0.1;
    nCellsBetweenLevels 3;
    features ( {f'{{ file "hull.eMesh"; level {a.hull_level}; }}' if not a.no_hull else ""} );
    refinementSurfaces {{ {f"hull {{ level ({a.hull_level} {a.hull_level}); patchInfo {{ type wall; inGroups (wall); }} }}" if not a.no_hull else ""} }}
    resolveFeatureAngle 30;
    refinementRegions
    {{
{regs_bg}    }}
    locationInMesh {fmt(loc_bg)};
    allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; nFeatureSnapIter 10; implicitFeatureSnap false; explicitFeatureSnap true; multiRegionFeatureSnap false; }}
addLayersControls {{ {LAYERS} }}
meshQualityControls {{ {MESHQ} }}
writeFlags ();
mergeTolerance 1e-6;
""")
# cellZones：合并前在每个子网格里把全部单元做成 zone（背景 c0，腿 c1..cN），mergeMeshes 会把 zone 带过来，
# 不依赖 regionToCell（合并后背景与部件在空间上重叠，findCell 会选错区域）
def zone_dict(name):
    return hdr("dictionary", "topoSetDict") + f"""
actions
(
    {{ name {name}; type cellSet; action new; source boxToCell; box (-100 -100 -100) (100 100 100); }}
    {{ name {name}; type cellZoneSet; action new; source setToCellZone; set {name}; }}
);
"""
W(os.path.join(CASE, "system", "topoSetDict"), zone_dict("c0"))
for i, L in enumerate(legs, 1):
    W(os.path.join(CASE, f"paddleMesh_{L}", "system", "topoSetDict"), zone_dict(f"c{i}"))
sf = "defaultFieldValues ( volScalarFieldValue zoneID 123" + (" volScalarFieldValue alpha.water 0" if FS else "") + " );\nregions\n(\n"
for i in range(len(legs) + 1):
    sf += f"    zoneToCell {{ zone c{i}; fieldValues ( volScalarFieldValue zoneID {i} ); }}\n"
if FS:
    sf += f"    boxToCell {{ box (-100 -100 -100) (100 100 {WL}); fieldValues ( volScalarFieldValue alpha.water 1 ); }}\n"
    sf += f"    boxToFace {{ box (-100 -100 -100) (100 100 {WL}); fieldValues ( volScalarFieldValue alpha.water 1 ); }}\n"
W(os.path.join(CASE, "system", "setFieldsDict"), hdr("dictionary", "setFieldsDict") + "\n" + sf + ");\n")

# controlDict + 受力
hull_c = (hmin + hmax) / 2 if not a.no_hull else np.array([0.0, 0.0, WL])
CofR = f"({hull_c[0]:.4f} 0 {WL:.4f})"
rho_fo = "rhoInf          998.8;" if FS else f"rho             rhoInf;\n        rhoInf          {RHO};"
def fo(name, patches):
    return f"""    forces_{name}
    {{
        type            forces;
        libs            ("libforces.so");
        patches         ({patches});
        {rho_fo}
        CofR            {CofR};      // 力矩参考点：船体中点、水线高度
        writeControl    timeStep;
        writeInterval   1;
        log             false;
    }}
"""
fos = "".join(fo(L, f"paddle_{L}") for L in legs)
if not a.no_hull: fos += fo("hull", "hull")
fos += fo("all", " ".join([f"paddle_{L}" for L in legs] + ([] if a.no_hull else ["hull"])))
W(os.path.join(CASE, "system", "controlDict"), hdr("dictionary", "controlDict") + f"""
application     {"overInterDyMFoam" if FS else "overPimpleDyMFoam"};
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
{"maxAlphaCo      " + str(a.maxAlphaCo) + ";" if FS else ""}
maxDeltaT       {4*dt0};

functions
{{
{fos}    vorticity1
    {{
        type            vorticity;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
    }}
    zoneOut                   // 每个写出时刻把 zoneID 也写出来，后处理里用它把背景网格和桨部件网格分开（画自由面）
    {{
        type            writeObjects;
        libs            ("libutilityFunctionObjects.so");
        objects         (zoneID);
        writeControl    writeTime;
    }}
    maxU                      // 每步报告最大速度及其位置：Δt 被卡住时一眼看出是哪里的伪速度
    {{
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        mode            magnitude;
        fields          (U);
        location        yes;
        writeControl    timeStep;
        writeInterval   1;
        log             true;
    }}
}}
""")
if a.limitU > 0:
    W(os.path.join(CASE, "system", "fvOptions"), hdr("dictionary", "fvOptions") + f"""
// 速度上限：桨尖 ~0.4 m/s，水里真实速度 < 1 m/s；空气区和桨面薄片单元里的伪速度可到几 m/s，会把 Δt 卡在 0.1 ms
limitU
{{
    type            limitVelocity;
    active          yes;
    selectionMode   all;
    max             {a.limitU};
}}
""")
turbDiv = "" if a.turb == "laminar" else "\n    div(phi,k)      Gauss limitedLinear 1;\n    div(phi,omega)  Gauss limitedLinear 1;"
if FS:
    divs = f"""    div(rhoPhi,U)   Gauss limitedLinearV 1;
    div(U)          Gauss linear;
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;{turbDiv}"""
    extra = "oversetInterpolationSuppressed { grad(p_rgh); surfaceIntegrate(phiHbyA); }\nfluxRequired { default no; p_rgh; pcorr; alpha.water; }"
else:
    divs = f"""    div(phi,U)      Gauss limitedLinearV 1;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;{turbDiv}"""
    extra = "fluxRequired { default no; p; }"
W(os.path.join(CASE, "system", "fvSchemes"), hdr("dictionary", "fvSchemes") + f"""
ddtSchemes {{ default Euler; }}
gradSchemes {{ default Gauss linear; }}
divSchemes
{{
    default         none;
{divs}
}}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default corrected; }}
// 洞切：搜索盒 = 四腿运动包络 + 3 cm，体素 ≈{vox_eff*1000:.2f} mm（必须 ≲ 桨厚/4）；v2606 若报 searchBoxDivisions 类型错误，改成每个 zone 一组
oversetInterpolation
{{
    method              {a.overset};
    searchBox           {fmt(sb_min)} {fmt(sb_max)};
    searchBoxDivisions  ({sbd[0]} {sbd[1]} {sbd[2]});
}}
{extra}
wallDist {{ method meshWave; }}
""")
if FS:
    W(os.path.join(CASE, "system", "fvSolution"), hdr("dictionary", "fvSolution") + """
solvers
{
    "cellDisplacement.*" { solver PCG; preconditioner DIC; tolerance 1e-06; relTol 0; maxIter 300; }
    "alpha.water.*"
    {
        nAlphaCorr      3;
        nAlphaSubCycles 2;
        cAlpha          1;
        icAlpha         0;
        MULESCorr       yes;
        nLimiterIter    5;
        alphaApplyPrevCorr no;
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }
    "pcorr.*"   { solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0; }
    p_rgh       { solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.05; }
    p_rghFinal  { $p_rgh; relTol 0; }
    "(U|k|omega).*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; }
}
PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 0;
    ddtCorr             yes;
    correctPhi          no;
    massFluxInterpolation no;
    moveMeshOuterCorrectors no;
    turbOnFinalIterOnly no;
    oversetAdjustPhi    no;
    checkMeshCourantNo  yes;
}
relaxationFactors { equations { ".*" 1; } }
cache { grad(U); }
""")
else:
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
    momentumPredictor false; nOuterCorrectors 1; nCorrectors 3; nNonOrthogonalCorrectors 0;
    ddtCorr true; correctPhi false; oversetAdjustPhi false; checkMeshCourantNo yes;
}
relaxationFactors { equations { ".*" 1; } }
""")
def hier_n(n):
    best = (n, 1, 1)
    for nx_ in range(n, 0, -1):
        for ny_ in range(n, 0, -1):
            if n % (nx_ * ny_) == 0:
                nz_ = n // (nx_ * ny_)
                if (max(nx_, ny_, nz_), -nx_) < (max(best), -best[0]): best = (nx_, ny_, nz_)
    return best
hn = hier_n(a.np)
W(os.path.join(CASE, "system", "decomposeParDict"), hdr("dictionary", "decomposeParDict") + f"""
// overset 并行用 hierarchical（scotch 在 overset 上触发过 "Attempt to cast type patch to type lduInterface"）
numberOfSubdomains {a.np};
method          hierarchical;
hierarchicalCoeffs {{ n ({hn[0]} {hn[1]} {hn[2]}); order xyz; }}
""")

# ---- constant
os.makedirs(os.path.join(CASE, "constant"), exist_ok=True)
zones = ""
for i, (L, l) in enumerate(legs.items(), 1):
    zones += f"""    c{i}
    {{
        solidBodyMotionFunction tabulated6DoFMotion;
        tabulated6DoFMotionCoeffs
        {{
            CofG            ({l['E0'][0]:.6f} {l['y']:.6f} {l['E0'][1]:.6f});
            timeDataFileName "<constant>/6DoF_{L}.dat";
        }}
    }}
"""
    ry = -np.degrees(l["psi"] - l["psi0"]); trans = l["E"] - l["E0"]
    with open(os.path.join(CASE, "constant", f"6DoF_{L}.dat"), "w") as f:
        f.write(f"{nT}\n(\n")
        for k in range(nT): f.write(f"({t[k]:.6f} (({trans[k,0]:.6f} 0 {trans[k,1]:.6f}) (0 {ry[k]:.5f} 0)))\n")
        f.write(")\n")
W(os.path.join(CASE, "constant", "dynamicMeshDict"), hdr("dictionary", "dynamicMeshDict", "constant") + f"""
dynamicFvMesh       dynamicOversetFvMesh;
dynamicOversetFvMeshCoeffs {{ }}

// v2606 若报 "solver" 关键字：改成 motionSolver multiSolidBodyMotionSolver; 并把 Coeffs 内容保留
solver          multiSolidBodyMotionSolver;
multiSolidBodyMotionSolverCoeffs
{{
{zones}}}
""")
if FS:
    W(os.path.join(CASE, "constant", "transportProperties"), hdr("dictionary", "transportProperties", "constant") + f"""
phases (water air);
water {{ transportModel Newtonian; nu {NU}; rho {RHO}; }}
air   {{ transportModel Newtonian; nu 1.48e-05; rho 1; }}
sigma           0.07;
""")
    W(os.path.join(CASE, "constant", "g"), hdr("uniformDimensionedVectorField", "g", "constant") + "\ndimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n")
    W(os.path.join(CASE, "constant", "hRef"), hdr("uniformDimensionedScalarField", "hRef", "constant") + f"\ndimensions [0 1 0 0 0 0 0];\nvalue {WL};\n")
else:
    W(os.path.join(CASE, "constant", "transportProperties"), hdr("dictionary", "transportProperties", "constant") + f"\ntransportModel Newtonian;\nnu {NU};\n")
W(os.path.join(CASE, "constant", "turbulenceProperties"), hdr("dictionary", "turbulenceProperties", "constant") +
  ("\nsimulationType laminar;\n" if a.turb == "laminar" else "\nsimulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }\n"))

# ---- 0.orig
def OVS(bc):                       # 每个 overset patch 一条精确名条目
    return "".join(f"    overset_{L} {{ {bc} }}\n" for L in legs)
Uin = f"(-{a.U:.4f} 0 0)" if a.U > 0 else "(0 0 0)"
WALLS = '"(paddle_.*|hull)"'
TOPN = "atmosphere" if FS else "top"
if FS:
    W(os.path.join(CASE, "0.orig", "U"), hdr("volVectorField", "U", "0") + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform {Uin};
boundaryField
{{
{OVS("type overset; patchType overset; value uniform (0 0 0);")}    {WALLS} {{ type movingWallVelocity; value uniform (0 0 0); }}
    inlet        {{ type fixedValue; value uniform {Uin}; }}
    outlet       {{ type outletPhaseMeanVelocity; alpha alpha.water; Umean {a.U:.4f}; value uniform {Uin}; }}
    atmosphere   {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    "(sides|bottom)" {{ type slip; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "p_rgh"), hdr("volScalarField", "p_rgh", "0") + f"""
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
{OVS("type overset;")}    {WALLS} {{ type fixedFluxPressure; value uniform 0; }}
    inlet        {{ type fixedFluxPressure; value uniform 0; }}
    outlet       {{ type zeroGradient; }}
    atmosphere   {{ type totalPressure; p0 uniform 0; U U; phi phi; rho rho; psi none; gamma 1; value uniform 0; }}
    "(sides|bottom)" {{ type fixedFluxPressure; value uniform 0; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "alpha.water"), hdr("volScalarField", "alpha.water", "0") + f"""
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
{OVS("type overset; value uniform 1;")}    {WALLS} {{ type zeroGradient; }}
    inlet        {{ type fixedValue; value uniform 0; }}      // setFields 的 boxToFace 会把水线以下改成 1
    outlet       {{ type variableHeightFlowRate; lowerBound 0; upperBound 1; value uniform 0; }}
    atmosphere   {{ type inletOutlet; inletValue uniform 0; value uniform 0; }}
    "(sides|bottom)" {{ type zeroGradient; }}
}}
""")
else:
    W(os.path.join(CASE, "0.orig", "U"), hdr("volVectorField", "U", "0") + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform {Uin};
boundaryField
{{
{OVS("type overset;")}    {WALLS} {{ type movingWallVelocity; value uniform (0 0 0); }}
    inlet        {{ type fixedValue; value uniform {Uin}; }}
    outlet       {{ type inletOutlet; inletValue uniform (0 0 0); value uniform {Uin}; }}
    "(sides|bottom|top)" {{ type slip; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "p"), hdr("volScalarField", "p", "0") + f"""
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
{OVS("type overset;")}    {WALLS} {{ type zeroGradient; }}
    inlet        {{ type zeroGradient; }}
    outlet       {{ type fixedValue; value uniform 0; }}
    "(sides|bottom|top)" {{ type zeroGradient; }}
}}
""")
W(os.path.join(CASE, "0.orig", "zoneID"), hdr("volScalarField", "zoneID", "0") + f"""
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    ".*"         {{ type zeroGradient; }}      // 通配写在前面：OpenFOAM 正则按倒序匹配，后面的精确名优先
{OVS("type overset; value uniform 0;")}}}
""")
if a.turb != "laminar":
    k_in = 1.5 * (0.02 * max(a.U, 0.1))**2; om_in = np.sqrt(k_in) / (0.09**0.25 * 0.07 * 0.03)
    for name, val, wallbc in (("k", k_in, "kqRWallFunction"), ("omega", om_in, "omegaWallFunction")):
        W(os.path.join(CASE, "0.orig", name), hdr("volScalarField", name, "0") + f"""
dimensions {"[0 2 -2 0 0 0 0]" if name == "k" else "[0 0 -1 0 0 0 0]"};
internalField uniform {val:.3e};
boundaryField
{{
{OVS("type overset;")}    {WALLS} {{ type {wallbc}; value uniform {val:.3e}; }}
    inlet {{ type fixedValue; value uniform {val:.3e}; }}
    outlet {{ type inletOutlet; inletValue uniform {val:.3e}; value uniform {val:.3e}; }}
    {TOPN} {{ type inletOutlet; inletValue uniform {val:.3e}; value uniform {val:.3e}; }}
    "(sides|bottom)" {{ type slip; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "nut"), hdr("volScalarField", "nut", "0") + f"""
dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    ".*" {{ type calculated; value uniform 0; }}
{OVS("type overset;")}    {WALLS} {{ type nutkWallFunction; value uniform 0; }}
}}
""")

# ---- Allrun
merge = "\n".join(f"runApplication -s {L} mergeMeshes . paddleMesh_{L} -overwrite" for L in legs)
pm = "\n".join(f"""( cd paddleMesh_{L} || exit 1
  runApplication blockMesh
  runApplication surfaceFeatureExtract
  runApplication snappyHexMesh -overwrite
  runApplication topoSet
  runApplication checkMesh )""" for L in legs)
W(os.path.join(CASE, "Allrun"), f"""#!/bin/sh
cd "${{0%/*}}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 1 四块桨的部件网格
{pm}

# 2 背景网格（船体贴体 + 区域加密）
{"runApplication surfaceFeatureExtract" if not a.no_hull else ""}
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication topoSet

# 3 合并（cellZones c0..cN 随子网格带入）→ zoneID{" / alpha.water" if FS else ""}
{merge}
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
for d in paddleMesh_*; do rm -rf $d/constant/polyMesh $d/constant/extendedFeatureEdgeMesh $d/constant/triSurface/*.eMesh $d/log.*; done
rm -rf constant/extendedFeatureEdgeMesh constant/triSurface/*.eMesh postProcessing 0
""")
os.chmod(os.path.join(CASE, "Allclean"), 0o755)

# ---- 记录
json.dump(dict(fan=a.fan, U=a.U, cycles=a.cycles, T=T, maxCo=a.maxCo, maxAlphaCo=a.maxAlphaCo, limitU=a.limitU, overset=a.overset, DUTY=ld.DUTY, SWEEP=ld.SWEEP, PHI_EXT=ld.PHI_EXT, PHI_RET=ld.PHI_RET,
               FAN_DEG=ld.FAN_DEG, FAN_CLOSE=ld.FAN_CLOSE, CLOSED_W=ld.CLOSED_W, turb=a.turb, bg=a.bg, comp=a.comp, base=base,
               WL=WL, hull=(None if a.no_hull else os.path.basename(HULL)), free_surface=FS, gait=a.gait,
               legs={L: dict(phase=l["ph"], dx=l["dx"], y=l["y"], E0=l["E0"].tolist(), psi0_deg=float(np.degrees(l["psi0"]))) for L, l in legs.items()},
               dt0=dt0, tip_speed_max=float(tip_speed), psi0_deg=float(np.degrees(legs[LEGS[0]]["psi0"])), CofR=[float(hull_c[0]), 0.0, WL],
               min_gap=MIN_GAP, fan_clip=clip_info, gaps={L: [float(l["gap"][0]), float(l["gap"][1])] for L, l in legs.items()},
               paddle=a.paddle, head_nrm=list(HEAD_NRM), stem_nrm=list(STEM_NRM), thk_open=THK_OPEN, paddle_level=a.paddle_level,
               outlines={L: [[float(p[0]), float(p[1])] for p in l["outline"]] for L, l in legs.items()}),
          open(os.path.join(CASE, "gait.json"), "w"), indent=1, ensure_ascii=False)

print(f"算例：{CASE}   脚蹼 {a.fan}   U = {a.U} m/s   {a.cycles} 周期 → endTime {endTime:.3f} s   {'VOF 自由面' if FS else '单相'}   船体 {'无' if a.no_hull else os.path.basename(HULL)}")
print("腿：" + ", ".join(f"{L}(相位 {l['ph']}, dx {l['dx']:+.4f}, y {l['y']:+.4f})" for L, l in legs.items()) + f"   步态 {a.gait}")
print(f"部件盒 {n_s}×{n_w}×{n_n} = {n_s*n_w*n_n} 单元 × {len(legs)}（{a.comp*1000:.2f} mm，桨面 1 级加密）")
print(f"背景域 x[{dom_min[0]:.3f},{dom_max[0]:.3f}] y[{dom_min[1]:.3f},{dom_max[1]:.3f}] z[{dom_min[2]:.3f},{dom_max[2]:.3f}]  blockMesh {nx}×{ny}×{nz}（{base*1000:.0f} mm），加密到 {a.bg*1000:.0f} mm")
print(f"洞切搜索盒 {fmt(sb_min)}–{fmt(sb_max)} 体素 {sbd[0]}×{sbd[1]}×{sbd[2]}（≈{vox_eff*1000:.2f} mm）")
print(f"桨尖最大速度 {tip_speed:.2f} m/s → Δt0 = {dt0} s，maxCo {a.maxCo}{f'，maxAlphaCo {a.maxAlphaCo}' if FS else ''}，limitVelocity {a.limitU} m/s，overset {a.overset}；每周期 {a.writes} 帧")
print("跑：cd", CASE, "&& ./Allrun")
