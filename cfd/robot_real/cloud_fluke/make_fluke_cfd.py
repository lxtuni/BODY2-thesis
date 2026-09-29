#!/usr/bin/env python3
"""
升力型尾鳍（及翼型基准）单体非定常 CFD 算例生成器 —— overset + overPimpleDyMFoam
==================================================================================
与 make_leg_cfd.py 同一套流程（部件网格 paddleMesh + 背景网格 → mergeMeshes → tabulated6DoFMotion），
区别只在几何与运动：
  几何  --stl FILE        任意封闭 STL，坐标原点 = 转轴（尾鳍销轴），弦线在 x–z 面内、角度 --ref-deg
        --naca 0012,C,B,F  矩形 NACA 翼：弦 C、展 B，转轴在前缘后 F·C（翼型基准 V1）
  运动  --motion CSV      fluke_traj.py 的输出：t_s, hx_m, hz_m, fin_deg（Flare 约定：hx 朝船尾为正、hz 向上、
                          fin_deg = 弦线绝对角，0 = 水平朝船尾，正 = 尾缘上翘）
补丁名仍叫 paddle、部件网格目录仍叫 paddleMesh，所以 progress.sh / resume_case.sh / compare_verify.py 都照用。

    python3 make_fluke_cfd.py --case ~/run/fluke_L0 --motion fluke/motion_L0.csv --stl fluke/A2_fin.stl --U 0.2229
    python3 make_fluke_cfd.py --case ~/run/fluke_V1 --motion fluke/motion_V1.csv --naca 0012,0.1,0.6,0.3333 --origin 0,0 --U 0.4 \
                              --comp 0.002 --bg 0.004 --vox 0.002 --Lref 0.3 --cycles 3

坐标：机体系 x 前进、y 左、z 上；--origin = 运动表原点（腿：髋轴 = O2 (0.06345, 0.004)）。来流从 +x 吹向 −x。
"""
import argparse, os, sys, json, struct, shutil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True)
ap.add_argument("--motion", required=True, help="fluke_traj.py 输出的一个周期运动表")
gg = ap.add_mutually_exclusive_group(required=True)
gg.add_argument("--stl", help="封闭 STL（m），原点 = 转轴")
gg.add_argument("--naca", help="四位 NACA 矩形翼：XXXX,弦 C[m],展 B[m],转轴位置 F（前缘后 F·C）")
ap.add_argument("--ref-deg", type=float, default=-25.0, help="STL 里弦线的绝对角（Flare 约定）；flare-live 的 LF_fin.stl = −25（零位姿态）")
ap.add_argument("--origin", default="0.06345,0.004", help="运动表原点在机体系的 x,z [m]（默认髋轴 O2）")
ap.add_argument("--area", type=float, default=None, help="参考面积 [m²]（写进 gait.json 供后处理；默认 STL = 11.5e-4（A2 平面形），NACA = C·B）")
ap.add_argument("--U", type=float, default=0.0)
ap.add_argument("--cycles", type=float, default=4.0)
ap.add_argument("--bg", type=float, default=0.0025, help="背景网格在运动包络内的单元尺寸 [m]")
ap.add_argument("--comp", type=float, default=0.0012, help="部件网格基准单元尺寸 [m]（表面再细化 1 级）")
ap.add_argument("--vox", type=float, default=None, help="洞切体素 [m]（默认 min(最薄处/4.5, comp/2)；尾鳍梢部只有 ~2.7 mm 厚）")
ap.add_argument("--tmin", type=float, default=0.0027, help="几何最薄处厚度 [m]，只用于默认体素")
ap.add_argument("--Lref", type=float, default=0.262, help="背景域尺度 [m]（域 = 包络 ± 1–1.4 Lref）")
ap.add_argument("--turb", default="laminar", choices=["laminar", "kOmegaSST"])
ap.add_argument("--np", type=int, default=8)
ap.add_argument("--nsub", type=int, default=800, help="运动表每周期采样点数")
ap.add_argument("--maxCo", type=float, default=0.8)
ap.add_argument("--ramp", type=float, default=0.5, help="起步平滑段长度 [周期]：t=0 位姿、速度 0 起步，smoothstep 过渡到运动表（0 = 关）；只影响第 1 周期，平均取最后 2 个周期不受影响")
ap.add_argument("--limitU", type=float, default=0.0, help="fvOptions limitVelocity 上限 [m/s]（0 = 不加）；只压重叠边界/薄片单元里的伪速度，真实流速 < 1 m/s")
ap.add_argument("--nsamp", type=int, default=40, help="截面 / 尾鳍表面采样：每周期帧数（最后 2 个周期）")
ap.add_argument("--iso", action="store_true", help="最后一个周期另存 Q 准则等值面。⚠ v2606 实测：sampledIsoSurface 在重叠网格上段错误（getIsoFields），不要用；三维涡结构用保留的三维流场在 ParaView 里算 Q")
ap.add_argument("--qstar", type=float, default=2.0, help="Q 等值面取值 Q·c²/U_rel²（c = 弦长，U_rel = 来流 + 表面最大速度）")
ap.add_argument("--functions-only", action="store_true", help="只把 controlDict 的 functions 块换成新的（给已在运行的算例加输出；runTimeModifiable 会自动读入），其他文件不动")
ap.add_argument("--writes", type=int, default=20, help="每周期写出帧数")
a = ap.parse_args()

CASE = os.path.abspath(os.path.expanduser(a.case))
RHO, NU = 1000.0, 1.0e-6
Y0 = 0.0
OX, OZ = [float(v) for v in a.origin.split(",")]

# ------------------------------------------------------------------ 几何（体坐标：原点 = 转轴）
def read_stl(fn):
    b = open(fn, "rb").read()
    if b[:5] == b"solid" and b"facet" in b[:400]:
        v = [list(map(float, l.split()[1:4])) for l in b.decode(errors="ignore").splitlines() if l.strip().startswith("vertex")]
        return np.array(v, float).reshape(-1, 3, 3)
    n = struct.unpack("<I", b[80:84])[0]
    rec = np.frombuffer(b[84:84 + 50 * n], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
    return rec["v"].astype(float)

def naca_stl(code, C, B, F, n=80):
    t = int(code[-2:]) / 100.0
    be = np.linspace(0, np.pi, n); xc = 0.5 * (1 - np.cos(be))
    yt = 5 * t * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc**2 + 0.2843 * xc**3 - 0.1036 * xc**4)
    # 截面在 (x, z)：前缘在 +x（迎着来流），x = F·C − xc·C；上表面 z = +yt·C
    up = np.c_[F * C - xc * C, yt * C]; lo = np.c_[F * C - xc[::-1] * C, -yt[::-1] * C]
    poly = np.vstack([up, lo[1:-1]])                   # 闭合环（逆时针与否下面统一定向）
    m = len(poly); ys = (-B / 2, B / 2); tris = []
    for i in range(m):
        j = (i + 1) % m
        p0, p1 = poly[i], poly[j]
        a0, a1 = (p0[0], ys[0], p0[1]), (p1[0], ys[0], p1[1]); b0, b1 = (p0[0], ys[1], p0[1]), (p1[0], ys[1], p1[1])
        tris += [[a0, a1, b1], [a0, b1, b0]]
    c = poly.mean(0)
    for y in ys:
        for i in range(m):
            j = (i + 1) % m
            tris.append([(c[0], y, c[1]), (poly[i][0], y, poly[i][1]), (poly[j][0], y, poly[j][1])])
    tri = np.array(tris, float)
    # 外法向定向：法向与“面心 − 体心”同向
    cen = np.array([c[0], 0.0, c[1]]); nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    flip = np.einsum("ij,ij->i", nrm, tri.mean(1) - cen) < 0; tri[flip] = tri[flip][:, ::-1]
    return tri

def Ry(beta_deg):
    b = np.radians(beta_deg); c, s = np.cos(b), np.sin(b)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])   # 绕 +y 右手转：把 Flare 弦角 c 增加 β（见下）

if a.stl:
    tri_b = read_stl(os.path.expanduser(a.stl)); REF = a.ref_deg; AREA = a.area or 11.5e-4; GEOM = dict(stl=os.path.abspath(os.path.expanduser(a.stl)))
else:
    code, C, B, F = a.naca.split(","); C, B, F = float(C), float(B), float(F)
    tri_b = naca_stl(code, C, B, F); REF = 0.0; AREA = a.area or C * B; GEOM = dict(naca=code, chord=C, span=B, pivot_frac=F)
# 体坐标 → “弦线水平朝船尾”坐标（chord 系）：绕 y 转 −REF
tri_c = np.einsum("ij,nkj->nki", Ry(-REF), tri_b)
# 自检：弦向 = 世界 (−cos c, 0, sin c)；Ry(β) 把 c 变成 c+β
v = np.array([-np.cos(np.radians(10)), 0, np.sin(np.radians(10))]); v2 = Ry(15) @ v
assert np.allclose(v2, [-np.cos(np.radians(25)), 0, np.sin(np.radians(25))]), "Ry 符号约定错"

# ------------------------------------------------------------------ 运动表
D = np.genfromtxt(a.motion, delimiter=",", names=True)
tc = D["t_s"]; T = float(tc[-1] + (tc[1] - tc[0]))
nT = int(a.nsub * a.cycles * 1.02) + 1                        # 运动表比 endTime 多 2 %：自适应步长最后一步可能越过 endTime
t = np.linspace(0, a.cycles * T * 1.02, nT); tm = t % T
ip = lambda y: np.interp(tm, tc, y, period=T)
fin = np.degrees(np.unwrap(np.radians(D["fin_deg"])))
H = np.stack([OX - ip(D["hx_m"]), np.full(nT, Y0), OZ + ip(D["hz_m"])], 1)     # 铰链世界坐标
FIN = ip(fin)                                                                   # 弦线绝对角（Flare 约定）
if a.ramp > 0:                                                                  # 起步平滑：P = P0 + s(t)(P − P0)，s = smoothstep → t=0 位姿不变、速度 0
    x_ = np.clip(t / (a.ramp * T), 0, 1); s_ = x_ * x_ * (3 - 2 * x_)
    H = H[0] + s_[:, None] * (H - H[0]); FIN = FIN[0] + s_ * (FIN - FIN[0])
H0, FIN0 = H[0], FIN[0]
Rk = [Ry(f_) for f_ in FIN]                                                     # chord 系 → 世界（绕铰链）

# t=0 位姿下的 STL（世界）
tri = H0 + np.einsum("ij,nkj->nki", Rk[0], tri_c)

def write_stl(tri, fn):
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
    rec = np.zeros(len(tri), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])); rec["n"] = nrm; rec["v"] = tri
    open(fn, "wb").write(b"\0" * 80 + struct.pack("<I", len(tri)) + rec.tobytes())

# ------------------------------------------------------------------ 部件网格盒（chord 系包围盒 + 余量，随尾鳍转）
pc = tri_c.reshape(-1, 3); cmin, cmax = pc.min(0), pc.max(0)
mx, my, mz = 0.016, 0.012, 0.016 if a.stl else 0.02
lo = cmin - [mx, my, mz]; hi = cmax + [mx, my, mz]
corners_c = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
def to_world(k, p): return H[k] + Rk[k] @ p
V = [to_world(0, np.array(p)) for p in [(lo[0], lo[1], lo[2]), (hi[0], lo[1], lo[2]), (hi[0], hi[1], lo[2]), (lo[0], hi[1], lo[2]),
                                       (lo[0], lo[1], hi[2]), (hi[0], lo[1], hi[2]), (hi[0], hi[1], hi[2]), (lo[0], hi[1], hi[2])]]
v0, v1, v3, v4 = V[0], V[1], V[3], V[4]
if np.dot(np.cross(v1 - v0, v3 - v0), v4 - v0) < 0: V = V[4:] + V[:4]
n_x = int(round((hi[0] - lo[0]) / a.comp)); n_y = int(round((hi[1] - lo[1]) / a.comp)); n_z = int(round((hi[2] - lo[2]) / a.comp))
lim = to_world(0, np.array([hi[0] - 0.004, hi[1] - 0.004, hi[2] - 0.004]))       # 部件盒内、物体外

# 运动包络、最大表面速度
env_min = np.full(3, np.inf); env_max = np.full(3, -np.inf)
samp = pc[np.random.default_rng(0).choice(len(pc), min(400, len(pc)), replace=False)]
Pk = []
for k in range(nT):
    wk = H[k] + samp @ Rk[k].T; Pk.append(wk)
    ck = H[k] + corners_c @ Rk[k].T
    env_min = np.minimum(env_min, ck.min(0)); env_max = np.maximum(env_max, ck.max(0))
Pk = np.array(Pk); vmax = float(np.linalg.norm(np.gradient(Pk, t, axis=0), axis=2).max())

pad_in = 0.03
ref2 = (env_min - pad_in, env_max + pad_in); ref1 = (env_min - 3 * pad_in, env_max + 3 * pad_in)
L = a.Lref
dom_min = np.array([env_min[0] - 1.4 * L, min(env_min[1] - 0.9 * L, -0.25), env_min[2] - 1.0 * L])
dom_max = np.array([env_max[0] + 1.0 * L, max(env_max[1] + 0.9 * L, 0.25), max(env_max[2] + 0.12, 0.06)])
sb_min = env_min - 0.03; sb_max = env_max + 0.03
vox = a.vox if a.vox else min(a.tmin / 4.5, a.comp / 2)
sbd = [int(min(240, max(40, np.ceil((sb_max[i] - sb_min[i]) / vox)))) for i in range(3)]
vox_eff = max((sb_max[i] - sb_min[i]) / sbd[i] for i in range(3))
base = a.bg * 4
nx, ny, nz = [max(8, int(round((dom_max[i] - dom_min[i]) / base))) for i in range(3)]

trans = H - H0; ry = FIN - FIN0
dt0 = float(f"{0.5 * a.comp / 2 / max(vmax, 0.05) * a.maxCo:.1e}")
endTime = a.cycles * T; writeInterval = T / a.writes

# ------------------------------------------------------------------ 写文件（与 make_leg_cfd.py 相同的字典，只换几何/运动）
def W(path, txt):
    if a.functions_only: return                                    # 给运行中的算例加输出：不碰其他文件
    os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w", encoding="utf-8").write(txt)
def hdr(cls, obj, loc="system"):
    return f"FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    location    \"{loc}\";\n    object      {obj};\n}}\n"
if a.functions_only:
    if not os.path.exists(os.path.join(CASE, "system", "controlDict")): sys.exit(f"!! {CASE} 不是算例目录")
else:
    if os.path.exists(CASE): print(f"!! {CASE} 已存在，先删除"); shutil.rmtree(CASE)
    os.makedirs(CASE)
fmt = lambda v: "(" + " ".join(f"{x:.6g}" for x in v) + ")"

PM = os.path.join(CASE, "paddleMesh")
if not a.functions_only:
    os.makedirs(os.path.join(PM, "constant", "triSurface"))
    write_stl(tri, os.path.join(PM, "constant", "triSurface", "paddle.stl"))
W(os.path.join(PM, "system", "blockMeshDict"), hdr("dictionary", "blockMeshDict") + f"""
scale 1;
vertices
(
{chr(10).join('    ' + fmt(v) for v in V)}
);
blocks ( hex (0 1 2 3 4 5 6 7) ({n_x} {n_y} {n_z}) simpleGrading (1 1 1) );
edges ();
boundary
(
    oversetPatch {{ type overset; faces ( (0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7) ); }}
);
mergePatchPairs ();
""")
W(os.path.join(PM, "system", "surfaceFeatureExtractDict"), hdr("dictionary", "surfaceFeatureExtractDict") + """
paddle.stl
{
    extractionMethod    extractFromSurface;
    includedAngle       150;
    subsetFeatures { nonManifoldEdges no; openEdges yes; }
    writeObj            no;
}
""")
MQC = ("maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80; minVol 1e-14; minTetQuality 1e-15; "
       "minArea -1; minTwist 0.02; minDeterminant 0.001; minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;")
ALC = ("relativeSizes true; layers { } expansionRatio 1.2; finalLayerThickness 0.4; minThickness 0.1; nGrow 0; featureAngle 60; nRelaxIter 3; "
       "nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; "
       "minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;")
W(os.path.join(PM, "system", "snappyHexMeshDict"), hdr("dictionary", "snappyHexMeshDict") + f"""
castellatedMesh true;
snap            true;
addLayers       false;
geometry {{ paddle.stl {{ type triSurfaceMesh; name paddle; }} }}
castellatedMeshControls
{{
    maxLocalCells 3000000; maxGlobalCells 6000000; minRefinementCells 10; maxLoadUnbalance 0.1; nCellsBetweenLevels 3;
    features ( {{ file "paddle.eMesh"; level 1; }} );
    refinementSurfaces {{ paddle {{ level (1 1); patchInfo {{ type wall; inGroups (wall); }} }} }}
    resolveFeatureAngle 30;
    refinementRegions {{ }}
    locationInMesh {fmt(lim)};
    allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 50; nRelaxIter 5; nFeatureSnapIter 10; implicitFeatureSnap false; explicitFeatureSnap true; multiRegionFeatureSnap false; }}
addLayersControls {{ {ALC} }}
meshQualityControls {{ {MQC} }}
writeFlags ( scalarLevels );
mergeTolerance 1e-6;
""")
CD_MESH = """
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
"""
FVS_MESH = """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
W(os.path.join(PM, "system", "controlDict"), hdr("dictionary", "controlDict") + CD_MESH)
W(os.path.join(PM, "system", "fvSchemes"), hdr("dictionary", "fvSchemes") + FVS_MESH)
W(os.path.join(PM, "system", "fvSolution"), hdr("dictionary", "fvSolution") + "\n")
W(os.path.join(PM, "system", "meshQualityDict"), hdr("dictionary", "meshQualityDict") + "\n")

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
    outlet {{ type patch; faces ((0 4 7 3)); }}
    sides  {{ type patch; faces ((0 1 5 4) (3 7 6 2)); }}
    bottom {{ type patch; faces ((0 3 2 1)); }}
    top    {{ type patch; faces ((4 5 6 7)); }}    // 滑移刚盖（单相，远离尾鳍 ≥ 12 cm，等价无限深水 = Flare 同条件）
);
mergePatchPairs ();
""")
W(os.path.join(CASE, "system", "snappyHexMeshDict"), hdr("dictionary", "snappyHexMeshDict") + f"""
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
    maxLocalCells 4000000; maxGlobalCells 8000000; minRefinementCells 10; maxLoadUnbalance 0.1; nCellsBetweenLevels 4;
    features ();
    refinementSurfaces {{ }}
    resolveFeatureAngle 30;
    refinementRegions {{ ref1 {{ mode inside; levels ((1e15 1)); }} ref2 {{ mode inside; levels ((1e15 2)); }} }}
    locationInMesh ({dom_min[0]+0.05:.4f} {dom_min[1]+0.05:.4f} {dom_min[2]+0.05:.4f});
    allowFreeStandingZoneFaces true;
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; }}
addLayersControls {{ {ALC.replace('expansionRatio 1.2', 'expansionRatio 1.0').replace('finalLayerThickness 0.4', 'finalLayerThickness 0.3')} }}
meshQualityControls {{ {MQC} }}
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
# ---- function objects：力（每步）+ 涡量/Q + 采样面（尾鳍表面、展向中截面、Q 等值面；vtk 二进制，小文件）
#      重叠网格：截面 / 等值面按 cellZone 分开取（c0 背景、c1 尾鳍部件），后处理时把背景那份在部件盒内的部分裁掉，避免重影
CHORD = float(pc[:, 0].max() - pc[:, 0].min()) if a.stl else float(a.naca.split(",")[1])
URel = a.U + vmax; QLEV = a.qstar * (URel / CHORD) ** 2
YMID = float(0.5 * (pc[:, 1].min() + pc[:, 1].max())) + Y0
t_samp = T / a.nsamp; t_s2 = max(0.0, (a.cycles - 2) * T); t_s1 = max(0.0, (a.cycles - 1) * T)
def samp(name, fields, surfs, tstart):
    return f"""    {name}
    {{
        type            surfaces;
        libs            ("libsampling.so");
        writeControl    adjustableRunTime;
        writeInterval   {t_samp:.6f};
        timeStart       {tstart:.6f};
        surfaceFormat   vtk;
        formatOptions   {{ vtk {{ format binary; }} }}
        interpolationScheme cellPoint;
        fields          ({fields});
        surfaces
        {{
{surfs}
        }}
        errors          warn;
    }}
"""
plane = lambda z: f"""            mid_{z} {{ type cuttingPlane; planeType pointAndNormal; pointAndNormalDict {{ point (0 {YMID:.6f} 0); normal (0 1 0); }} zones ({z}); interpolate false; }}"""
iso = lambda z: f"""            q_{z} {{ type isoSurface; isoField Q; isoValue {QLEV:.4g}; zones ({z}); interpolate true; }}"""
FUNCS = f"""functions
{{
    forces
    {{
        type            forces;
        libs            ("libforces.so");
        patches         (paddle);
        rho             rhoInf;
        rhoInf          {RHO};
        CofR            ({OX:.5f} 0 {OZ:.5f});     // 力矩参考点 = 运动表原点（腿 = 髋轴 O2）；铰链力矩在后处理里换算
        writeControl    timeStep;
        writeInterval   1;
        log             false;
    }}
    vorticity {{ type vorticity; libs ("libfieldFunctionObjects.so"); executeControl timeStep; writeControl writeTime; errors warn; }}
    Q         {{ type Q;         libs ("libfieldFunctionObjects.so"); executeControl timeStep; writeControl writeTime; errors warn; }}
    // 尾鳍表面压力（最后 2 个周期，每周期 {a.nsamp} 帧）→ 弦向 Cp、压力中心
{samp("finSurf", "p", "            fin { type patch; patches (paddle); interpolate false; }", t_s2)}    // 展向中截面 y = {YMID*1000:.1f} mm（最后 2 个周期）→ 尾迹、LEV、涡量云图 / 动画
{samp("midPlane", "U p vorticity", plane("c0") + chr(10) + plane("c1"), t_s2)}"""
if a.iso:
    FUNCS += f"""    // Q 等值面 Q = {QLEV:.4g} 1/s²（Q* = {a.qstar}，c = {CHORD*1000:.1f} mm，U_rel = {URel:.2f} m/s；最后 1 个周期）→ 三维涡结构
{samp("qIso", "vorticity p", iso("c0") + chr(10) + iso("c1"), t_s1)}"""
FUNCS += "}\n"
if a.functions_only:
    import re as _re
    cd = os.path.join(CASE, "system", "controlDict"); txt = open(cd).read(); j = txt.index("functions")
    if not os.path.exists(cd + ".bak"): shutil.copyfile(cd, cd + ".bak")
    open(cd, "w").write(txt[:j] + FUNCS)
    print(f"已替换 {cd} 的 functions 块（原文件 → controlDict.bak）：采样每周期 {a.nsamp} 帧，从 t = {t_s2:.3f} s 起；iso = {a.iso}（Q = {QLEV:.4g}）")
    sys.exit(0)

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
{FUNCS}""")
turbDiv = "" if a.turb == "laminar" else "\n    div(phi,k)      Gauss limitedLinear 1;\n    div(phi,omega)  Gauss limitedLinear 1;"
W(os.path.join(CASE, "system", "fvSchemes"), hdr("dictionary", "fvSchemes") + f"""
ddtSchemes {{ default Euler; }}
gradSchemes {{ default Gauss linear; }}
divSchemes
{{
    default         none;
    div(phi,U)      Gauss limitedLinearV 1;{turbDiv}
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default corrected; }}
// 洞切：搜索盒限制在运动包络内，体素 ≈ {vox_eff*1000:.2f} mm（必须明显小于最薄处 {a.tmin*1000:.1f} mm，否则全域变洞）
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
    best = (n, 1, 1)
    for nx_ in range(n, 0, -1):
        for ny_ in range(n, 0, -1):
            if n % (nx_ * ny_) == 0:
                nz_ = n // (nx_ * ny_)
                if (max(nx_, ny_, nz_), -nx_) < (max(best), -best[0]): best = (nx_, ny_, nz_)
    return best
hn = hier_n(a.np)
W(os.path.join(CASE, "system", "decomposeParDict"), hdr("dictionary", "decomposeParDict") + f"""
numberOfSubdomains {a.np};
method          hierarchical;           // scotch 在重叠网格上会触发 lduInterface 转型错误（v2606 实测），用 hierarchical
hierarchicalCoeffs {{ n ({hn[0]} {hn[1]} {hn[2]}); order xyz; }}
""")
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
        CofG            ({H0[0]:.6f} {H0[1]:.6f} {H0[2]:.6f});     // 转轴（尾鳍销轴）t=0 位置
        timeDataFileName "<constant>/6DoF.dat";
    }}
}}
""")
with open(os.path.join(CASE, "constant", "6DoF.dat"), "w") as fh:
    fh.write(f"{nT}\n(\n")
    for k in range(nT):
        fh.write(f"({t[k]:.6f} (({trans[k,0]:.6f} 0 {trans[k,2]:.6f}) (0 {ry[k]:.5f} 0)))\n")
    fh.write(")\n")
if a.limitU > 0:
    W(os.path.join(CASE, "system", "fvOptions"), hdr("dictionary", "fvOptions") + f"""
// 速度上限：尾鳍表面最大 {vmax:.2f} m/s + 来流 {a.U} m/s；重叠边界 / 薄片单元里的伪速度可到几 m/s，会把 Δt 卡在 0.1 ms
limitU {{ type limitVelocity; active yes; selectionMode all; max {a.limitU}; }}
""")
W(os.path.join(CASE, "constant", "transportProperties"), hdr("dictionary", "transportProperties", "constant") + f"\ntransportModel  Newtonian;\nnu              {NU};\n")
if a.turb == "laminar":
    W(os.path.join(CASE, "constant", "turbulenceProperties"), hdr("dictionary", "turbulenceProperties", "constant") + "\nsimulationType laminar;\n")
else:
    W(os.path.join(CASE, "constant", "turbulenceProperties"), hdr("dictionary", "turbulenceProperties", "constant") + "\nsimulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }\n")
Uin = f"(-{a.U:.4f} 0 0)" if a.U > 0 else "(0 0 0)"
W(os.path.join(CASE, "0.orig", "U"), hdr("volVectorField", "U", "0") + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform {Uin};
boundaryField
{{
    oversetPatch {{ type overset; }}
    paddle       {{ type movingWallVelocity; value uniform (0 0 0); }}
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
    paddle       { type zeroGradient; }
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
    for fld, dim, iv, wall in (("k", "[0 2 -2 0 0 0 0]", k_in, "kqRWallFunction"), ("omega", "[0 0 -1 0 0 0 0]", om_in, "omegaWallFunction")):
        W(os.path.join(CASE, "0.orig", fld), hdr("volScalarField", fld, "0") + f"""
dimensions {dim};
internalField uniform {iv:.3e};
boundaryField
{{
    oversetPatch {{ type overset; }}
    paddle {{ type {wall}; value uniform {iv:.3e}; }}
    inlet {{ type fixedValue; value uniform {iv:.3e}; }}
    outlet {{ type inletOutlet; inletValue uniform {iv:.3e}; value uniform {iv:.3e}; }}
    "(sides|bottom|top)" {{ type slip; }}
}}
""")
    W(os.path.join(CASE, "0.orig", "nut"), hdr("volScalarField", "nut", "0") + """
dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField { oversetPatch { type overset; } paddle { type nutkWallFunction; value uniform 0; } ".*" { type calculated; value uniform 0; } }
""")
W(os.path.join(CASE, "Allrun"), """#!/bin/sh
cd "${0%/*}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
( cd paddleMesh || exit 1
  runApplication blockMesh
  runApplication surfaceFeatureExtract
  runApplication snappyHexMesh -overwrite
  runApplication checkMesh )
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication mergeMeshes . paddleMesh -overwrite
runApplication topoSet
restore0Dir
runApplication setFields
runApplication -s merged checkMesh
runApplication decomposePar
[ -n "$MESH_ONLY" ] && exit 0            # run_fluke.sh：先只做网格，求解由它带洞切检查启动
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

mj = json.load(open(a.motion.rsplit(".", 1)[0] + ".json")) if os.path.exists(a.motion.rsplit(".", 1)[0] + ".json") else {}
json.dump(dict(kind="fluke", ramp=a.ramp, U=a.U, cycles=a.cycles, T=T, AREA=AREA, origin=[OX, OZ], H0=H0.tolist(), FIN0=float(FIN0), REF=REF,
               motion=os.path.abspath(a.motion), motion_meta=mj, **GEOM, turb=a.turb, bg=a.bg, comp=a.comp, vox=vox_eff, dt0=dt0, vmax=vmax),
          open(os.path.join(CASE, "gait.json"), "w"), indent=1, ensure_ascii=False)
np.savetxt(os.path.join(CASE, "kinematics_check.csv"), np.column_stack([t, H[:, 0], H[:, 2], FIN, ry]),
           delimiter=",", header="t,Hx,Hz,fin_deg,ry_deg", comments="", fmt="%.7f")
pmin, pmax = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
print(f"算例：{CASE}   U = {a.U} m/s   {a.cycles} 个周期 (T = {T:.4f} s) → endTime {endTime:.3f} s   参考面积 {AREA*1e4:.2f} cm²")
print(f"几何 STL：{len(tri)} 三角形，t=0 包围盒 x[{pmin[0]:.4f},{pmax[0]:.4f}] y[{pmin[1]:.4f},{pmax[1]:.4f}] z[{pmin[2]:.4f},{pmax[2]:.4f}]")
print(f"部件网格盒 {n_x}×{n_y}×{n_z} = {n_x*n_y*n_z} 单元（{a.comp*1000:.2f} mm，表面再细化 1 级）")
print(f"运动包络 x[{env_min[0]:.3f},{env_max[0]:.3f}] y[{env_min[1]:.3f},{env_max[1]:.3f}] z[{env_min[2]:.3f},{env_max[2]:.3f}]")
print(f"背景域 x[{dom_min[0]:.3f},{dom_max[0]:.3f}] y[{dom_min[1]:.3f},{dom_max[1]:.3f}] z[{dom_min[2]:.3f},{dom_max[2]:.3f}]  blockMesh {nx}×{ny}×{nz}，加密到 {a.bg*1000:.1f} mm")
print(f"表面最大速度 {vmax:.2f} m/s → 初始 Δt = {dt0} s，maxCo {a.maxCo}；每周期写 {a.writes} 帧")
print(f"洞切搜索盒 {fmt(sb_min)}–{fmt(sb_max)}，体素 {sbd[0]}×{sbd[1]}×{sbd[2]}（≈{vox_eff*1000:.2f} mm）")
print("跑：cd", CASE, "&& ./Allrun")
