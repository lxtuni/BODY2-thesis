#!/usr/bin/env pvbatch
"""
单腿 overset 算例可视化（pvbatch）：
    VIEW=side  y=0 切面 |U| + 桨板（俯视机构平面，看划水/回收的尾迹）
    VIEW=3d    透视：桨板着色压力 + |U| 切面 + 涡量等值面（看脚蹼张开态的三维涡）
    xvfb-run -a pvbatch paddle_video.py <case> <frames_dir>        （用 paddle_video.sh 调）
环境变量：VIEW, UMAX(0.5), VORT(等值面涡量 1/s, 默认 40), RES(1400x900), SKIP(每几帧取一帧, 默认 1)
          LEG=1 叠画整条腿的连杆骨架（摇臂/曲柄/从动杆/销轴，按运动学同步；桨板是 CFD 里的真实物体），WL 水面参照高度
"""
import os, sys, json, math
from paraview.simple import *

case, outdir = sys.argv[1], sys.argv[2]
VIEW = os.environ.get("VIEW", "side"); UMAX = float(os.environ.get("UMAX", 0.5)); VORT = float(os.environ.get("VORT", 40))
RES = [int(v) for v in os.environ.get("RES", "1400x900").split("x")]; SKIP = int(os.environ.get("SKIP", 1))
os.makedirs(outdir, exist_ok=True)
gait = json.load(open(os.path.join(case, "gait.json"))) if os.path.exists(os.path.join(case, "gait.json")) else {}
T = gait.get("T", 1.25); DUTY = gait.get("DUTY", 0.5); fan = gait.get("fan", "?")
O2 = (0.0634, 0.0, 0.004)
LEG = int(os.environ.get("LEG", 1))          # 1 = 叠画整条腿的连杆骨架（按运动学同步），0 = 只画桨板
WL_Z = float(os.environ.get("WL", -0.010))   # 水线（画一个半透明面做参照；单相算例里它不是计算边界）
ld = None
if LEG:
    try:
        for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W"]:
            if k in gait: os.environ[k] = str(gait[k])
        if "psi0_deg" in gait: os.environ["PSI_START"] = str(gait["psi0_deg"])
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import importlib.util
        spec = importlib.util.spec_from_file_location("ld", os.path.join(os.path.dirname(os.path.abspath(__file__)), "leg_dynamics.py"))
        ld = importlib.util.module_from_spec(spec); spec.loader.exec_module(ld)
    except Exception as e:
        print("连杆叠画不可用（leg_dynamics.py 导入失败）:", e); ld = None

def set_if(obj, name, val):
    try:
        if name in obj.ListProperties(): setattr(obj, name, val); return True
    except Exception: pass
    return False

foam = os.path.join(case, "case.foam"); open(foam, "a").close()
r = OpenFOAMReader(FileName=foam); r.MeshRegions = ["internalMesh"]; r.CellArrays = ["U", "p", "cellTypes", "vorticity"]
set_if(r, "Createcelltopointfiltereddata", 1); set_if(r, "CreateCellToPoint", 1); set_if(r, "Decomposepolyhedra", 0)
r.UpdatePipelineInformation(); times = list(r.TimestepValues) if len(r.TimestepValues) else [0.0]
pad = OpenFOAMReader(FileName=foam); pad.CellArrays = ["p", "U"]
try: avail = list(pad.GetProperty("MeshRegions").GetAvailable())
except Exception: avail = ["patch/paddle"]
pad.MeshRegions = [x for x in ("patch/paddle", "patch/links") if x in avail] or ["patch/paddle"]   # --leg 算例带 links 壁面
set_if(pad, "Createcelltopointfiltereddata", 1); set_if(pad, "CreateCellToPoint", 1)

# 去掉洞单元（cellTypes: 0 计算, 1 插值, 2 洞）
th = Threshold(Input=r); th.Scalars = ["CELLS", "cellTypes"]
if not (set_if(th, "LowerThreshold", -0.5) and set_if(th, "UpperThreshold", 1.5)): th.ThresholdRange = [-0.5, 1.5]
set_if(th, "ThresholdMethod", "Between")
sl = Slice(Input=th); sl.SliceType = "Plane"; sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [0, 1, 0]
set_if(sl, "Triangulatetheslice", 0)

view = GetActiveViewOrCreate("RenderView"); view.ViewSize = RES; view.Background = [0.96, 0.97, 0.98]
set_if(view, "UseColorPaletteForBackground", 0); set_if(view, "OrientationAxesVisibility", 0)
lutU = GetColorTransferFunction("U"); lutU.ApplyPreset("Viridis (matplotlib)", True)
set_if(lutU, "AutomaticRescaleRangeMode", "Never"); lutU.RescaleTransferFunction(0.0, UMAX)
lutP = GetColorTransferFunction("p"); lutP.ApplyPreset("Cool to Warm", True)

dS = Show(sl, view); ColorBy(dS, ("POINTS", "U", "Magnitude")); dS.LookupTable = lutU
bar = GetScalarBar(lutU, view); bar.Title = "|U| [m/s]"; bar.ComponentTitle = ""; bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
bar.WindowLocation = "Any Location"; bar.Position = [0.88, 0.25]; bar.ScalarBarLength = 0.5; dS.SetScalarBarVisibility(view, True)
dP = Show(pad, view)
if VIEW == "3d":
    ColorBy(dP, ("POINTS", "p")); dP.LookupTable = lutP; dP.SetScalarBarVisibility(view, False)
    dS.Opacity = 0.6
    # 涡量：优先用求解器写出的 vorticity 场（controlDict 里的 vorticity 函数对象），没有才用 ParaView 算
    src = th
    r.UpdatePipeline(times[-1])
    have_vort = "vorticity" in list(r.CellData.keys())
    if not have_vort:
        src = Gradient(Input=th) if "Gradient" in dir() else GradientOfUnstructuredDataSet(Input=th)
        set_if(src, "ScalarArray", ["POINTS", "U"]); set_if(src, "ComputeVorticity", 1); set_if(src, "VorticityArrayName", "vorticity")
    calc = Calculator(Input=src); calc.ResultArrayName = "vortMag"; calc.Function = "mag(vorticity)"
    iso = Contour(Input=calc); iso.ContourBy = ["POINTS", "vortMag"]; iso.Isosurfaces = [VORT]
    dI = Show(iso, view); ColorBy(dI, ("POINTS", "U", "Magnitude")); dI.LookupTable = lutU; dI.Opacity = 0.55
else:
    try: ColorBy(dP, None)
    except Exception: dP.ColorArrayName = ["POINTS", ""]
    dP.DiffuseColor = [0.55, 0.16, 0.08]; dP.AmbientColor = dP.DiffuseColor

txt = Text(); dT = Show(txt, view); dT.Color = [0.1, 0.1, 0.1]; dT.FontSize = 20; dT.WindowLocation = "Upper Left Corner"

# ---- 整条腿的连杆骨架（3D 圆管），每帧按运动学更新 ----
import numpy as np
def hexrgb(h): return [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)]
LINKS = {   # 名称: (颜色, 半径, y 偏移)  —— 全部放在 y<0 一侧（相机在 -y），否则会被 y=0 的切面挡住
    "crank1":   ("#fd8d3c", 0.0045, -0.010),
    "rod":      ("#756bb1", 0.0025, -0.008),
    "rocker":   ("#31a354", 0.0045, -0.006),
    "crank2":   ("#e6550d", 0.0045, -0.004),
    "follower": ("#636363", 0.0028, -0.002),
}
PIN_Y = -0.006
tubes, pins = {}, {}
if ld is not None:
    for name, (col, rad, yo) in LINKS.items():
        pl = PolyLineSource(Points=[0, yo, 0, 0.01, yo, 0]); tb = Tube(Input=pl); tb.Radius = rad
        set_if(tb, "NumberofSides", 16)
        d_ = Show(tb, view); d_.DiffuseColor = hexrgb(col); d_.AmbientColor = hexrgb(col)
        try: ColorBy(d_, None)
        except Exception: d_.ColorArrayName = ["POINTS", ""]
        tubes[name] = (pl, yo)
    for name in ("O1", "O2", "A", "B", "C", "D", "E"):
        sp = Sphere(Radius=0.0032); set_if(sp, "ThetaResolution", 16); set_if(sp, "PhiResolution", 16)
        d_ = Show(sp, view); d_.DiffuseColor = [0.12, 0.12, 0.12]; d_.AmbientColor = [0.12, 0.12, 0.12]
        try: ColorBy(d_, None)
        except Exception: d_.ColorArrayName = ["POINTS", ""]
        pins[name] = sp
    # 水面参照
    wlp = Plane(); wlp.Origin = [-0.25, -0.2, WL_Z]; wlp.Point1 = [0.35, -0.2, WL_Z]; wlp.Point2 = [-0.25, 0.2, WL_Z]
    dW = Show(wlp, view); dW.DiffuseColor = [0.55, 0.75, 0.95]; dW.Opacity = 0.18
    try: ColorBy(dW, None)
    except Exception: dW.ColorArrayName = ["POINTS", ""]

def update_leg(t):
    """按时间 t 更新连杆与销轴位置（机构平面 y=0）"""
    if ld is None: return
    psi, phi, _ = ld.gait(np.array([t])); P = ld.pose(float(psi[0]), float(phi[0]))
    O1 = ld.O1; O2_ = ld.O2
    th1 = ld.crank1_from_rocker(float(psi[0]) + np.pi)
    A = O1 + ld.L_CRANK1 * ld.rot(th1) if np.isfinite(th1) else P["B"]
    J = {"O1": O1, "O2": O2_, "A": A, "B": P["B"], "C": P["C"], "D": P["D"], "E": P["E"]}
    segs = {"crank1": ("O1", "A"), "rod": ("A", "B"), "rocker": ("B", "C"), "crank2": ("O2", "D"), "follower": ("C", "E")}
    for name, (a, b) in segs.items():
        pl, yo = tubes[name]; pa, pb = J[a], J[b]
        pl.Points = [float(pa[0]), yo, float(pa[1]), float(pb[0]), yo, float(pb[1])]
    for name, sp in pins.items():
        p = J[name]; sp.Center = [float(p[0]), PIN_Y, float(p[1])]

# 相机：以髋轴 O2 下方的扫掠区为中心
fx, fz = O2[0] - 0.005, O2[2] - 0.045
CAM = {"side": dict(pos=[fx, -1.0, fz], foc=[fx, 0, fz], up=[0, 0, 1], par=1, scale=0.095),
       "3d":   dict(pos=[fx + 0.14, -0.26, fz + 0.11], foc=[fx, 0, fz], up=[0, 0, 1], par=0, scale=0.1)}[VIEW]
def apply_cam():
    view.CameraPosition = CAM["pos"]; view.CameraFocalPoint = CAM["foc"]; view.CameraViewUp = CAM["up"]
    view.CameraParallelProjection = CAM["par"]; view.CameraParallelScale = CAM["scale"]
    if not CAM["par"]: view.CameraViewAngle = 30

Render(view)                       # 第一次 Render 会自动重置相机，先渲染一次再逐帧设相机
n = 0
for i, t in enumerate(times):
    if i % SKIP: continue
    view.ViewTime = t
    tau = (t / T) % 1.0 if T > 0 else 0
    phase = "POWER stroke" if tau < DUTY else "RECOVERY"
    txt.Text = f"t = {t:6.3f} s   t/T = {t/T:5.2f}   {phase}   fan: {fan}   U = {gait.get('U',0)} m/s"
    update_leg(t)
    apply_cam(); lutU.RescaleTransferFunction(0.0, UMAX); Render(view)
    SaveScreenshot(os.path.join(outdir, f"frame.{n:04d}.png"), view, ImageResolution=RES); n += 1
print(f"{n} frames -> {outdir}")
