#!/usr/bin/env pvbatch
"""
整机 overset 算例可视化（pvbatch）：
    VIEW=3d    透视：船体（灰）+ 四块桨（压力着色）+ 自由液面（alpha.water = 0.5 等值面，半透明）+ 水下 |U| 切面
    VIEW=side  侧视 y = +y_leg 切面（左腿平面）|U| + 桨 + 船体轮廓
    VIEW=top   俯视：自由液面高程着色（看兴波与桨拍出的水花）
    VIEW=rear  艉视透视：从船尾后下方看四块桨与尾迹（|U| 切面 z = 水线 − 5 cm + 自由面）
    （full_video.sh <case> quad 会把 3d / side / top / rear 四个视角拼成一个 2×2 视频）
    xvfb-run -a pvbatch full_video.py <case> <frames_dir>     （用 full_video.sh 调）
环境变量：VIEW, UMAX(0.5), RES(1600x900), SKIP(1)
"""
import os, sys, json
from paraview.simple import *

case, outdir = sys.argv[1], sys.argv[2]
VIEW = os.environ.get("VIEW", "3d"); UMAX = float(os.environ.get("UMAX", 0.5))
RES = [int(v) for v in os.environ.get("RES", "1600x900").split("x")]; SKIP = int(os.environ.get("SKIP", 1))
os.makedirs(outdir, exist_ok=True)
gait = json.load(open(os.path.join(case, "gait.json"))) if os.path.exists(os.path.join(case, "gait.json")) else {}
T = gait.get("T", 1.25); DUTY = gait.get("DUTY", 0.5); fan = gait.get("fan", "?"); U0 = gait.get("U", 0.0)
WL = gait.get("WL", -0.010); FS = gait.get("free_surface", True); legs = gait.get("legs", {})
y_leg = abs(next(iter(legs.values()))["y"]) if legs else 0.075

def set_if(obj, name, val):
    try:
        if name in obj.ListProperties(): setattr(obj, name, val); return True
    except Exception: pass
    return False

foam = os.path.join(case, "case.foam"); open(foam, "a").close()
r = OpenFOAMReader(FileName=foam); r.MeshRegions = ["internalMesh"]
r.CellArrays = ["U", "p_rgh" if FS else "p", "cellTypes", "zoneID"] + (["alpha.water"] if FS else [])
set_if(r, "Createcelltopointfiltereddata", 1); set_if(r, "CreateCellToPoint", 1); set_if(r, "Decomposepolyhedra", 0)
r.UpdatePipelineInformation(); times = list(r.TimestepValues) if len(r.TimestepValues) else [0.0]
r.UpdatePipeline(times[-1])
have_zone = "zoneID" in list(r.CellData.keys())          # 算例若每步写了 zoneID（writeObjects），自由面只取背景网格，避免部件盒处的重叠伪影
try: avail = list(r.GetProperty("MeshRegions").GetAvailable())
except Exception: avail = []
pads = [x for x in avail if x.startswith("patch/paddle_")]; hulls = [x for x in avail if x == "patch/hull"]
pad = OpenFOAMReader(FileName=foam); pad.MeshRegions = pads or ["patch/paddle_FL"]; pad.CellArrays = ["p_rgh" if FS else "p"]
set_if(pad, "Createcelltopointfiltereddata", 1); set_if(pad, "CreateCellToPoint", 1)
hull = None
if hulls:
    hull = OpenFOAMReader(FileName=foam); hull.MeshRegions = hulls; hull.CellArrays = []

th0 = Threshold(Input=r); th0.Scalars = ["CELLS", "cellTypes"]
if not (set_if(th0, "LowerThreshold", -0.5) and set_if(th0, "UpperThreshold", 1.5)): th0.ThresholdRange = [-0.5, 1.5]
set_if(th0, "ThresholdMethod", "Between")
# 只画机器人附近（远场是 16 mm 粗网格，自由面在那里分辨不出来）：x ∈ [-0.45, 0.30], y ∈ ±0.25
def box_clip(src, xr=(-0.45, 0.30), yr=(-0.25, 0.25)):
    cx = Calculator(Input=src); cx.ResultArrayName = "cx"; cx.Function = "coordsX"
    tx = Threshold(Input=cx); tx.Scalars = ["POINTS", "cx"]
    if not (set_if(tx, "LowerThreshold", xr[0]) and set_if(tx, "UpperThreshold", xr[1])): tx.ThresholdRange = list(xr)
    set_if(tx, "ThresholdMethod", "Between")
    cy = Calculator(Input=tx); cy.ResultArrayName = "cy"; cy.Function = "coordsY"
    ty = Threshold(Input=cy); ty.Scalars = ["POINTS", "cy"]
    if not (set_if(ty, "LowerThreshold", yr[0]) and set_if(ty, "UpperThreshold", yr[1])): ty.ThresholdRange = list(yr)
    set_if(ty, "ThresholdMethod", "Between")
    return ty
th = box_clip(th0)
if have_zone:                                             # 背景网格（zoneID = 0）：画自由面用
    bgz = Threshold(Input=th); bgz.Scalars = ["CELLS", "zoneID"]
    if not (set_if(bgz, "LowerThreshold", -0.5) and set_if(bgz, "UpperThreshold", 0.5)): bgz.ThresholdRange = [-0.5, 0.5]
    set_if(bgz, "ThresholdMethod", "Between")
else: bgz = th

view = GetActiveViewOrCreate("RenderView"); view.ViewSize = RES; view.Background = [0.96, 0.97, 0.98]
set_if(view, "UseColorPaletteForBackground", 0); set_if(view, "OrientationAxesVisibility", 0)
lutU = GetColorTransferFunction("U"); lutU.ApplyPreset("Viridis (matplotlib)", True)
set_if(lutU, "AutomaticRescaleRangeMode", "Never"); lutU.RescaleTransferFunction(0.0, UMAX)
pname = "p_rgh" if FS else "p"
lutP = GetColorTransferFunction(pname); lutP.ApplyPreset("Cool to Warm", True)

dP = Show(pad, view); ColorBy(dP, ("POINTS", pname)); dP.LookupTable = lutP; dP.SetScalarBarVisibility(view, False)
if hull is not None:
    dH = Show(hull, view)
    try: ColorBy(dH, None)
    except Exception: dH.ColorArrayName = ["POINTS", ""]
    dH.DiffuseColor = [0.78, 0.78, 0.76]; dH.AmbientColor = dH.DiffuseColor; dH.Opacity = 0.85 if VIEW == "3d" else 1.0

if VIEW == "side":
    sl = Slice(Input=th); sl.SliceType = "Plane"; sl.SliceType.Origin = [0, y_leg, 0]; sl.SliceType.Normal = [0, 1, 0]
    dS = Show(sl, view); ColorBy(dS, ("POINTS", "U", "Magnitude")); dS.LookupTable = lutU
    bar = GetScalarBar(lutU, view); bar.Title = "|U| [m/s]"; bar.ComponentTitle = ""; bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.WindowLocation = "Any Location"; bar.Position = [0.88, 0.25]; bar.ScalarBarLength = 0.5; dS.SetScalarBarVisibility(view, True)
    if FS:
        fs = Contour(Input=bgz); fs.ContourBy = ["POINTS", "alpha.water"]; fs.Isosurfaces = [0.5]
        dF = Show(fs, view); dF.DiffuseColor = [0.3, 0.55, 0.9]; dF.Opacity = 0.35
        try: ColorBy(dF, None)
        except Exception: dF.ColorArrayName = ["POINTS", ""]
elif VIEW == "top":
    if FS:
        fs = Contour(Input=bgz); fs.ContourBy = ["POINTS", "alpha.water"]; fs.Isosurfaces = [0.5]
        el = Calculator(Input=fs); el.ResultArrayName = "eta"; el.Function = f"coordsZ-({WL})"
        dF = Show(el, view); ColorBy(dF, ("POINTS", "eta")); lutE = GetColorTransferFunction("eta"); lutE.ApplyPreset("Cool to Warm", True)
        set_if(lutE, "AutomaticRescaleRangeMode", "Never"); lutE.RescaleTransferFunction(-0.008, 0.008); dF.LookupTable = lutE
        bar = GetScalarBar(lutE, view); bar.Title = "wave elevation [m]"; bar.ComponentTitle = ""; bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]; dF.SetScalarBarVisibility(view, True)
else:  # 3d / rear
    sl = Slice(Input=th); sl.SliceType = "Plane"; sl.SliceType.Origin = [0, 0, WL - 0.05]; sl.SliceType.Normal = [0, 0, 1]
    dS = Show(sl, view); ColorBy(dS, ("POINTS", "U", "Magnitude")); dS.LookupTable = lutU; dS.Opacity = 0.7
    bar = GetScalarBar(lutU, view); bar.Title = "|U| [m/s]"; bar.ComponentTitle = ""; bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.WindowLocation = "Any Location"; bar.Position = [0.88, 0.2]; bar.ScalarBarLength = 0.45; dS.SetScalarBarVisibility(view, True)
    if FS:
        fs = Contour(Input=bgz); fs.ContourBy = ["POINTS", "alpha.water"]; fs.Isosurfaces = [0.5]
        dF = Show(fs, view); dF.DiffuseColor = [0.35, 0.6, 0.95]; dF.Opacity = 0.35
        try: ColorBy(dF, None)
        except Exception: dF.ColorArrayName = ["POINTS", ""]

txt = Text(); dT = Show(txt, view); dT.Color = [0.1, 0.1, 0.1]; dT.FontSize = 20; dT.WindowLocation = "Upper Left Corner"
CAM = {"3d":   dict(pos=[0.36, -0.42, 0.24], foc=[-0.04, 0.0, -0.04], up=[0, 0, 1], par=0, scale=0.2),
       "side": dict(pos=[-0.06, -1.2, -0.03], foc=[-0.06, 0, -0.03], up=[0, 0, 1], par=1, scale=0.115),
       "top":  dict(pos=[-0.07, 0, 1.2], foc=[-0.07, 0, WL], up=[0, 1, 0], par=1, scale=0.135),
       "rear": dict(pos=[-0.60, -0.26, -0.10], foc=[-0.04, 0.0, -0.04], up=[0, 0, 1], par=0, scale=0.2)}[VIEW]
def apply_cam():
    view.CameraPosition = CAM["pos"]; view.CameraFocalPoint = CAM["foc"]; view.CameraViewUp = CAM["up"]
    view.CameraParallelProjection = CAM["par"]; view.CameraParallelScale = CAM["scale"]
    if not CAM["par"]: view.CameraViewAngle = 28
Render(view)                       # 第一次 Render 会自动重置相机，先渲染一次再逐帧设相机

n = 0
for i, t in enumerate(times):
    if i % SKIP: continue
    view.ViewTime = t; tau = (t / T) % 1.0
    ph = " ".join(f"{L}:{'P' if ((tau + v['phase']) % 1.0) < DUTY else 'R'}" for L, v in legs.items())
    txt.Text = f"t = {t:6.3f} s   t/T = {t/T:5.2f}   {ph}   fan: {fan}   U = {U0} m/s"
    apply_cam(); lutU.RescaleTransferFunction(0.0, UMAX); Render(view)
    SaveScreenshot(os.path.join(outdir, f"frame.{n:04d}.png"), view, ImageResolution=RES); n += 1
print(f"{n} frames -> {outdir}")
