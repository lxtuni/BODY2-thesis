#!/usr/bin/env pvbatch
"""
裸船体算例的动画帧（ParaView 批处理）—— 供 make_hull_video.sh 调用
    pvbatch make_hull_video.py <算例目录> <帧输出目录>
视角 VIEW=...（环境变量）：
    side   船体中纵剖面 y = 0 的竖直切面 |U|（只画水相）+ 自由液面线 + 船体，正侧视（正交投影，+x 艏在右）   ← 默认
    top    自由液面波高（俯视，正交投影）+ 船体甲板
其他环境变量：
    U=0.2           来流速度（决定 |U| 色标上限与波高色标范围的默认值）
    UMAX=           |U| 色标上限 m/s（默认 1.6·U）
    BAND=           波高着色范围 ±BAND m（默认 0.8·U²/g）
    WL=             水线 m（默认读 constant/hRef）
    RES=1280x380    分辨率（top 视角默认 640x380）
    SKIP=1          跳过前几帧（t = 0 的初始场没有信息）
    DECOMPOSED=1    直接读 processor* 分块结果（没 reconstruct 时）
    LABEL=          左上角文字（ASCII；中文标签由 ffmpeg drawtext 加）
"""
import os, re, sys
from paraview.simple import *  # noqa

case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd()
out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(case, "frames")
os.makedirs(out, exist_ok=True)

def read_href(case):
    try:
        m = re.search(r"^\s*value\s+([-+0-9.eE]+)\s*;", open(os.path.join(case, "constant", "hRef")).read(), re.M)
        return float(m.group(1)) if m else None
    except OSError:
        return None

def env_float(name, default):
    v = os.environ.get(name, ""); return float(v) if v.strip() else default

VIEW = os.environ.get("VIEW", "side").lower()
U0 = env_float("U", 0.2)
WL = env_float("WL", read_href(case) if read_href(case) is not None else -0.010)
UMAX = env_float("UMAX", 1.6 * U0)
BAND = env_float("BAND", 0.8 * U0 ** 2 / 9.81)
RES = [int(v) for v in os.environ.get("RES", "640x380" if VIEW == "top" else "1280x380").lower().split("x")]
SKIP = int(os.environ.get("SKIP", 1))
DECOMPOSED = os.environ.get("DECOMPOSED", "0") == "1"
LABEL = os.environ.get("LABEL", "")
INK = [0.06, 0.13, 0.17]

foam = os.path.join(case, "case.foam"); open(foam, "a").close()

def set_if_exists(proxy, names, value):
    for nm in names:
        if nm in proxy.ListProperties():
            setattr(proxy, nm, value); return True
    return False

# ---------------- 读数据 ----------------
rWater = OpenFOAMReader(FileName=foam)
if DECOMPOSED: set_if_exists(rWater, ["CaseType"], "Decomposed Case")
rWater.MeshRegions = ["internalMesh"]; rWater.CellArrays = ["alpha.water", "U"]
set_if_exists(rWater, ["CreateCellToPoint", "Createcelltopointfiltereddata"], 1)
avail = list(rWater.MeshRegions.Available)
hull_region = next((r for r in avail if r == "hull" or r.endswith("/hull")), None)
if hull_region is None:
    print("!! 找不到船体 patch（hull），可用区域：", avail); sys.exit(1)
rHull = OpenFOAMReader(FileName=foam)
if DECOMPOSED: set_if_exists(rHull, ["CaseType"], "Decomposed Case")
rHull.MeshRegions = [hull_region]; rHull.CellArrays = ["U"]
times = [t for t in rWater.TimestepValues]
if not times:
    print("!! 没有时间步 —— 算例还没跑，或者还没 reconstruct"); sys.exit(1)
water = MergeBlocks(Input=rWater); hull = MergeBlocks(Input=rHull)
hull.UpdatePipeline(times[-1]); b = hull.GetDataInformation().GetBounds()
cx = 0.5 * (b[0] + b[1]); L = max(b[1] - b[0], 1e-3); zmin = b[4]
print(f"算例 {case}\n{len(times)} 帧  视角 {VIEW}  船长 {L:.3f} m  龙骨 {zmin:.4f} m  水线 {WL} m  |U|max {UMAX:.3f}  波高 ±{BAND*1e3:.1f} mm")

# 感兴趣区域（船体偏右，尾流在左）：侧视 艏前 0.25 L / 艉后 0.75 L；俯视 艏前 0.4 L / 艉后 1.0 L
X0, X1 = (cx - 1.5 * L, cx + 0.9 * L) if VIEW == "top" else (cx - 1.25 * L, cx + 0.75 * L)
def clip_box(src, lo, hi):
    c = Clip(Input=src); c.ClipType = "Box"; c.ClipType.Position = lo; c.ClipType.Length = [hi[i] - lo[i] for i in range(3)]
    set_if_exists(c, ["Invert"], 1); set_if_exists(c, ["Crinkleclip", "Crinkle clip"], 0); return c

def water_only(src):
    t = Threshold(Input=src); t.Scalars = ["POINTS", "alpha.water"]
    if set_if_exists(t, ["LowerThreshold"], 0.5):
        set_if_exists(t, ["UpperThreshold"], 1.5); set_if_exists(t, ["ThresholdMethod"], "Between")
    else:
        t.ThresholdRange = [0.5, 1.5]
    return t

def style_bar(bar, title, pos, horizontal=False, length=0.33):
    bar.Title = title; bar.ComponentTitle = ""; bar.TitleFontSize = 13; bar.LabelFontSize = 11
    bar.ScalarBarLength = length; bar.TitleColor = INK; bar.LabelColor = INK; bar.RangeLabelFormat = "%-#6.3g"
    if horizontal:
        try: bar.Orientation = "Horizontal"
        except Exception: pass
    try: bar.WindowLocation = pos
    except Exception: pass

view = GetActiveViewOrCreate("RenderView"); view.ViewSize = RES
try: view.UseColorPaletteForBackground = 0
except Exception: pass
view.Background = [0.955, 0.96, 0.955]; view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1

def show_hull(view):
    d = Show(hull, view)
    try: ColorBy(d, None)
    except Exception: d.ColorArrayName = ["POINTS", ""]
    d.DiffuseColor = [0.42, 0.44, 0.46]; d.AmbientColor = [0.42, 0.44, 0.46]; d.Specular = 0.35
    return d

if VIEW == "top":
    fs = Contour(Input=water, ContourBy=["POINTS", "alpha.water"], Isosurfaces=[0.5]); fs.ComputeScalars = 0
    el = Calculator(Input=fs, ResultArrayName="elevation", Function=f"(coordsZ-({WL}))*1000")      # 相对静水线的波高 [mm]
    fsc = clip_box(el, [X0, -0.9 * L, WL - 0.2 * L], [X1, 0.9 * L, WL + 0.2 * L])
    d = Show(fsc, view); ColorBy(d, ("POINTS", "elevation"))
    lut = GetColorTransferFunction("elevation"); lut.ApplyPreset("Cool to Warm", True); lut.RescaleTransferFunction(-BAND * 1e3, BAND * 1e3)
    d.SetScalarBarVisibility(view, True); style_bar(GetScalarBar(lut, view), "wave elevation [mm]", "Lower Center", horizontal=True, length=0.45)
    show_hull(view)
    W = 0.5 * (X1 - X0); Hh = W * RES[1] / RES[0]                 # 视野半宽 / 半高
    CAM = dict(focal=[0.5 * (X0 + X1), 0.0, WL], pos=[0.5 * (X0 + X1), 0.0, WL + 3 * L], up=[0, 1, 0], scale=Hh)
else:
    sl = Slice(Input=water); sl.SliceType = "Plane"; sl.SliceType.Origin = [cx, 0.0, WL]; sl.SliceType.Normal = [0, 1, 0]
    slw = clip_box(water_only(sl), [X0, -0.1, zmin - 0.4 * L], [X1, 0.1, WL + 0.25 * L])
    d = Show(slw, view); ColorBy(d, ("POINTS", "U", "Magnitude"))
    lut = GetColorTransferFunction("U"); lut.ApplyPreset("Viridis (matplotlib)", True); lut.VectorMode = "Magnitude"; lut.RescaleTransferFunction(0.0, UMAX)
    d.SetScalarBarVisibility(view, True); style_bar(GetScalarBar(lut, view), "|U| [m/s]", "Lower Right Corner", horizontal=True, length=0.22)
    # 自由液面线：切面上 alpha = 0.5 的等值线
    fl = Contour(Input=clip_box(sl, [X0, -0.1, zmin - 0.4 * L], [X1, 0.1, WL + 0.25 * L]), ContourBy=["POINTS", "alpha.water"], Isosurfaces=[0.5]); fl.ComputeScalars = 0
    dl = Show(fl, view)
    try: ColorBy(dl, None)
    except Exception: dl.ColorArrayName = ["POINTS", ""]
    dl.DiffuseColor = [0.05, 0.25, 0.55]; dl.AmbientColor = dl.DiffuseColor; dl.LineWidth = 3.0; dl.RenderLinesAsTubes = 1
    # 静水线参考（细虚线用一条细线代替）
    ref = Line(); ref.Point1 = [X0, 0.0, WL]; ref.Point2 = [X1, 0.0, WL]; ref.Resolution = 1
    dr = Show(ref, view); dr.DiffuseColor = [0.55, 0.55, 0.55]; dr.AmbientColor = dr.DiffuseColor; dr.LineWidth = 1.0
    show_hull(view)
    W = 0.5 * (X1 - X0); Hh = W * RES[1] / RES[0]; zc = (WL + 0.18 * L) - Hh          # 上边留 0.18 L 的空气放标签，其余给水下
    CAM = dict(focal=[0.5 * (X0 + X1), 0.0, zc], pos=[0.5 * (X0 + X1), -3 * L, zc], up=[0, 0, 1], scale=Hh)

def apply_camera():      # 第一次 Render 会自动 ResetCamera，所以每帧渲染前都重新设置相机
    view.CameraParallelProjection = 1; view.CameraFocalPoint = CAM["focal"]; view.CameraPosition = CAM["pos"]
    view.CameraViewUp = CAM["up"]; view.CameraParallelScale = CAM["scale"]

if LABEL:
    label = Text(Text=LABEL); dLabel = Show(label, view); dLabel.FontSize = 15; dLabel.Color = INK
    try: dLabel.WindowLocation = "Upper Left Corner"
    except Exception: pass
step = Text(Text=""); dStep = Show(step, view); dStep.FontSize = 13; dStep.Color = INK
try: dStep.WindowLocation = "Upper Right Corner" if VIEW == "top" else "Upper Left Corner"
except Exception: pass

n = 0
for i, t in enumerate(times):
    if i < SKIP: continue
    view.ViewTime = t; step.Text = f"LTS step {t:g} / {times[-1]:g}" if VIEW == "top" else f"LTS step {t:g} / {times[-1]:g}     flow -x, bow at right"
    apply_camera(); Render(view); apply_camera(); Render(view); SaveScreenshot(os.path.join(out, f"frame.{n:04d}.png"), view, ImageResolution=RES); n += 1
    print(f"  帧 {n:3d}  t = {t:g}")
print(f"完成：{n} 帧 → {out}")
