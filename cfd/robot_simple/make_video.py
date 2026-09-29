#!/usr/bin/env pvbatch
"""
CFD 结果动画帧生成器（ParaView 批处理，不用点鼠标）
==================================================
用法（Ubuntu 终端，算例已跑完并 reconstruct）：
    pvbatch make_video.py <算例目录> [帧输出目录]

视角 VIEW=...（环境变量）：
    flow      水下总览：桨板平面的竖直切面 |U| + 流线 + 半透明水面 + 船体，侧前方 3/4 视角   ← 默认
    slice     竖直切面 |U| 的正侧视图（正交投影，定量看尾流/边界层）
    topuw     桨板深度的水平切面 |U|，俯视（看每片桨板后面的尾流）
    pressure  船体表面动压 p_rgh + 半透明水面，前下方视角（看驻点和低压区）
    surface   自由液面波高（俯前方视角）

其他环境变量：
    WL=0.012        水线高度 m（默认读 constant/hRef）
    BAND=0.004      波高着色范围 ±BAND m
    UMAX=           |U| 着色上限 m/s（默认取数据里的最大值）
    YSLICE=         竖直切面的 y 位置（默认 -0.39×船宽，正好切过桨板）
    ZSLICE=         水平切面的 z 位置（默认在水下部分 30% 深处）
    ROI=1.0         画面裁剪范围（倍船长）
    SKIP=0          跳过前几帧（起步冲击）
    RES=1920x1080   分辨率
    DECOMPOSED=1    直接读 processor* 分块结果（没 reconstruct 时用）
"""
import os
import re
import sys

from paraview.simple import *  # noqa

try:
    print("ParaView", GetParaViewVersion())
except Exception:
    pass

# ---------------------------------------------------------------- 参数
case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd()
out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(case, "frames")
os.makedirs(out, exist_ok=True)


def read_href(case):
    try:
        txt = open(os.path.join(case, "constant", "hRef")).read()
        m = re.search(r"^\s*value\s+([-+0-9.eE]+)\s*;", txt, re.M)
        return float(m.group(1)) if m else None
    except OSError:
        return None


def env_float(name, default):
    v = os.environ.get(name, "")
    return float(v) if v.strip() else default


VIEW = os.environ.get("VIEW", "flow").lower()
WL = env_float("WL", read_href(case) or 0.0)
BAND = env_float("BAND", 0.004)
ROI = env_float("ROI", 1.0)
SKIP = int(os.environ.get("SKIP", 0))
RES = [int(v) for v in os.environ.get("RES", "1920x1080").lower().split("x")]
DECOMPOSED = os.environ.get("DECOMPOSED", "0") == "1"
INK = [0.06, 0.13, 0.17]

foam = os.path.join(case, "case.foam")
open(foam, "a").close()


def set_if_exists(proxy, names, value):
    """不同 ParaView 版本属性名不一样：能设就设，设不了用默认值。"""
    for nm in names:
        if nm in proxy.ListProperties():
            setattr(proxy, nm, value)
            return True
    return False


# ---------------------------------------------------------------- 读数据
rWater = OpenFOAMReader(FileName=foam)
if DECOMPOSED:
    set_if_exists(rWater, ["CaseType"], "Decomposed Case")
rWater.MeshRegions = ["internalMesh"]
rWater.CellArrays = ["alpha.water", "U", "p_rgh"]
set_if_exists(rWater, ["CreateCellToPoint", "Createcelltopointfiltereddata"], 1)

avail = list(rWater.MeshRegions.Available)
hull_region = next((r for r in avail if r.endswith("/hull")), None)
if hull_region is None:
    print("!! 找不到船体 patch（hull），可用区域：", avail)
    sys.exit(1)

rHull = OpenFOAMReader(FileName=foam)
if DECOMPOSED:
    set_if_exists(rHull, ["CaseType"], "Decomposed Case")
rHull.MeshRegions = [hull_region]
rHull.CellArrays = ["p_rgh", "U"]
set_if_exists(rHull, ["CreateCellToPoint", "Createcelltopointfiltereddata"], 1)

times = list(rWater.TimestepValues)
if not times:
    print("!! 没有时间步 —— 算例还没跑，或者还没 reconstruct")
    sys.exit(1)

water = MergeBlocks(Input=rWater)
hull = MergeBlocks(Input=rHull)

# 船体包围盒 → 尺度
hull.UpdatePipeline(times[0])
b = hull.GetDataInformation().GetBounds()
cx = (b[0] + b[1]) / 2
L = max(b[1] - b[0], 1e-3)
# 半模型的数据只在 y ≤ 0；如果 y 两侧都有数据，就是全模型，不做镜像
FULLMODEL = b[3] > 0.1 * (b[3] - b[2])
B = max(b[3] - b[2], 1e-3) * (1 if FULLMODEL else 2)
print("全模型（不镜像）" if FULLMODEL else "半模型（沿 y=0 镜像显示）")
zmin, zmax = b[4], b[5]
D = max(WL - zmin, 1e-3)                 # 水下深度
YS = env_float("YSLICE", -0.39 * B)
ZS = env_float("ZSLICE", zmin + 0.3 * D)

# |U| 上限：默认取最后一帧的最大值
water.UpdatePipeline(times[-1])
try:
    umax_data = water.PointData["U"].GetRange(-1)[1]
except Exception:
    umax_data = 0.3
UMAX = env_float("UMAX", umax_data)

print(f"算例 {case}\n{len(times)} 帧：{times[0]:g} … {times[-1]:g}   视角 {VIEW}")
print(f"船长 {L:.3f} m  船宽 {B:.3f} m  水线 {WL} m  水下深度 {D:.3f} m")
print(f"竖直切面 y = {YS:.4f}  水平切面 z = {ZS:.4f}  |U| 上限 {UMAX:.3f} m/s")


# ---------------------------------------------------------------- 工具
def mirror_y(src, name):
    """半模型沿 y=0 镜像。兼容 ≤5.13 的 Reflect(Plane) / 6.0 的 ReflectionPlane，
    都不行退回 Transform(y 缩放 -1) + 翻绕向 + 拼接。全模型直接原样返回。"""
    if FULLMODEL:
        return src
    if os.environ.get("MIRROR", "reflect") == "reflect":
        try:
            r = Reflect(Input=src)
            props = r.ListProperties()
            if "ReflectionPlane" in props:
                r.ReflectionPlane = "Y Max"
            elif "Plane" in props:
                r.Plane = "Y Max"
            else:
                raise AttributeError("no plane property")
            set_if_exists(r, ["CopyInput"], 1)
            r.UpdatePipeline()
            return r
        except Exception as e:
            print(f"  ({name}) Reflect 不可用（{e}），改用 Transform 镜像")
    tr = Transform(Input=src)
    tr.Transform = "Transform"
    tr.Transform.Scale = [1.0, -1.0, 1.0]
    rs = ReverseSense(Input=ExtractSurface(Input=tr))
    rs.ReverseCells = 1
    rs.ReverseNormals = 0
    return AppendGeometry(Input=[ExtractSurface(Input=src), rs])


def clip_roi(src, ymirror=True):
    """裁到船体周围（尾流方向多留一些）。"""
    yw = (0.9 + 0.4 * ROI) * L
    lo = [cx - (1.8 + ROI) * L, -yw, zmin - 0.35 * L]
    hi = [cx + (0.8 + 0.4 * ROI) * L, yw if ymirror else 0.0, WL + 0.3 * L]
    c = Clip(Input=src)
    c.ClipType = "Box"
    c.ClipType.Position = lo
    c.ClipType.Length = [hi[i] - lo[i] for i in range(3)]
    set_if_exists(c, ["Invert"], 1)
    set_if_exists(c, ["Crinkleclip", "Crinkle clip"], 0)
    return c


def water_only(src):
    """只保留 alpha.water ≥ 0.5 的部分（切面上把空气去掉）。"""
    t = Threshold(Input=src)
    t.Scalars = ["POINTS", "alpha.water"]
    if set_if_exists(t, ["LowerThreshold"], 0.5):
        set_if_exists(t, ["UpperThreshold"], 1.5)
        set_if_exists(t, ["ThresholdMethod"], "Between")
    else:
        t.ThresholdRange = [0.5, 1.5]
    return t


def free_surface():
    s = Contour(Input=water, ContourBy=["POINTS", "alpha.water"], Isosurfaces=[0.5])
    s.ComputeScalars = 0
    full = mirror_y(s, "water")
    e = Calculator(Input=full, ResultArrayName="elevation", Function="coordsZ")
    return clip_roi(e)


def vslice():
    s = Slice(Input=water)
    s.SliceType = "Plane"
    s.SliceType.Origin = [cx, YS, WL]
    s.SliceType.Normal = [0, 1, 0]
    return clip_roi(water_only(s), ymirror=False)


def hslice(mirror=True):
    s = Slice(Input=water)
    s.SliceType = "Plane"
    s.SliceType.Origin = [cx, 0, ZS]
    s.SliceType.Normal = [0, 0, 1]
    s2 = mirror_y(s, "hslice") if mirror else s
    return clip_roi(s2, ymirror=mirror)


def streamlines():
    st = StreamTracer(Input=water, SeedType="Line")
    st.Vectors = ["POINTS", "U"]
    st.MaximumStreamlineLength = 6 * L
    st.SeedType.Point1 = [cx + 0.75 * L, -1.0 * B, ZS]
    st.SeedType.Point2 = [cx + 0.75 * L, YS - 0.03 * B, ZS]
    st.SeedType.Resolution = 22
    tube = Tube(Input=st)
    tube.Radius = 0.006 * L
    return tube


def style_bar(bar, title, pos):
    bar.Title = title
    bar.ComponentTitle = ""
    bar.TitleFontSize = 14
    bar.LabelFontSize = 12
    bar.ScalarBarLength = 0.3
    bar.TitleColor = INK
    bar.LabelColor = INK
    bar.RangeLabelFormat = "%-#6.3g"
    try:
        bar.WindowLocation = pos
    except Exception:
        pass


def show_umag(src, view, title="|U| [m/s]", bar_pos="Lower Right Corner"):
    d = Show(src, view)
    ColorBy(d, ("POINTS", "U", "Magnitude"))
    lut = GetColorTransferFunction("U")
    lut.ApplyPreset("Viridis (matplotlib)", True)
    lut.VectorMode = "Magnitude"
    lut.RescaleTransferFunction(0.0, UMAX)
    d.SetScalarBarVisibility(view, True)
    style_bar(GetScalarBar(lut, view), title, bar_pos)
    return d


def show_hull(view, color_by_p=False):
    hf = mirror_y(hull, "hull")
    d = Show(hf, view)
    if color_by_p:
        ColorBy(d, ("POINTS", "p_rgh"))
        lut = GetColorTransferFunction("p_rgh")
        lut.ApplyPreset("Cool to Warm", True)
        hf.UpdatePipeline(times[-1])
        try:
            r = hf.PointData["p_rgh"].GetRange()
            m = max(abs(r[0]), abs(r[1]))
            lut.RescaleTransferFunction(-m, m)
        except Exception:
            pass
        d.SetScalarBarVisibility(view, True)
        style_bar(GetScalarBar(lut, view), "p_rgh [m2/s2]", "Lower Right Corner")
    else:
        try:
            ColorBy(d, None)
        except Exception:
            d.ColorArrayName = ["POINTS", ""]
        d.DiffuseColor = [0.35, 0.38, 0.40]
        d.AmbientColor = [0.35, 0.38, 0.40]
    d.Specular = 0.4
    return d


def show_surface(view, opacity=1.0, with_bar=True):
    fs = free_surface()
    d = Show(fs, view)
    ColorBy(d, ("POINTS", "elevation"))
    lut = GetColorTransferFunction("elevation")
    lut.ApplyPreset("Cool to Warm", True)
    lut.RescaleTransferFunction(WL - BAND, WL + BAND)
    d.Opacity = opacity
    d.Specular = 0.3
    if with_bar:
        d.SetScalarBarVisibility(view, True)
        style_bar(GetScalarBar(lut, view), "wave elevation z [m]", "Lower Right Corner")
    return d


# ---------------------------------------------------------------- 视图
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = RES
try:
    view.UseColorPaletteForBackground = 0
except Exception:
    pass
view.Background = [0.93, 0.945, 0.945]
view.OrientationAxesVisibility = 0
CAM = {"parallel": 0, "angle": 28}


def apply_camera():
    view.CameraParallelProjection = CAM["parallel"]
    view.CameraFocalPoint = CAM["focal"]
    view.CameraPosition = CAM["pos"]
    view.CameraViewUp = CAM["up"]
    if CAM["parallel"]:
        view.CameraParallelScale = CAM["pscale"]
    else:
        view.CameraViewAngle = CAM["angle"]

if VIEW == "surface":
    show_surface(view, 1.0, True)
    show_hull(view)
    focal = [cx - 0.5 * L, 0.0, WL]
    CAM['focal'] = focal
    CAM['pos'] = [focal[0] + 2.2 * L, focal[1] - 2.7 * L, focal[2] + 1.6 * L]
    CAM['up'] = [0, 0, 1]
    CAM['angle'] = 28
    caption = "free surface, colored by wave elevation"

elif VIEW == "slice":
    show_umag(vslice(), view)
    show_hull(view)
    CAM['parallel'] = 1
    focal = [cx - 0.3 * L, YS, WL - 0.16 * L]
    CAM['focal'] = focal
    CAM['pos'] = [focal[0], focal[1] - 5 * L, focal[2]]
    CAM['up'] = [0, 0, 1]
    CAM['pscale'] = 0.45 * L
    caption = f"vertical slice y = {YS:.3f} m (through the paddles), |U|"

elif VIEW == "topuw":
    show_umag(hslice(True), view)
    show_hull(view)
    CAM['parallel'] = 1
    focal = [cx - 0.3 * L, 0.0, ZS]
    CAM['focal'] = focal
    CAM['pos'] = [focal[0], focal[1], focal[2] + 5 * L]
    CAM['up'] = [0, 1, 0]
    CAM['pscale'] = 0.7 * L
    caption = f"horizontal slice z = {ZS:.3f} m (paddle depth), |U|"

elif VIEW == "pressure":
    show_hull(view, color_by_p=True)
    show_surface(view, 0.25, False)
    focal = [cx - 0.05 * L, -0.1 * B, WL - 0.5 * D]
    CAM['focal'] = focal
    CAM['pos'] = [focal[0] + 1.8 * L, focal[1] - 1.6 * L, focal[2] - 0.8 * L]
    CAM['up'] = [0, 0, 1]
    CAM['angle'] = 28
    caption = "hull surface pressure p_rgh (dynamic part), seen from below"

else:  # flow
    d_sl = show_umag(vslice(), view)
    d_sl.Opacity = 0.85
    show_umag(streamlines(), view)
    show_hull(view)
    show_surface(view, 0.30, False)
    focal = [cx - 0.3 * L, -0.2 * B, WL - 0.3 * D]
    CAM['focal'] = focal
    CAM['pos'] = [focal[0] + 2.4 * L, focal[1] - 3.1 * L, focal[2] + 0.75 * L]
    CAM['up'] = [0, 0, 1]
    CAM['angle'] = 28
    caption = "underwater: |U| on the paddle-plane slice + streamlines, free surface translucent"

label = Text(Text="")
dLabel = Show(label, view)
dLabel.FontSize = 18
dLabel.Color = INK
try:
    dLabel.WindowLocation = "Upper Left Corner"
except Exception:
    pass

# ---------------------------------------------------------------- 逐帧输出
n = 0
for i, t in enumerate(times):
    if i < SKIP:
        continue
    view.ViewTime = t
    label.Text = f"LTS step {t:g}   |   flow along -x   |   {caption}"
    apply_camera()
    Render(view)
    fn = os.path.join(out, f"frame.{n:04d}.png")
    SaveScreenshot(fn, view, ImageResolution=RES)
    n += 1
    print(f"  帧 {n:3d}/{len(times) - SKIP}  t = {t:g}")

print(f"\n完成：{n} 帧存在 {out}")
