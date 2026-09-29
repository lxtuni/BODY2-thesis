# pvbatch render_hulls.py <out_dir> : 三个船体的两个视角（斜下方看船底、斜上方看甲板）
import sys, os
from paraview.simple import *
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
files = [("M1_current", "hull_M1_current.stl"), ("M2_roundbottom", "hull_M2_roundbottom.stl"), ("M3_streamlined", "hull_M3_streamlined.stl")]
view = GetActiveViewOrCreate("RenderView"); view.ViewSize = [900, 600]; view.Background = [0.97, 0.97, 0.96]
try: view.UseColorPaletteForBackground = 0
except Exception: pass
view.OrientationAxesVisibility = 1
for name, fn in files:
    src = STLReader(FileNames=[os.path.join(os.path.dirname(os.path.abspath(__file__)), fn)])
    d = Show(src, view); d.DiffuseColor = [0.80, 0.80, 0.78]; d.AmbientColor = d.DiffuseColor
    try: ColorBy(d, None)
    except Exception: pass
    # 水面：z = -0.010 的半透明平面
    pl = Plane(); pl.Origin = [-0.20, -0.12, -0.010]; pl.Point1 = [0.20, -0.12, -0.010]; pl.Point2 = [-0.20, 0.12, -0.010]
    dp = Show(pl, view); dp.DiffuseColor = [0.3, 0.55, 0.9]; dp.Opacity = 0.25
    Render(view)                                   # 第一次 Render 会重置相机，先渲染一次
    for tag, pos, foc in (("below", [0.35, -0.30, -0.25], [0.0, 0.0, -0.02]), ("above", [0.35, -0.30, 0.22], [0.0, 0.0, -0.01])):
        view.CameraPosition = pos; view.CameraFocalPoint = foc; view.CameraViewUp = [0, 0, 1]; view.CameraViewAngle = 30
        Render(view); SaveScreenshot(os.path.join(out, f"{name}_{tag}.png"), view, ImageResolution=[900, 600])
    Hide(src, view); Hide(pl, view); Delete(pl); Delete(src)
print("done")
