import cadquery as cq, json, numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.gp import gp_Pln, gp_Pnt, gp_Dir
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GCPnts import GCPnts_UniformAbscissa
hull = cq.importers.importStep("hull_only.step").solids().vals()[0]
out = {}
for x in (-60, -40, -21, 0, 20, 60):
    sec = BRepAlgoAPI_Section(hull.wrapped, gp_Pln(gp_Pnt(x, 0, 0), gp_Dir(1, 0, 0)))
    sec.Build(); sh = sec.Shape()
    pts = []
    ex = TopExp_Explorer(sh, TopAbs_EDGE)
    while ex.More():
        c = BRepAdaptor_Curve(cq.Edge(ex.Current()).wrapped)
        ga = GCPnts_UniformAbscissa(c, 1.0)
        if ga.IsDone():
            for i in range(1, ga.NbPoints() + 1):
                p = c.Value(ga.Parameter(i)); pts.append((p.Y(), p.Z()))
        ex.Next()
    pts = np.array(pts)
    # 每个 |z| 带的最低 y = 船底外轮廓
    prof = []
    for z in np.arange(0, 58.5, 2.5):
        m = np.abs(np.abs(pts[:, 1]) - z) < 1.3
        prof.append((float(z), float(pts[m, 0].min()) if m.any() else None, float(pts[m, 0].max()) if m.any() else None))
    out[str(x)] = prof
    print(x, "n", len(pts), " bottom:", [(z, None if lo is None else round(lo, 1)) for z, lo, hi in prof])
json.dump(out, open("hull_profile.json", "w"), indent=1)
