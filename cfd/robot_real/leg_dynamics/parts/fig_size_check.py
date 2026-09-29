#!/usr/bin/env python3
"""尺寸核算图：船体横截面（后腿尾鳍所在 x ≈ −21）+ 尾鳍在冲程顶/底端的包络。左：A1 展 86 + 对称摆动（顶到船底）；右：A2 展 60 + 下移摆动。"""
import json, math, numpy as np, cadquery as cq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.gp import gp_Pln, gp_Pnt, gp_Dir
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GCPnts import GCPnts_UniformAbscissa
hull = cq.importers.importStep("hull_only.step").solids().vals()[0]
def section(x):
    sec = BRepAlgoAPI_Section(hull.wrapped, gp_Pln(gp_Pnt(x, 0, 0), gp_Dir(1, 0, 0))); sec.Build()
    pts = []; ex = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while ex.More():
        c = BRepAdaptor_Curve(cq.Edge(ex.Current()).wrapped); ga = GCPnts_UniformAbscissa(c, 0.5)
        if ga.IsDone():
            for i in range(1, ga.NbPoints() + 1): p = c.Value(ga.Parameter(i)); pts.append((p.Z(), p.Y()))
        ex.Next()
    return np.array(pts)
pts = section(-21.0)
A_Y = 4.0; Z_LEG = 70.0; WL = 11.0
def fluke_box(b, c, th_top_deg, th_bot_deg, R=101.0, ZC=-52.0):
    """返回冲程顶端与底端尾鳍（展 b、弦 c，与杆同向、绕枢轴）在 y–z 截面里的矩形（z 展向, y 竖直）"""
    out = []
    for th in (th_top_deg, th_bot_deg):
        yJ = A_Y + ZC + R * math.sin(math.radians(th))
        t = 0.21 * c / 2   # 半厚
        out.append((Z_LEG - b / 2, yJ - t - 2.0, b, 2 * t + 4.0))   # +2 mm 摆幅/相位余量
    return out
fig, axs = plt.subplots(1, 2, figsize=(12, 6.2), sharey=True)
cases = [("A1：展 86 mm，对称摆动 ±18°（杆水平上下）", 86.0, 28.7, 18.0, -18.0, "#c0392b"),
         ("A2：展 60 mm，摆动 −18° ± 18°（顶端杆水平）", 60.0, 25.6, 0.0, -36.0, "#1f77b4")]
for ax, (title, b, c, tt, tb, col) in zip(axs, cases):
    ax.scatter(pts[:, 0], pts[:, 1], s=1.5, c="k", label="船体截面 x = −21（后腿尾鳍处）")
    ax.axhline(WL, color="#2a9d8f", ls="--", lw=1); ax.text(-135, WL + 1.5, "压载后水线 +11", color="#2a9d8f", fontsize=8)
    ax.axhline(-35, color="gray", ls=":", lw=1); ax.text(-135, -34, "船底 −35", color="gray", fontsize=8)
    for sgn in (1, -1):
        for i, (z0, y0, w, h) in enumerate(fluke_box(b, c, tt, tb)):
            z0 = sgn * z0 if sgn > 0 else -(z0 + w)
            ax.add_patch(plt.Rectangle((z0, y0), w, h, fc=col, alpha=0.35 if i == 0 else 0.15, ec=col, lw=1))
        ax.plot([sgn * Z_LEG] * 2, [A_Y, A_Y - 52 - 101 * math.sin(math.radians(abs(tt)))], color="0.5", lw=1)   # 腿平面示意
        ax.plot(sgn * Z_LEG, A_Y, "o", color="0.3", ms=4)
    yt = A_Y - 52 + 101 * math.sin(math.radians(tt))
    ax.annotate(f"冲程顶端 y ≈ {yt:.0f}", (Z_LEG - b / 2 - 2, yt + 8), ha="right", fontsize=8, color=col)
    ax.annotate(f"底端 y ≈ {A_Y - 52 + 101 * math.sin(math.radians(tb)):.0f}", (Z_LEG + b / 2 + 2, A_Y - 52 + 101 * math.sin(math.radians(tb))), fontsize=8, color=col)
    ax.set_title(title, fontsize=10); ax.set_aspect("equal"); ax.set_xlim(-140, 140); ax.set_ylim(-125, 30); ax.grid(alpha=0.25)
    ax.set_xlabel("z 横向 / mm（机器人总宽 ≈ %d）" % round(2 * (Z_LEG + b / 2)))
axs[0].set_ylabel("y 竖直 / mm（髋轴 y = 4）")
axs[0].text(0, -60, "尾鳍内侧 18 mm 在冲程顶端\n进入船底 → 干涉", ha="center", color="#c0392b", fontsize=9)
axs[1].text(0, -82, "顶端离船底 ≥ 10 mm\n内侧尖离船侧 −— 全程无干涉", ha="center", color="#1f77b4", fontsize=9)
plt.tight_layout(); plt.savefig("fig_size_check.png", dpi=150); print("ok")
