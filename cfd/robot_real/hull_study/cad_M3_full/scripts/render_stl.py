#!/usr/bin/env python3
"""简易软件渲染（画家算法 + 法向光照）：cad_M3_full/BODY2_M3_hull_print.stl → fig_M3_render.png（水下侧视 / 俯视 / 尖端正视）"""
import os, numpy as np, trimesh
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib import font_manager
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "cad_M3_full") if os.path.isdir(os.path.join(HERE, "cad_M3_full")) else (os.path.dirname(HERE) if os.path.basename(HERE) == "scripts" else HERE)
m = trimesh.load(os.path.join(OUT, "BODY2_M3_hull_print.stl"), force="mesh")
V, F = m.vertices, m.faces
Nf = m.face_normals; Ns = m.vertex_normals[m.faces].mean(1); Ns /= np.linalg.norm(Ns, axis=1, keepdims=True)
smooth = (Ns * Nf).sum(1) > np.cos(np.radians(15)); N = np.where(smooth[:, None], Ns, Nf)      # 光滑区用顶点法向平均着色，棱边处用面法向
def R(axis, deg):
    t = np.radians(deg); c, s = np.cos(t), np.sin(t)
    return {"x": np.array([[1, 0, 0], [0, c, -s], [0, s, c]]), "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[axis]
def render(ax, Rm, title, light=(0.3, 0.8, 0.5), base=(0.30, 0.55, 0.85)):
    P = V @ Rm.T; Nn = N @ Rm.T; Nfn = Nf @ Rm.T                       # 视图坐标：屏幕 x = P[:,0]，屏幕 y = P[:,1]，深度 = P[:,2]（朝观察者为正）
    tri = P[F]; depth = tri[:, :, 2].mean(1); front = Nfn[:, 2] > 0
    order = np.argsort(depth); order = order[front[order]]
    L = np.array(light) / np.linalg.norm(light); shade = np.clip(Nn[order] @ L, 0, 1) * 0.75 + 0.25
    cols = np.clip(np.array(base)[None] * shade[:, None], 0, 1)
    pc = PolyCollection(tri[order][:, :, :2], facecolors=cols, edgecolors="none", antialiased=False); ax.add_collection(pc)
    ax.set_xlim(P[:, 0].min() - 5, P[:, 0].max() + 5); ax.set_ylim(P[:, 1].min() - 5, P[:, 1].max() + 5); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=11)
fig, axs = plt.subplots(2, 2, figsize=(16, 11))
# 约定：Rm 左乘 = 在屏幕坐标里继续转；R("x",−90) 后 depth = −y_model（观察者在船底下方），R("x",90) 观察者在上方，R("y",90) 观察者在尖端 −x
render(axs[0, 0], R("x", 25) @ R("y", -35) @ R("x", -90), "从船底斜看（水下新型线 M3：连续曲率，无折角）", light=(0.4, 0.6, 0.7))
render(axs[0, 1], R("x", 25) @ R("y", 35) @ R("x", 90), "从上方斜看（原上部：腿槽 / 法兰 / 孔位保持不变）", light=(0.3, 0.6, 0.75))
render(axs[1, 0], R("x", -12), "侧视（尖端在左；对接面 y = −10 以下为新水下壳）", light=(0.3, 0.5, 0.8))
render(axs[1, 1], R("y", 90), "尖端正视（Spantenriss 方向）", light=(0.5, 0.6, 0.6))
plt.tight_layout(); fn = os.path.join(OUT, "fig_M3_render.png"); plt.savefig(fn, dpi=100); print("→", fn)
