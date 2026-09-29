import trimesh, numpy as np, matplotlib.pyplot as plt
from render_parts import draw
def prep(fn):
    m = trimesh.load(fn); v = m.vertices.copy(); m2 = trimesh.Trimesh(np.c_[-v[:, 1], v[:, 0], v[:, 2]], m.faces); return m2
orig = prep("fan_foot_orig.stl")
A1, A2 = prep("fluke_straight_A1.stl"), prep("fluke_straight_A2.stl")
fig, axs = plt.subplots(2, 3, figsize=(18, 9), facecolor="#fcfcfb")
views = [("俯视（尾鳍平面）", 0, 0), ("侧视", -90, 0), ("等轴", -30, -35)]
for row, (name, m, col) in enumerate([("A1（展 86，18.5 cm²，112 mm 长，10.3 g）", A1, "#c0392b"), ("A2（展 60，11.5 cm²，98 mm 长，6.9 g）", A2, "#1f77b4")]):
    for ax, (ttl, el, az) in zip(axs[row], views):
        draw(ax, [(orig, "#bbbbbb"), (m, col)], el, az)
        ax.set_title(f"{name} — {ttl}", fontsize=10)
        if el == 0: ax.set_xlim(-20, 105); ax.set_ylim(-50, 50)
        elif az == 0: ax.set_xlim(-20, 105); ax.set_ylim(-14, 14)
fig.suptitle("灰 = 原扇形足（同一位置叠放，30 宽 × 92.5 长）；同比例", fontsize=11)
fig.savefig("fluke_A1_vs_A2.png", dpi=110, facecolor="#fcfcfb", bbox_inches="tight"); print("ok")
