import trimesh, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib import font_manager
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
def rot(el, az):
    a, e = np.radians(az), np.radians(el)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(e), -np.sin(e)], [0, np.sin(e), np.cos(e)]])
    return Rx @ Rz
def draw(ax, meshes, el, az, zoom=None, edge=False):
    R = rot(el, az); light = np.array([0.3, -0.5, 0.8]); light /= np.linalg.norm(light)
    allp, allc, alld = [], [], []
    for m, col in meshes:
        v = m.vertices[:, [1, 2, 0]].copy(); v[:, 2] *= -1          # 显示坐标：X=y(后) Y=z(横) Z=−x(上)
        p = v @ R.T; tri = p[m.faces]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); n /= (np.linalg.norm(n, axis=1)[:, None] + 1e-12)
        shade = 0.55 + 0.45 * np.clip(n @ (R @ light), 0, 1)
        base = np.array(matplotlib.colors.to_rgb(col))
        allp.append(tri[:, :, :2]); allc.append(shade[:, None] * base[None, :]); alld.append(tri[:, :, 2].mean(axis=1))
    P = np.concatenate(allp); C = np.concatenate(allc); D = np.concatenate(alld)
    o = np.argsort(D)
    ax.add_collection(PolyCollection(P[o], facecolors=C[o], edgecolors="k" if edge else "none", linewidths=0.1))
    ax.set_aspect("equal"); ax.autoscale(); ax.set_xticks([]); ax.set_yticks([])
    if zoom: ax.set_xlim(*zoom[0]); ax.set_ylim(*zoom[1])
if __name__ == "__main__":
    import sys
    ver = sys.argv[1] if len(sys.argv) > 1 else "A"
    meshes = [(trimesh.load(f"strut_{ver}.stl"), "#8a8f96"), (trimesh.load(f"fluke_{ver}.stl"), "#c2185b"), (trimesh.load("flexure_t1.0.stl"), "#d4a017")]
    for m, _ in meshes: print(m.metadata.get("file_name"), "watertight", m.is_watertight, "vol", round(m.volume / 1e3, 2), "tris", len(m.faces))
    fig, axs = plt.subplots(1, 4, figsize=(20, 8), facecolor="#fcfcfb")
    for ax, (ttl, el, az) in zip(axs, [("等轴视图", -30, -35), ("侧视（横向看）", -90, 0), ("后视（沿来流看）", -90, 90), ("俯视（从上看）", 0, 0)]):
        draw(ax, meshes, el, az); ax.set_title(ttl, fontsize=11)
    fig.suptitle(f"月牙尾鳍脚 三件套（版本 {ver}）：灰 = 安装片+支杆+铰座（PETG），红 = 尾鳍（PETG），金 = TPU 柔性铰片", fontsize=12)
    fig.savefig(f"fluke_foot_views_{ver}.png", dpi=110, facecolor="#fcfcfb", bbox_inches="tight")
    # 铰座局部
    fig2, axs2 = plt.subplots(1, 3, figsize=(18, 6), facecolor="#fcfcfb")
    sub = []
    for m, col in meshes:
        keep = (m.vertices[m.faces].mean(axis=1)[:, 0] > 85)
        sub.append((trimesh.Trimesh(m.vertices, m.faces[keep]), col))
    for ax, (ttl, el, az, ed) in zip(axs2, [("铰座侧视（限位爪、柔性片、根部夹块）", -90, 0, True), ("铰座等轴", -35, -40, False), ("铰座后视", -90, 90, True)]):
        draw(ax, sub, el, az, edge=ed); ax.set_title(ttl, fontsize=11)
    fig2.savefig(f"fluke_foot_hinge_detail_{ver}.png", dpi=110, facecolor="#fcfcfb", bbox_inches="tight")
