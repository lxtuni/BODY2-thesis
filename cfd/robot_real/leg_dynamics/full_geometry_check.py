#!/usr/bin/env python3
"""
整机算例的几何核对：把真正进入仿真的模型（船体 + 四块桨，含外移与裁扇）画出来，
并算出整个周期里每块桨（含 5.6 mm 厚度）离船体的最近距离，确认不碰船。
    python3 full_geometry_check.py <case目录> [--hull hull_body_only.stl] [--out results_full] [--video]
输出：<out>/geom_<fan>_views.png（三视图 + 间隙曲线）、geom_<fan>_montage.png（周期内 6 个时刻的俯视/侧视）、
      geom_<fan>_clearance.json；--video 再出 geom_<fan>.mp4（需要 ffmpeg）。
"""
import argparse, json, os, sys, struct, subprocess, shutil
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon

def _cjk():
    import glob as g
    c = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]
    c += g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + g.glob("/mnt/c/Windows/Fonts/simhei.ttf")
    for f in c:
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return
        except Exception: pass
_cjk(); plt.rcParams["axes.unicode_minus"] = False
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
LEGCOL = {"FL": "#2a78d6", "FR": "#eb6834", "RL": "#3a9c6e", "RR": "#8e6bd6"}
HULLCOL = "#c9c8c3"

ap = argparse.ArgumentParser()
ap.add_argument("case"); ap.add_argument("--hull", default=None); ap.add_argument("--out", default=None)
ap.add_argument("--video", action="store_true"); ap.add_argument("--frames", type=int, default=100)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
CASE = os.path.abspath(os.path.expanduser(a.case))
g = json.load(open(os.path.join(CASE, "gait.json")))
for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W"]:
    if k in g: os.environ[k] = str(g[k])
if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
import leg_dynamics as ld
T = g["T"]; DUTY = g["DUTY"]; WL = g.get("WL", -0.010); fan = g["fan"]; clip = g.get("fan_clip", {}); legs = g["legs"]
OUT = os.path.abspath(a.out) if a.out else os.path.join(os.environ.get("COPY_TO_FULL", "/mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results_full"))
if not os.path.isdir(os.path.dirname(OUT.rstrip("/"))): OUT = os.path.abspath("results_full")
os.makedirs(OUT, exist_ok=True)
HULL = a.hull or os.path.join(CASE, "constant", "triSurface", "hull.stl")
if not os.path.exists(HULL): HULL = os.path.join(HERE, "..", "hull_body_only.stl")

# ---------------- 船体
def read_stl(fn):
    b = open(fn, "rb").read()
    if b[:5] == b"solid" and b"facet" in b[:400]:
        return np.array([[float(x) for x in ln.split()[1:4]] for ln in b.decode(errors="ignore").splitlines() if ln.strip().startswith("vertex")]).reshape(-1, 3, 3)
    n = struct.unpack("<I", b[80:84])[0]
    return np.frombuffer(b[84:84 + n * 50], dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]))["v"].reshape(-1, 3, 3).astype(float)
TT = read_stl(HULL); HV = TT.reshape(-1, 3)
uu, vv = np.meshgrid(np.linspace(0, 1, 4), np.linspace(0, 1, 4)); uu = uu.ravel(); vv = vv.ravel(); k_ = uu + vv <= 1; uu, vv = uu[k_], vv[k_]
pc = (TT[:, None, 0] + uu[None, :, None] * (TT[:, None, 1] - TT[:, None, 0]) + vv[None, :, None] * (TT[:, None, 2] - TT[:, None, 0])).reshape(-1, 3)
try:
    from scipy.spatial import cKDTree; tree = cKDTree(pc)
    def nearest(S): return tree.query(S)[0]
except Exception:
    P_ = pc[:: max(1, len(pc) // 20000)]
    def nearest(S):
        S = np.asarray(S); d = np.empty(len(S))
        for i in range(0, len(S), 256):
            blk = S[i:i + 256]; d[i:i + 256] = np.sqrt(((blk[:, None, :] - P_[None, :, :]) ** 2).sum(-1)).min(1)
        return d
def hull_outline(axis_pair, zmax=None, zmin=None):
    """船体在某个投影面上的外轮廓（凸包；船体在三个投影面上都近似凸）"""
    V = HV
    if zmax is not None: V = V[V[:, 2] <= zmax]
    if zmin is not None: V = V[V[:, 2] >= zmin]
    P = V[:, axis_pair]
    try:
        from scipy.spatial import ConvexHull
        return P[ConvexHull(P).vertices]
    except Exception:
        c = P.mean(0); ang = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]); r = np.hypot(P[:, 0] - c[0], P[:, 1] - c[1])
        bins = np.linspace(-np.pi, np.pi, 181); idx = np.digitize(ang, bins); out = []
        for b_ in range(1, len(bins)):
            m = idx == b_
            if m.any(): out.append(P[m][r[m].argmax()])
        return np.array(out)
def hull_section_x(x0, half=0.004):
    """船体在 x = x0 附近的横截面轮廓（y–z）"""
    V = HV[np.abs(HV[:, 0] - x0) < half]
    if len(V) < 10: return None
    P = V[:, [1, 2]]; c = P.mean(0); ang = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]); r = np.hypot(P[:, 0] - c[0], P[:, 1] - c[1])
    bins = np.linspace(-np.pi, np.pi, 121); idx = np.digitize(ang, bins); out = []
    for b_ in range(1, len(bins)):
        m = idx == b_
        if m.any(): out.append(P[m][r[m].argmax()])
    return np.array(out)

# ---------------- 桨（与 make_full_cfd.py 完全一致的轮廓 + 裁剪）
def paddle_outline(fan_open):
    S, H, SW, HW = ld.STEM, ld.HEAD, ld.STEM_W, ld.HEAD_W
    if fan_open:
        ang = np.deg2rad(ld.FAN_DEG); th = np.linspace(ang / 2, -ang / 2, 48)
        fanp = [(ld.FAN_R * np.sin(x), S + ld.FAN_R * np.cos(x)) for x in th]
        return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + fanp + [(-SW / 2, S)]
    k = np.linspace(0, 1, 24)
    sh = [(HW / 2 * np.sin(np.pi / 2 * q), S + 0.025 * (1 - np.cos(np.pi / 2 * q))) for q in k]
    right = [(w, z) for w, z in sh if w > SW / 2]
    return [(-SW / 2, 0), (SW / 2, 0), (SW / 2, S)] + right + [(HW / 2, S + H), (-HW / 2, S + H)] + [(-w, z) for w, z in right[::-1]] + [(-SW / 2, S)]
def clipped(fan_open, lat_min, lat_max):
    pts = paddle_outline(fan_open)
    if fan_open:
        S = ld.STEM; out = []
        for (w, z) in pts:
            if w < lat_min or w > lat_max:
                wc = min(max(w, lat_min), lat_max); zc = S + np.sqrt(max(ld.FAN_R ** 2 - wc ** 2, 0.0)); out.append((wc, zc))
            else: out.append((w, z))
        res = []
        for p in out:
            if res and abs(p[0] - res[-1][0]) < 1e-9 and abs(p[1] - res[-1][1]) < 1e-9: continue
            res.append(p)
        ang = np.deg2rad(ld.FAN_DEG); R = ld.FAN_R; hw0 = R * np.sin(ang / 2)
        def side_pt(latc):
            k = abs(latc) / hw0; return (latc, S + k * R * np.cos(ang / 2))
        fixed = []
        for i, p in enumerate(res):
            if lat_max < hw0 - 1e-9 and abs(p[0] - lat_max) < 1e-9 and i > 0 and res[i - 1][1] <= S + 1e-9: fixed.append(side_pt(lat_max))
            fixed.append(p)
            if lat_min > -hw0 + 1e-9 and abs(p[0] - lat_min) < 1e-9 and i + 1 < len(res) and res[i + 1][1] <= S + 1e-9: fixed.append(side_pt(lat_min))
        return np.array(fixed)
    return np.array([(min(max(w, lat_min), lat_max), z) for (w, z) in pts])
hw_full = max(abs(p[0]) for p in paddle_outline(fan == "open"))
OUTL = {}
if "outlines" in g:                                   # 生成器写入的、真正用来造 STL 的轮廓（含 CAD 桨头 / 裁扇）
    for L in legs: OUTL[L] = np.array(g["outlines"][L], dtype=float)
else:
    for L, l in legs.items():
        sgn_in = -1.0 if l["y"] > 0 else 1.0
        hw_in = clip.get(L, {}).get("hw_in", hw_full)
        lat_min, lat_max = (-hw_in, hw_full) if sgn_in < 0 else (-hw_full, hw_in)
        OUTL[L] = clipped(fan == "open", lat_min, lat_max)
NRM_LO = min(g.get("head_nrm", [-ld.THK / 2])[0], g.get("stem_nrm", [-ld.THK / 2])[0])
NRM_HI = max(g.get("head_nrm", [0, ld.THK / 2])[1], g.get("stem_nrm", [0, ld.THK / 2])[1])
THK_TXT = f"{(g['head_nrm'][1]-g['head_nrm'][0])*1000:.1f} mm 桨头板 + {(g['stem_nrm'][1]-g['stem_nrm'][0])*1000:.1f} mm 杆" if "head_nrm" in g else f"{ld.THK*1000:.1f} mm"

nf = a.frames; tt = np.linspace(0, T, nf, endpoint=False)
def leg_pose(L, t):
    l = legs[L]; psi, phi, pw = ld.gait(np.array([t + l["phase"] * T])); P = ld.pose(float(psi[0]), float(phi[0]))
    E = P["E"] + np.array([l["dx"], 0.0]); u = P["u"]; n = np.array([-u[1], u[0]])
    return E, u, n, bool(pw[0])
def plate_points(L, t, dense=True):
    """桨板两个表面上的三维点（用来量距离）和轮廓多边形（用来画）"""
    E, u, n, pw = leg_pose(L, t); o = OUTL[L]; Y = legs[L]["y"]
    poly3 = np.array([(E[0] + s * u[0], Y + lat, E[1] + s * u[1]) for lat, s in o])      # 中面轮廓（x,y,z）
    if not dense: return poly3, pw
    # 面内采样：每个 s 站位上用多边形与直线 s = const 的交点得到真实横向范围
    S = []
    ss = o[:, 1]
    def lat_range(s_):
        xs = []
        for i in range(len(o)):
            (l1, s1), (l2, s2) = o[i], o[(i + 1) % len(o)]
            if (s1 - s_) * (s2 - s_) <= 0 and s1 != s2:
                xs.append(l1 + (l2 - l1) * (s_ - s1) / (s2 - s1))
        return (min(xs), max(xs)) if xs else (0.0, 0.0)
    for s_ in np.linspace(1e-4, ss.max() - 1e-4, 30):
        lo, hi = lat_range(s_)
        for la in np.linspace(lo, hi, 7):
            for off in (NRM_LO, NRM_HI):
                p = E + s_ * u + off * n; S.append((p[0], Y + la, p[1]))
    return np.array(S), poly3, pw

# ---------------- 间隙曲线
dist = {L: np.zeros(nf) for L in legs}; where = {L: None for L in legs}
for i, t in enumerate(tt):
    for L in legs:
        S, _, _ = plate_points(L, t); d = nearest(S); j = int(d.argmin()); dist[L][i] = d[j]
        if where[L] is None or d[j] < where[L][0]: where[L] = (float(d[j]), float(t / T), S[j].tolist())
summary = dict(case=os.path.basename(CASE), fan=fan, hull=os.path.basename(HULL), y_leg={L: legs[L]["y"] for L in legs},
               fan_clip=clip, min_gap_required=g.get("min_gap"), min_clearance_mm={L: round(where[L][0] * 1000, 2) for L in legs},
               at_tau={L: round(where[L][1], 3) for L in legs}, at_point={L: where[L][2] for L in legs},
               all_clear=bool(min(where[L][0] for L in legs) > 0.0), thickness_included=True)
json.dump(summary, open(os.path.join(OUT, f"geom_{fan}_clearance.json"), "w"), indent=1, ensure_ascii=False)
print(f"[{fan}] 整周期最小间隙（含板厚 {THK_TXT}）：" + "  ".join(f"{L} {where[L][0]*1000:.1f} mm (tau {where[L][1]:.2f})" for L in legs) + f"   要求 ≥ {g.get('min_gap', 0.01)*1000:.0f} mm")

# ---------------- 图 1：三视图 + 间隙曲线
top = hull_outline([0, 1], zmax=0.0)          # 水下部分的俯视外轮廓
topwl = hull_outline([0, 1], zmin=WL - 0.003, zmax=WL + 0.003)
side = hull_outline([0, 2])
rear = hull_outline([1, 2])
def draw_views(axs, t, title_extra=""):
    axT, axS, axR = axs
    for ax in axs: ax.set_facecolor(SURF); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.set_aspect("equal"); ax.grid(color=GRID, lw=0.6)
    axT.add_patch(Polygon(top * 1000, closed=True, fc=HULLCOL, ec=INK2, lw=0.8, alpha=0.9)); axT.plot(topwl[:, 0] * 1000, topwl[:, 1] * 1000, color=INK2, lw=0.6, ls="--")
    axS.add_patch(Polygon(side * 1000, closed=True, fc=HULLCOL, ec=INK2, lw=0.8, alpha=0.9)); axS.axhline(WL * 1000, color="#2a78d6", lw=0.8, ls="--")
    axR.add_patch(Polygon(rear * 1000, closed=True, fc=HULLCOL, ec=INK2, lw=0.8, alpha=0.9)); axR.axhline(WL * 1000, color="#2a78d6", lw=0.8, ls="--")
    for L in legs:
        S, poly3, pw = plate_points(L, t); col = LEGCOL[L]; al = 0.85 if pw else 0.45
        axT.add_patch(Polygon(poly3[:, [0, 1]] * 1000, closed=True, fc=col, ec=col, alpha=al, lw=0.8))
        axS.add_patch(Polygon(poly3[:, [0, 2]] * 1000, closed=True, fc=col, ec=col, alpha=al, lw=0.8))
        axR.add_patch(Polygon(poly3[:, [1, 2]] * 1000, closed=True, fc=col, ec=col, alpha=al, lw=0.8))
        d = nearest(S); j = int(d.argmin())
        axT.text(poly3[:, 0].mean() * 1000, poly3[:, 1].mean() * 1000 + (14 if legs[L]["y"] > 0 else -14), f"{L} {d[j]*1000:.0f} mm", ha="center", fontsize=8, color=col)
    axT.set_xlim(-150, 160); axT.set_ylim(-140, 140); axT.set_xlabel("x [mm]（艏 +x）"); axT.set_ylabel("y [mm]"); axT.set_title(f"俯视{title_extra}", loc="left", fontsize=10)
    axS.set_xlim(-150, 160); axS.set_ylim(-120, 40); axS.set_xlabel("x [mm]"); axS.set_ylabel("z [mm]"); axS.set_title("侧视（左舷腿实色 = 划水相）", loc="left", fontsize=10)
    axR.set_xlim(-140, 140); axR.set_ylim(-120, 40); axR.set_xlabel("y [mm]"); axR.set_ylabel("z [mm]"); axR.set_title("艉视", loc="left", fontsize=10)
fig = plt.figure(figsize=(16, 9.5), facecolor=SURF)
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1])
axs = [fig.add_subplot(gs[0, i]) for i in range(3)]
draw_views(axs, 0.0, "  t = 0（FL/RR 划水开始，FR/RL 划水结束）")
ax = fig.add_subplot(gs[1, :]); ax.set_facecolor(SURF); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.grid(axis="y", color=GRID)
for L in legs: ax.plot(tt / T, dist[L] * 1000, color=LEGCOL[L], lw=2, label=f"{L}  最小 {where[L][0]*1000:.1f} mm @ t/T={where[L][1]:.2f}")
req = g.get("min_gap", 0.010) * 1000
ax.axhline(req, color=INK2, lw=1, ls=":"); ax.text(0.005, req + 0.6, f"overset 要求 ≥ {req:.0f} mm", fontsize=9, color=INK2)
ax.axhline(0, color=INK, lw=0.8); ax.set_xlim(0, 1); ax.set_ylim(0, max(45, max(dist[L].max() for L in legs) * 1000 * 1.05))
ax.set_xlabel("t / T（FL、RR 的相位；FR、RL 相差半周期）"); ax.set_ylabel("桨板（含板厚）到船体最近距离 [mm]"); ax.legend(frameon=False, fontsize=9, ncol=4, loc="upper right")
ax.set_title("整个周期内每块桨到船体的最近距离", loc="left", fontsize=10)
clipmsg = "；".join(f"{L} 扇内侧裁到 {clip[L]['hw_in']*1000:.0f} mm" for L in clip) if clip else "无裁剪"
pad_txt = "收拢桨 = STEP 零件 扇形足 的平面轮廓与板厚" if (fan == "closed" and g.get("paddle") == "cad") else ("张开扇 = 120° 扇形模型" if fan == "open" else "理想化轮廓")
fig.suptitle(f"进入仿真的模型（{fan}，{pad_txt}）：前腿 y = ±{abs(legs['FL']['y'])*1000:.0f} mm，后腿 y = ±{abs(legs['RL']['y'])*1000:.0f} mm；{clipmsg}；船体 {os.path.basename(HULL)}，水线 {WL*1000:.0f} mm", fontsize=11)
plt.tight_layout(); fn1 = os.path.join(OUT, f"geom_{fan}_views.png"); plt.savefig(fn1, dpi=110, facecolor=SURF); plt.close(); print("图:", fn1)

# ---------------- 图 2：周期内 6 个时刻的俯视 + 侧视
taus = [0.0, 0.15, 0.3, 0.45, 0.62, 0.85]
fig, axs = plt.subplots(2, 6, figsize=(22, 7.5), facecolor=SURF)
for j, tau in enumerate(taus):
    t = tau * T
    for ax in axs[:, j]: ax.set_facecolor(SURF); [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.set_aspect("equal"); ax.grid(color=GRID, lw=0.5)
    axs[0, j].add_patch(Polygon(top * 1000, closed=True, fc=HULLCOL, ec=INK2, lw=0.6)); axs[1, j].add_patch(Polygon(side * 1000, closed=True, fc=HULLCOL, ec=INK2, lw=0.6)); axs[1, j].axhline(WL * 1000, color="#2a78d6", lw=0.7, ls="--")
    dmin = 1e9
    for L in legs:
        S, poly3, pw = plate_points(L, t); col = LEGCOL[L]; al = 0.85 if pw else 0.45
        axs[0, j].add_patch(Polygon(poly3[:, [0, 1]] * 1000, closed=True, fc=col, ec=col, alpha=al, lw=0.6))
        if legs[L]["y"] > 0: axs[1, j].add_patch(Polygon(poly3[:, [0, 2]] * 1000, closed=True, fc=col, ec=col, alpha=al, lw=0.6))
        dmin = min(dmin, nearest(S).min())
    axs[0, j].set_xlim(-150, 160); axs[0, j].set_ylim(-140, 140); axs[0, j].set_title(f"t/T = {tau:.2f}   最近 {dmin*1000:.1f} mm", fontsize=9.5, loc="left")
    axs[1, j].set_xlim(-150, 160); axs[1, j].set_ylim(-120, 40)
    if j == 0: axs[0, j].set_ylabel("俯视 y [mm]"); axs[1, j].set_ylabel("侧视 z [mm]（只画左舷腿）")
    axs[1, j].set_xlabel("x [mm]")
fig.suptitle(f"周期内 6 个时刻（{fan}；实色 = 划水相，淡色 = 回收相）", fontsize=11)
plt.tight_layout(); fn2 = os.path.join(OUT, f"geom_{fan}_montage.png"); plt.savefig(fn2, dpi=100, facecolor=SURF); plt.close(); print("图:", fn2)

# ---------------- 视频
if a.video:
    fr = os.path.join(OUT, f"_geomframes_{fan}"); shutil.rmtree(fr, ignore_errors=True); os.makedirs(fr)
    for i, t in enumerate(tt):
        fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), facecolor=SURF); draw_views(list(axs), t, f"  t/T = {t/T:.2f}")
        plt.tight_layout(); plt.savefig(os.path.join(fr, f"f{i:04d}.png"), dpi=90, facecolor=SURF); plt.close(fig)
    mp4 = os.path.join(OUT, f"geom_{fan}.mp4")
    r = subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "12", "-i", os.path.join(fr, "f%04d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", mp4])
    if r.returncode == 0: shutil.rmtree(fr); print("视频:", mp4)
    else: print("ffmpeg 失败，帧在", fr)
