#!/usr/bin/env python3
"""
BODY2 船体形状比选：三个模型的几何生成、静水力学与三视图
  M1  现船体（hull_body_only.stl 原样；本脚本只做测量与画图）
  M2  现船体的平底改成圆滑底：甲板轮廓、船长、最大深度不变，横剖面改为从甲板边到龙骨的半椭圆
  M3  低阻方案：甲板轮廓不变（架子/腿的安装位置不变），水下改为流线型圆舭船身——半椭圆剖面 + 抛物型龙骨线（两端上翘），
      中部吃水按"水线相同 (z = −0.010) 且排水量与 M1 相同"确定
用法：python3 hull_shapes.py [--out DIR]   → M2.stl / M3.stl（M1 复制自原 STL）、hull_hydrostatics.json、fig_hull_3views.png
"""
import os, sys, json, struct, argparse
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

def _cjk():
    import glob as g
    c = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]
    c += g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + g.glob("/mnt/c/Windows/Fonts/simhei.ttf")
    for f in c:
        try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return
        except Exception: pass
_cjk(); plt.rcParams["axes.unicode_minus"] = False
INK, INK2, SURF, BLUE, RED, ORANGE, GREEN, GRAY, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6", "#e34948", "#eb6834", "#3a9d5d", "#a8a7a1", "#e6e5e1"

ap = argparse.ArgumentParser()
ap.add_argument("--stl", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "robot_real", "hull_body_only.stl"))
ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
ap.add_argument("--wl", type=float, default=-0.010)
ap.add_argument("--n2", type=float, default=2.0, help="M2 剖面超椭圆指数（2 = 椭圆，越大越接近方底）")
ap.add_argument("--n3", type=float, default=2.0, help="M3 剖面超椭圆指数")
ap.add_argument("--p3", type=float, default=0.8, help="M3 龙骨线形状指数 T(x) ∝ (1−ξ²)^p（1 = 抛物线，0.5 = 椭圆，越小中部越平）")
a = ap.parse_args()
WL = a.wl; RHO = 998.8; OUT = a.out; os.makedirs(OUT, exist_ok=True)
STL = a.stl if os.path.exists(a.stl) else "/home/claude/robot_real/hull_body_only.stl"

# ---------------- 读原 STL，量外轮廓 ----------------
def read_stl(fn):
    data = open(fn, "rb").read(); n = struct.unpack("<I", data[80:84])[0]
    tri = np.frombuffer(data[84:84 + n * 50], dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]))
    return tri["v"].reshape(-1, 3, 3).astype(float)
def write_stl(fn, tris, name="hull"):
    tris = np.asarray(tris, dtype=np.float32); nrm = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-20)
    with open(fn, "wb") as f:
        f.write(name.encode().ljust(80, b"\0")); f.write(struct.pack("<I", len(tris)))
        rec = np.zeros(len(tris), dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")])); rec["n"] = nrm; rec["v"] = tris.reshape(-1, 9); f.write(rec.tobytes())
V = read_stl(STL)
def slice_plane(axis, val):
    d = V[:, :, axis] - val; keep = ~((d > 0).all(1) | (d < 0).all(1)); segs = []
    for t, dd in zip(V[keep], d[keep]):
        pts = []
        for i in range(3):
            j = (i + 1) % 3
            if (dd[i] > 0) != (dd[j] > 0): w = dd[i] / (dd[i] - dd[j]); pts.append(t[i] + w * (t[j] - t[i]))
        if len(pts) == 2: segs.append(pts)
    return np.array(segs)
def outer_hw(segs, z0):
    ys = []
    for p_, q_ in segs:
        da, db = p_[2] - z0, q_[2] - z0
        if (da > 0) != (db > 0): w = da / (da - db); ys.append(abs(p_[1] + w * (q_[1] - p_[1])))
    return max(ys) if ys else np.nan
Z_DECK = float(V[:, :, 2].max())                      # 0.0155
xs_m = np.round(np.arange(-0.125, 0.1351, 0.0025), 5); zs_m = np.round(np.arange(-0.035, Z_DECK + 1e-6, 0.0025), 5)
Wm = np.full((len(xs_m), len(zs_m)), np.nan); keel_m = np.full(len(xs_m), np.nan)
for i, x0 in enumerate(xs_m):
    s = slice_plane(0, x0)
    if len(s) == 0: continue
    keel_m[i] = s.reshape(-1, 3)[:, 2].min()
    for j, z0 in enumerate(zs_m): Wm[i, j] = outer_hw(s, z0)
# 甲板边缘轮廓（去掉架子槽口：槽口区间 z ≥ −2.5 mm 的宽度被架子墙 34 / 27 mm 截断，用 z = −5 mm 的轮廓 + 外飘斜率外推）
jz5 = int(np.argmin(np.abs(zs_m + 0.005))); jdeck = len(zs_m) - 1
notch = ((xs_m > -0.100) & (xs_m < -0.050)) | ((xs_m > 0.050) & (xs_m < 0.090))
w_deck = Wm[:, jdeck].copy()
slope = 0.36                                             # 甲板附近外飘（−5 → +15 mm 高度里半宽增加 ~7 mm）
w_deck[notch] = Wm[notch, jz5] + slope * (Z_DECK + 0.005)
ok = np.isfinite(w_deck)
# 光顺甲板轮廓（对 x 做样条），并让两端收尖
from numpy.polynomial import polynomial as Pn
xd, wd = xs_m[ok], w_deck[ok]
def smooth1(y, k=5):
    yp = np.pad(y, (k, k), mode="edge"); ker = np.ones(2 * k + 1) / (2 * k + 1); return np.convolve(yp, ker, "same")[k:-k]
wd_s = smooth1(wd, 3)
X_STERN, X_BOW = -0.128, 0.138                           # 程序化船体的两端（收尖点，略超出原 STL 的小端面）
def deck_hw(x):
    x = np.asarray(x, dtype=float); w = np.interp(x, xd, wd_s)
    # 端部过渡到 0
    w = np.where(x < xd[0], wd_s[0] * np.clip((x - X_STERN) / (xd[0] - X_STERN), 0, 1), w)
    w = np.where(x > xd[-1], wd_s[-1] * np.clip((X_BOW - x) / (X_BOW - xd[-1]), 0, 1), w)
    return w
keel_ok = np.isfinite(keel_m)
def keel_current(x):
    x = np.asarray(x, dtype=float); k = np.interp(x, xs_m[keel_ok], keel_m[keel_ok])
    k = np.where(x < xs_m[keel_ok][0], Z_DECK - (Z_DECK - keel_m[keel_ok][0]) * np.clip((x - X_STERN) / (xs_m[keel_ok][0] - X_STERN), 0, 1), k)
    k = np.where(x > xs_m[keel_ok][-1], Z_DECK - (Z_DECK - keel_m[keel_ok][-1]) * np.clip((X_BOW - x) / (X_BOW - xs_m[keel_ok][-1]), 0, 1), k)
    return np.minimum(k, Z_DECK - 1e-4)

# ---------------- 程序化船体：水线以上 = 原船体外表面（按 (x,z) 表采样，含架子槽口），水线以下 = 新的圆滑底 ----------------
Z_T = WL                                                  # 新底与原船壳的分界高度（= 水线：水线以上完全不动）
def w1(x, z):
    """原船体外半宽（表格双线性插值；表外为 nan）"""
    i = np.interp(x, xs_m, np.arange(len(xs_m))); j = np.interp(z, zs_m, np.arange(len(zs_m)))
    i0, j0 = int(np.clip(np.floor(i), 0, len(xs_m) - 2)), int(np.clip(np.floor(j), 0, len(zs_m) - 2)); fi, fj = i - i0, j - j0
    q = Wm[i0:i0 + 2, j0:j0 + 2]
    if np.isnan(q).any():
        q = np.where(np.isnan(q), np.nanmax(q) if np.isfinite(q).any() else np.nan, q)
    return (1 - fi) * (1 - fj) * q[0, 0] + fi * (1 - fj) * q[1, 0] + (1 - fi) * fj * q[0, 1] + fi * fj * q[1, 1]
def section_new(x, keel_fn, n_exp, nj=28):
    """一个站位的右舷外轮廓（从龙骨到甲板）：若原船体在此站位有水下部分，则 [新龙骨点] + 椭圆弧（到水线）+ 原船壳采样点（z > 水线）；
       否则完全用原船壳（z 从原龙骨到甲板）。表格范围之外的站位按甲板轮廓收尖。"""
    k = float(keel_fn(x)); k1 = float(keel_current(x)); pts = []
    if k1 < Z_T - 1e-6 and k < Z_T - 1e-6:
        wt = float(w1(x, Z_T)); wt = wt if np.isfinite(wt) else 0.0; D = Z_T - k
        for ph in np.linspace(np.pi / 2, 0, nj):                       # 龙骨 → 水线
            pts.append((wt * np.cos(ph) ** (2 / n_exp), Z_T - D * np.sin(ph) ** (2 / n_exp)))
        z0 = Z_T
    else:
        pts.append((0.0, k1)); z0 = k1
    for z in zs_m:
        if z > z0 + 1e-6 and z < Z_DECK - 1e-6:
            w = w1(x, z); pts.append((float(w) if np.isfinite(w) else pts[-1][0], float(z)))
    wd_ = w1(x, Z_DECK); pts.append((float(wd_) if np.isfinite(wd_) else pts[-1][0], Z_DECK))
    P = np.array(pts)
    # 表格范围外（船体两端的小端面之外）：按甲板轮廓比例收尖
    if x < xd[0] or x > xd[-1]:
        ref = float(w1(xd[0] if x < xd[0] else xd[-1], Z_DECK)); P[:, 0] *= float(deck_hw(x)) / max(ref, 1e-9)
    return P
def zipper(A, B, tris, xa, xb, mirror=False):
    """把两条从龙骨到甲板的折线 A(x=xa)、B(x=xb) 缝成三角形（按归一化弧长贪心）；B 只有一个点时为扇形"""
    def P3(P, x, k): y = P[k, 0] * (-1 if mirror else 1); return np.array([x, y, P[k, 1]])
    def emit(t_): tris.append(t_[::-1] if mirror else t_)
    if len(B) == 1:
        for i in range(len(A) - 1): emit([P3(A, xa, i), P3(B, xb, 0), P3(A, xa, i + 1)])
        return
    if len(A) == 1:
        for j in range(len(B) - 1): emit([P3(A, xa, 0), P3(B, xb, j), P3(B, xb, j + 1)])
        return
    def cum(P):
        d = np.hypot(np.diff(P[:, 0]), np.diff(P[:, 1])); c = np.concatenate([[0], np.cumsum(d)]); return c / max(c[-1], 1e-12)
    sa, sb = cum(A), cum(B); i = j = 0
    while i < len(A) - 1 or j < len(B) - 1:
        if j == len(B) - 1 or (i < len(A) - 1 and sa[i + 1] <= sb[j + 1]):
            emit([P3(A, xa, i), P3(B, xb, j), P3(A, xa, i + 1)]); i += 1
        else:
            emit([P3(A, xa, i), P3(B, xb, j), P3(B, xb, j + 1)]); j += 1
def build_hull(keel_fn, n_exp, nx=150):
    xs_in = np.linspace(X_STERN + 0.0015, X_BOW - 0.0015, nx)
    secs = [np.array([[0.0, Z_DECK]])] + [section_new(x, keel_fn, n_exp) for x in xs_in] + [np.array([[0.0, Z_DECK]])]   # 两端收成甲板高度的一个点
    xs = np.concatenate([[X_STERN], xs_in, [X_BOW]]); tris = []
    for i in range(len(xs) - 1):
        A, B = secs[i], secs[i + 1]
        zipper(A, B, tris, xs[i], xs[i + 1]); zipper(A, B, tris, xs[i], xs[i + 1], mirror=True)
        # 甲板封顶（用中线点）
        c0 = np.array([xs[i], 0, Z_DECK]); c1 = np.array([xs[i + 1], 0, Z_DECK])
        pa, pb = np.array([xs[i], A[-1, 0], Z_DECK]), np.array([xs[i + 1], B[-1, 0], Z_DECK]); qa, qb = pa * [1, -1, 1], pb * [1, -1, 1]
        tris.append([pa, c1, c0]); tris.append([pa, pb, c1]); tris.append([qa, c0, c1]); tris.append([qa, c1, qb])   # 甲板法向朝 +z
    tris = np.array(tris)
    ar = 0.5 * np.linalg.norm(np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1); tris = tris[ar > 1e-14]
    vol = np.einsum("ij,ij->i", tris[:, 0], np.cross(tris[:, 1], tris[:, 2])).sum() / 6
    if vol < 0: tris = tris[:, ::-1]
    return tris, secs[1:-1], xs_in

# ---------------- 静水力学（按剖面积分，水线以下） ----------------
def hydro_from_sections(xs, secs):
    """secs[i] = (y, z) 剖面外轮廓（右舷，从龙骨到甲板）；返回排水体积、湿表面、水线宽/长、最大剖面积、水线面惯量等"""
    A = np.zeros(len(xs)); Sper = np.zeros(len(xs)); bwl = np.zeros(len(xs)); zc = np.zeros(len(xs)); ymax = np.zeros(len(xs))
    for i, (y, z) in enumerate(secs):
        k = z <= WL + 1e-9
        if k.sum() < 2 and not (z.min() < WL): continue
        # 截到水线
        yy, zz = y[k], z[k]
        if (~k).any():
            j = int(np.argmax(~k));
            if j > 0:
                w_ = (WL - z[j - 1]) / (z[j] - z[j - 1]); yy = np.append(yy, y[j - 1] + w_ * (y[j] - y[j - 1])); zz = np.append(zz, WL)
        A[i] = 2 * np.trapezoid(yy, zz) if len(zz) > 1 else 0.0
        Sper[i] = 2 * np.sum(np.hypot(np.diff(yy), np.diff(zz)))
        bwl[i] = 2 * yy[-1] if len(yy) else 0.0; ymax[i] = yy.max() if len(yy) else 0.0
        zc[i] = (np.trapezoid(yy * zz, zz) / max(np.trapezoid(yy, zz), 1e-12)) if len(zz) > 1 else WL
    vol = np.trapezoid(A, xs); S = np.trapezoid(Sper, xs); Awp = np.trapezoid(bwl, xs); Iwp = np.trapezoid((bwl / 2) ** 3 * 2 / 3, xs)
    xb = np.trapezoid(A * xs, xs) / max(vol, 1e-12); KB = WL - np.trapezoid(A * zc, xs) / max(vol, 1e-12)   # 浮心在水线下深度
    wet = xs[A > 1e-9]; LWL = float(wet.max() - wet.min()) if len(wet) else 0.0
    BWL = float(bwl.max()); Amax = float(A.max()); T = None
    return dict(vol_cm3=vol * 1e6, disp_g=vol * RHO * 1e3, S_cm2=S * 1e4, Awp_cm2=Awp * 1e4, LWL_mm=LWL * 1e3, BWL_mm=BWL * 1e3, Amax_cm2=Amax * 1e4,
                BM_mm=Iwp / max(vol, 1e-12) * 1e3, KB_below_WL_mm=KB * 1e3, LCB_mm=xb * 1e3, Cp=vol / max(Amax * LWL, 1e-12), A=A, bwl=bwl)

# M1 剖面（原 STL 外轮廓）
secs1 = []
for i, x0 in enumerate(xs_m):
    zz = zs_m[np.isfinite(Wm[i])]; yy = Wm[i][np.isfinite(Wm[i])]
    if len(zz) == 0 or not np.isfinite(keel_m[i]): secs1.append((np.array([0.0]), np.array([Z_DECK]))); continue
    secs1.append((np.concatenate([[0.0], yy]), np.concatenate([[keel_m[i]], zz])))
H1 = hydro_from_sections(xs_m, secs1); T1 = WL - keel_m[keel_ok].min(); H1["T_mm"] = T1 * 1e3; H1["Cb"] = H1["vol_cm3"] * 1e-6 / (H1["LWL_mm"] * H1["BWL_mm"] * 1e-6 * T1)
H1["Cm"] = H1["Amax_cm2"] * 1e-4 / (H1["BWL_mm"] * 1e-3 * T1)

# M2：龙骨线不变、剖面半椭圆
tris2, S2, xs2 = build_hull(keel_current, a.n2)
def secs_of(S): return [(sec[:, 0], sec[:, 1]) for sec in S]       # 从龙骨到甲板
H2 = hydro_from_sections(xs2, secs_of(S2)); T2 = WL - keel_current(xs2).min(); H2["T_mm"] = T2 * 1e3
H2["Cb"] = H2["vol_cm3"] * 1e-6 / (H2["LWL_mm"] * H2["BWL_mm"] * 1e-6 * T2); H2["Cm"] = H2["Amax_cm2"] * 1e-4 / (H2["BWL_mm"] * 1e-3 * T2)
# M2 若要与 M1 同排水量，水线要下沉到哪里（等效吃水增加）
def vol_at(S, xs, wl):
    global WL
    old = WL; WL = wl; h = hydro_from_sections(xs, secs_of(S)); WL = old; return h["vol_cm3"]
wl_lo, wl_hi = -0.034, -0.005
for _ in range(40):
    mid = 0.5 * (wl_lo + wl_hi)
    if vol_at(S2, xs2, mid) < H1["vol_cm3"]: wl_lo = mid
    else: wl_hi = mid
H2["WL_equal_disp_mm"] = 0.5 * (wl_lo + wl_hi) * 1e3

# M3：抛物型龙骨线（水线长与 M1 相同、两端上翘到与 M1 同高），中部吃水 T_max 由"同水线、同排水量"解出
wetx = xs_m[(H1["A"] > 1e-9)]; XW0, XW1 = float(wetx.min()), float(wetx.max()); XC = 0.5 * (XW0 + XW1); AH = 0.5 * (XW1 - XW0)
def keel3_factory(Tmax, p=a.p3):
    def keel3(x):
        x = np.asarray(x, dtype=float); xi = (x - XC) / AH
        k = np.where(np.abs(xi) < 1, WL - Tmax * np.clip(1 - xi ** 2, 0, None) ** p, WL)
        # 水线以外：船体 = 原船体（两端上翘部分不变）
        k1 = keel_current(x); k = np.where((x < XW0) | (x > XW1), k1, np.minimum(k, k1 + 1e9))
        return np.minimum(k, Z_DECK - 1e-4)
    return keel3
lo, hi = 0.015, 0.09
for _ in range(30):
    Tm = 0.5 * (lo + hi); tris3, S3, xs3 = build_hull(keel3_factory(Tm), a.n3, nx=70); v3 = hydro_from_sections(xs3, secs_of(S3))["vol_cm3"]
    if v3 < H1["vol_cm3"]: lo = Tm
    else: hi = Tm
T3 = 0.5 * (lo + hi); tris3, S3, xs3 = build_hull(keel3_factory(T3), a.n3)
H3 = hydro_from_sections(xs3, secs_of(S3)); H3["T_mm"] = T3 * 1e3; H3["Cb"] = H3["vol_cm3"] * 1e-6 / (H3["LWL_mm"] * H3["BWL_mm"] * 1e-6 * T3); H3["Cm"] = H3["Amax_cm2"] * 1e-4 / (H3["BWL_mm"] * 1e-3 * T3)
H3["keel_min_mm"] = float(keel3_factory(T3)(xs3).min() * 1e3)
# 写 STL
write_stl(os.path.join(OUT, "hull_M2_roundbottom.stl"), tris2, "M2"); write_stl(os.path.join(OUT, "hull_M3_streamlined.stl"), tris3, "M3")
import shutil; shutil.copyfile(STL, os.path.join(OUT, "hull_M1_current.stl"))
for nm, tr in (("M2", tris2), ("M3", tris3)):
    vol = np.einsum("ij,ij->i", tr[:, 0], np.cross(tr[:, 1], tr[:, 2])).sum() / 6
    print(f"{nm}: {len(tr)} 三角形，封闭体积 {vol*1e6:.1f} cm³")
# 原 STL 湿表面（近似）
area1 = 0.5 * np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1); c1 = V.mean(1)
H1["S_cm2_tri"] = float(area1[c1[:, 2] < WL].sum() * 1e4)
def clean(h): return {k: (float(v) if isinstance(v, (float, np.floating, int)) else v) for k, v in h.items() if k not in ("A", "bwl")}
hyd = {"M1_current": clean(H1), "M2_roundbottom": clean(H2), "M3_streamlined": clean(H3), "WL": WL, "params": dict(n2=a.n2, n3=a.n3, p3=a.p3, T3_mm=T3 * 1e3)}
json.dump(hyd, open(os.path.join(OUT, "hull_hydrostatics.json"), "w"), indent=1, ensure_ascii=False)
for k_, h in hyd.items():
    if isinstance(h, dict) and "disp_g" in h: print(f"[{k_:<15}] 排水量 {h['disp_g']:.0f} g  湿表面 {h['S_cm2']:.0f} cm²  水线 L×B {h['LWL_mm']:.0f}×{h['BWL_mm']:.0f} mm  吃水 {h['T_mm']:.1f} mm  最大剖面 {h['Amax_cm2']:.1f} cm²  Cb {h['Cb']:.2f} Cp {h['Cp']:.2f} Cm {h['Cm']:.2f}  BM {h['BM_mm']:.0f} mm  KB(水线下) {h['KB_below_WL_mm']:.1f} mm  LCB {h['LCB_mm']:+.0f} mm" + (f"  同排水量水线 {h['WL_equal_disp_mm']:.1f} mm" if "WL_equal_disp_mm" in h else ""))

# ---------------- 三视图 ----------------
stations = [-0.11, -0.09, -0.07, -0.05, -0.03, -0.01, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11]
def contour_hw(S, xs, z0):
    """程序化船体在高度 z0 的半宽（无交点 → nan）"""
    out = np.full(len(xs), np.nan)
    for i, sec in enumerate(S):
        y, z = sec[:, 0], sec[:, 1]
        if z[0] > z0 or z[-1] < z0: continue
        out[i] = np.interp(z0, z, y)
    return out
def closed_outline(xs, hw):
    m = np.isfinite(hw); x_, w_ = xs[m], hw[m]
    return np.r_[x_, x_[::-1]] * 1e3, np.r_[w_, -w_[::-1]] * 1e3
fig, axs = plt.subplots(3, 3, figsize=(19, 12.5), facecolor=SURF, gridspec_kw=dict(width_ratios=[2.2, 2.2, 1.3]))
titles = ["M1  现船体（平底托盘，原 STL）", f"M2  现船体改圆滑底：水线以上不动，水线以下换成半椭圆底，龙骨线不变（同水线排水量 {H2['disp_g']:.0f} g）",
          f"M3  低阻方案：水线以上不动，水下为流线型圆舭 + 抛物型龙骨（同水线、同排水量 {H3['disp_g']:.0f} g，中部吃水 {H3['T_mm']:.0f} mm）"]
def style(ax):
    ax.set_facecolor(SURF); ax.set_aspect("equal"); ax.grid(color=GRID, lw=0.7)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
xx = np.linspace(X_STERN, X_BOW, 300)
LEGS = ((63.4, 68.8), (63.4, -68.8), (-80.1, 70), (-80.1, -70))
for row, (nm, H) in enumerate((("M1", H1), ("M2", H2), ("M3", H3))):
    S, xsP = (None, None) if nm == "M1" else ((S2, xs2) if nm == "M2" else (S3, xs3))
    # ---- 俯视
    ax = axs[row][0]
    if nm == "M1":
        for z0, col, lw in ((Z_DECK - 0.002, INK, 1.0), (WL, BLUE, 1.4), (-0.030, GRAY, 1.0)):
            for p_, q_ in slice_plane(2, z0): ax.plot([p_[0] * 1e3, q_[0] * 1e3], [p_[1] * 1e3, q_[1] * 1e3], color=col, lw=lw)
        ax.plot([], [], color=INK, label="甲板边（原 STL，含架子槽口与内壁）"); ax.plot([], [], color=BLUE, lw=1.4, label=f"水线 z = {WL*1e3:.0f} mm"); ax.plot([], [], color=GRAY, label="z = −30 mm 轮廓")
    else:
        ax.plot(*closed_outline(xsP, contour_hw(S, xsP, Z_DECK - 0.002)), color=INK, lw=1.0, label="甲板边（同 M1）")
        ax.plot(*closed_outline(xsP, contour_hw(S, xsP, WL)), color=BLUE, lw=1.4, label="水线（同 M1）")
        ax.plot(*closed_outline(xsP, contour_hw(S, xsP, -0.030)), color=GRAY, lw=1.0, label="z = −30 mm 轮廓")
        if nm == "M3": ax.plot(*closed_outline(xsP, contour_hw(S, xsP, -0.040)), color=GRAY, lw=0.8, ls="--", label="z = −40 mm 轮廓")
        for p_, q_ in slice_plane(2, -0.030): ax.plot([p_[0] * 1e3, q_[0] * 1e3], [p_[1] * 1e3, q_[1] * 1e3], color=GRAY, lw=0.6, ls=":", alpha=0.8)
        ax.plot([], [], color=GRAY, lw=0.6, ls=":", label="M1 的 z = −30 mm 轮廓")
    for xl, yl in LEGS: ax.plot(xl, yl, "x", color=RED, ms=8, mew=1.6)
    ax.plot([], [], "x", color=RED, label="腿的髋轴 x / 桨平面 y")
    style(ax); ax.set_xlim(-140, 150); ax.set_ylim(-80, 80); ax.set_xlabel("x [mm]（+x = 艏）"); ax.set_ylabel("y [mm]"); ax.set_title(titles[row] + "   ——  俯视", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    # ---- 侧视
    ax = axs[row][1]
    if nm == "M1":
        for p_, q_ in slice_plane(1, 0.0): ax.plot([p_[0] * 1e3, q_[0] * 1e3], [p_[2] * 1e3, q_[2] * 1e3], color=INK, lw=1.0)
        for p_, q_ in slice_plane(1, 0.040): ax.plot([p_[0] * 1e3, q_[0] * 1e3], [p_[2] * 1e3, q_[2] * 1e3], color=GRAY, lw=0.7)
        ax.plot([], [], color=INK, label="中纵剖面 y = 0（含内壁、盖板）"); ax.plot([], [], color=GRAY, lw=0.7, label="y = 40 mm 纵剖面")
    else:
        kf = keel_current if nm == "M2" else keel3_factory(T3)
        ax.plot(xx * 1e3, kf(xx) * 1e3, color=INK, lw=1.5, label="龙骨线 y = 0"); ax.plot([X_STERN * 1e3, X_BOW * 1e3], [Z_DECK * 1e3] * 2, color=INK, lw=1.0)
        ax.plot(xx * 1e3, keel_current(xx) * 1e3, color=GRAY, lw=0.9, ls="--", label="M1 龙骨线")
        zb = np.array([np.interp(0.040, sec[::-1, 0], sec[::-1, 1]) if (sec[:, 0].max() > 0.040 and sec[0, 1] < WL) else np.nan for sec in S])
        ax.plot(xsP * 1e3, zb * 1e3, color=GRAY, lw=0.8, label="y = 40 mm 纵剖面（船底）")
    ax.axhline(WL * 1e3, color=BLUE, lw=1.2, ls="--"); ax.text(-138, WL * 1e3 + 1.5, "WL", color=BLUE, fontsize=9)
    for xl in (63.4, -80.1): ax.plot(xl, 4.0, "x", color=RED, ms=8, mew=1.6)
    style(ax); ax.set_xlim(-140, 150); ax.set_ylim(-70, 25); ax.set_xlabel("x [mm]"); ax.set_ylabel("z [mm]"); ax.set_title("侧视（纵剖面）", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    # ---- 横剖面
    ax = axs[row][2]; cmap = plt.cm.viridis
    for x0 in stations:
        col = cmap((x0 + 0.12) / 0.24); sgn = 1 if x0 >= 0 else -1
        segs = slice_plane(0, x0)
        if nm == "M1":
            for p_, q_ in segs:
                if (p_[1] * sgn >= -1e-4) and (q_[1] * sgn >= -1e-4): ax.plot([p_[1] * 1e3, q_[1] * 1e3], [p_[2] * 1e3, q_[2] * 1e3], color=col, lw=1.0)
        else:
            i = int(np.argmin(np.abs(xsP - x0))); sec = S[i]
            ax.plot(sgn * sec[:, 0] * 1e3, sec[:, 1] * 1e3, color=col, lw=1.3)
            for p_, q_ in segs:
                if (p_[1] * sgn >= -1e-4) and (q_[1] * sgn >= -1e-4) and min(p_[2], q_[2]) < WL: ax.plot([p_[1] * 1e3, q_[1] * 1e3], [p_[2] * 1e3, q_[2] * 1e3], color=col, lw=0.6, ls="--", alpha=0.6)
    ax.axhline(WL * 1e3, color=BLUE, lw=1.2, ls="--"); ax.axvline(0, color=INK2, lw=0.6)
    ax.text(-66, 19, "艉部站位 x = −110…−10", fontsize=8, color=INK2); ax.text(6, 19, "艏部站位 x = +10…+110", fontsize=8, color=INK2)
    style(ax); ax.set_xlim(-70, 70); ax.set_ylim(-75, 25); ax.set_xlabel("y [mm]"); ax.set_ylabel("z [mm]"); ax.set_title("横剖面（每 20 mm 一站，颜色 = x；虚线 = M1 水下轮廓）", loc="left", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_hull_3views.png"), dpi=110, facecolor=SURF); plt.close()
print("图:", os.path.join(OUT, "fig_hull_3views.png"))
