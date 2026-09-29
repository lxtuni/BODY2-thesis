#!/usr/bin/env python3
"""deep_dive_optimal.py / make_optimal_anim.py 共用：读取 open/closed 算例 + 运动学、两态合成、准定常参考、绘图样式"""
import os, sys, json, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); LD = os.path.abspath(os.path.join(HERE, "..", ".."))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]:
    try: font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
INK, INK2, MUTED, SURF, GRID = "#0b0b0b", "#52514e", "#8a8f96", "#fcfcfb", "#e6e5e1"
COL = {"current": "#8a8f96", "V1": "#eb6834", "V2": "#7a5cc0", "V3": "#2a78d6", "V3c12": "#16a085", "V3c06": "#7fb800", "V3trap": "#c2185b", "opt": "#c2185b"}
NAME = {"current": "现步态", "V3": "V3  余弦 · 收拢 30 mm", "V3c12": "V3c12  余弦 · 收拢 12 mm", "V3c06": "V3c06  余弦 · 羽化 5.6 mm", "V3trap": "V3trap  梯形 · 收拢 12 mm"}
RHO = 998.8
def style(ax):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)

# ------------------------------------------------------------------ 读取 + 运动学
def load_ld(g):
    for k in ["T", "DUTY", "SWEEP", "PHI_EXT", "PHI_RET", "FAN_DEG", "FAN_CLOSE", "CLOSED_W", "SWITCH_FRAC", "VEL_LAW", "RAMP"]:
        if k in g and g[k] is not None: os.environ[k] = str(g[k])
        elif k in os.environ: del os.environ[k]
    if "psi0_deg" in g: os.environ["PSI_START"] = str(g["psi0_deg"])
    os.environ["OUT"] = "/tmp/claude-0/ld_tmp"; os.makedirs("/tmp/claude-0/ld_tmp", exist_ok=True)
    sp = importlib.util.spec_from_file_location("ld_%d" % abs(hash(json.dumps(g, sort_keys=True))), os.path.join(LD, "leg_dynamics.py")); ld = importlib.util.module_from_spec(sp); sp.loader.exec_module(ld); return ld

def load(R, N=2400):
    def rd(tag):
        d = np.loadtxt(os.path.join(R, f"force_{tag}.dat"), comments="#"); d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
        m = np.loadtxt(os.path.join(R, f"moment_{tag}.dat"), comments="#"); m = m[np.argsort(m[:, 0])]; _, i = np.unique(m[:, 0], return_index=True); m = m[i]
        return dict(t=d[:, 0], F=d[:, 1:4], tm=m[:, 0], M=m[:, 1:4], gait=json.load(open(os.path.join(R, f"gait_{tag}.json"))))
    c = {k: rd(k) for k in ("open", "closed") if os.path.exists(os.path.join(R, f"force_{k}.dat"))}
    g = (c.get("open") or c["closed"])["gait"]; T, DUTY = g["T"], g["DUTY"]
    ncyc = int(np.floor(min(x["t"][-1] for x in c.values()) / T + 1e-6)); t0 = (ncyc - 1) * T
    tt = np.linspace(t0, t0 + T, N, endpoint=False); tau = (tt - t0) / T
    F = {k: np.stack([np.interp(tt, c[k]["t"], c[k]["F"][:, i]) for i in range(3)], 1) for k in c}
    M = {k: np.stack([np.interp(tt, c[k]["tm"], c[k]["M"][:, i]) for i in range(3)], 1) for k in c}
    ld = load_ld(g)
    psi, phi, _ = ld.gait(tt); PP = [ld.pose(p, q) for p, q in zip(psi, phi)]
    E = np.array([p["E"] for p in PP]); u = np.array([p["u"] for p in PP]); n = np.array([p["n"] for p in PP])
    Ph = E + (ld.STEM + ld.HEAD / 2) * u; tip = E + (ld.STEM + ld.HEAD) * u
    vE = np.gradient(E, tt, axis=0); vP = np.gradient(Ph, tt, axis=0); om = -np.gradient(np.unwrap(psi), tt); rE = E - np.array(ld.O2)
    return dict(R=R, T=T, DUTY=DUTY, ncyc=ncyc, tt=tt, tau=tau, F=F, M=M, gait=g, ld=ld, psi=psi, phi=phi, E=E, u=u, n=n, Ph=Ph, tip=tip, vE=vE, vP=vP, om=om, rE=rE, PP=PP)

def composite(d, to, tc):
    """open 窗口 [to, tc)（按周期取模），其余 closed。返回 dict(F [N,3], M [N,3], ME_y, P [W], w)"""
    tau = d["tau"]; w = ((tau - to) % 1.0) < ((tc - to) % 1.0)
    Fc = np.where(w[:, None], d["F"]["open"], d["F"]["closed"]); Mc = np.where(w[:, None], d["M"]["open"], d["M"]["closed"])
    ME = Mc[:, 1] - (d["rE"][:, 1] * Fc[:, 0] - d["rE"][:, 0] * Fc[:, 2]); P = -(Fc[:, 0] * d["vE"][:, 0] + Fc[:, 2] * d["vE"][:, 1] + ME * d["om"])
    return dict(F=Fc, M=Mc, ME=ME, P=P, w=w)

def smooth(y, tau, wtau=0.02):
    k = max(1, int(round(wtau * len(tau)))); ker = np.ones(k) / k; yp = np.concatenate([y[-k:], y, y[:k]]); return np.convolve(yp, ker, "same")[k:-k]

def qs_force(d, fan, w_closed=None):
    """准定常参考（不含附加质量）：桨头 F = −½ρ C_N A |v_n| v_n n，v_n = 桨头中心速度的法向分量；细杆按 C_D 1.0、杆中点速度。
       open：扇 120°、R 42 mm（18.5 cm²，C_N 1.3）；closed：w × 42 mm（C_N 1.2）。返回与 d['tau'] 对齐的 F [N,2]（x,z）"""
    ld = d["ld"]; n = d["n"]; u = d["u"]; E = d["E"]; tt = d["tt"]
    if fan == "open": A, CN = 0.5 * ld.FAN_R ** 2 * np.deg2rad(ld.FAN_DEG), 1.3
    else: A, CN = (ld.CLOSED_W if w_closed is None else w_closed) * ld.HEAD, 1.2
    vn = np.einsum("ij,ij->i", d["vP"], n); Fh = (-0.5 * RHO * CN * A * np.abs(vn) * vn)[:, None] * n
    Ps = E + ld.STEM / 2 * u; vs = np.gradient(Ps, tt, axis=0); vsn = np.einsum("ij,ij->i", vs, n); Fs = (-0.5 * RHO * 1.0 * ld.STEM * ld.STEM_W * np.abs(vsn) * vsn)[:, None] * n
    return Fh + Fs

