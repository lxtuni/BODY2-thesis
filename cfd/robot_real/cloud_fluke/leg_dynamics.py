#!/usr/bin/env python3
"""
BODY2 单腿动力模型 —— 第 0 级（运动学）+ 第 1 级（准定常叶素 quasi-steady blade element）
=====================================================================================
几何：从 STEP 抠出的 FL 腿销孔位置（米，x 前进 / z 向上，机体坐标系）
    O1 上舵机轴 (servo_2)  → 曲柄1 (28 mm) → 链接杆 (35 mm) → 摇臂 B–O2–C 绕 O2 转   = 髋摆动 ψ
    O2 下舵机轴 (servo_1)  → 曲柄2 (39 mm) 驱动平行四边形 O2–C–E–D               = 伸展/蜷缩
    约束：平行四边形始终在摇臂左侧 —— 相对角 φ = θ2 − ψ ∈ (0°, 180°)；φ→0 桨板折起贴到摇臂上（蜷缩），
          φ≈90–110° 桨板伸展在最下方。整套机构严格在 x–z 平面内运动（脚蹼只沿 ±y 张开）。
    桨板刚性连接在 E、D 上，桨轴与摇臂杆平行（平行四边形保证），从 E 起 50 mm 细杆 + 42 mm 桨头
脚蹼：收拢态 = CAD 桨头（30 mm 宽）；展开态 = 以桨头肩部为扇心、半径 42 mm 的扇形（展开角 FAN_DEG 可调）
      划水相（桨向后扫）由线拉开；回收相无拉力自动收拢。
步态：由 (周期 T, 髋摆幅, 抬桨幅, 相位, 划水占比) 参数化 —— 拿到真实舵机曲线后直接替换 gait()。
力模型：每个叶素  dF = ½ρ C_n |w_n| w_n n̂ dA  +  附加质量 ρπw²/4 · a_n n̂ dr ；忽略切向摩擦。
输出：轨迹图、推力/阻力时间历程、周期平均净推力（U = 0 / 0.1 / 0.2 m/s）、动画、展开/收拢桨 STL。

    python3 leg_dynamics.py                      # 默认参数
    T=1.25 FAN_DEG=120 SWEEP=60 PHI_EXT=105 PHI_RET=40 python3 leg_dynamics.py
    PSI_START=120 SWEEP=60 python3 leg_dynamics.py     # 把行程居中到摇臂竖直（120→60°）
"""
import os, sys, struct
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
def _cjk_font():
    """找一个有中文的字体：系统 Noto CJK，或 WSL 里借 Windows 的微软雅黑/黑体"""
    import glob as _g, warnings
    cands = [f for f in font_manager.findSystemFonts() if "NotoSansCJK" in f and "Regular" in f]
    cands += _g.glob("/mnt/c/Windows/Fonts/msyh.ttc") + _g.glob("/mnt/c/Windows/Fonts/simhei.ttf") + _g.glob("/mnt/c/Windows/Fonts/simsun.ttc")
    cands += [f for f in font_manager.findSystemFonts() if any(k in f for k in ("wqy", "WenQuanYi", "SourceHanSans", "NotoSansSC"))]
    for f in cands:
        try:
            font_manager.fontManager.addfont(f); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); return True
        except Exception: pass
    warnings.filterwarnings("ignore", message="Glyph .* missing"); return False
_cjk_font()
plt.rcParams["axes.unicode_minus"] = False

E_ = lambda k, d: float(os.environ.get(k, d))
# ------------------------------------------------------------------ 参数
T_PER   = E_("T", 1.25)        # 划水周期 [s]（视频自相关估计 1.2–1.3 s）
SWEEP   = E_("SWEEP", 60.0)    # 髋摆动总幅度 [deg]（视频目测 50–60°）
PHI_EXT = E_("PHI_EXT", 105.0) # 划水相平行四边形张角 φ=θ2−ψ [deg]（CAD 位姿 109°；90° = 矩形）
PHI_RET = E_("PHI_RET", 40.0)  # 回收相蜷缩张角 [deg]（→0 = 桨板完全贴到摇臂上）；必须在 (PHI_MIN, PHI_MAX) 内
PHI_MIN, PHI_MAX = 12.0, 168.0 # 平行四边形只能在摇臂左侧：φ∈(0°,180°)，留出杆件厚度余量
DUTY    = E_("DUTY", 0.5)      # 划水相占周期比例
FAN_DEG = E_("FAN_DEG", 120.0) # 脚蹼展开角 [deg]
FAN_R   = E_("FAN_R", 0.042)   # 扇形半径 = 桨头长 [m]
CLOSED_W= E_("CLOSED_W", 0.030) # 收拢态桨头宽度 [m]（CAD 30 mm；真实折扇收拢可能只有 8–12 mm）
FAN_CLOSE=E_("FAN_CLOSE", -1)  # 脚蹼开始收拢的相位 tau（默认 -1 = 划水结束 DUTY；设 0.3 表示划水中途就松线）
VEL_LAW = os.environ.get("VEL_LAW", "cos")   # 相内速度律：cos = 余弦（起止为零的正弦速度）；trap = 梯形（RAMP 比例加减速 + 匀速，同峰值速率下冲量更大）
RAMP    = E_("RAMP", 0.15)                    # trap 的加/减速段占该相时长的比例
SWITCH_FRAC=E_("SWITCH_FRAC", 0.3)  # 回收相里蜷缩 / 伸展切换各占的时间比例（0.3 = 前 30 % 收、后 30 % 放；0.5 = 慢蜷缩，收完立刻放）
U_LIST  = [0.0, 0.1, 0.2]      # 机体前进速度 [m/s]
RHO, NU = 998.8, 1.09e-6
CN_HEAD, CN_FAN, CD_STEM = 1.2, 1.3, 1.0     # 法向力系数：桨头 / 展开扇 / 细杆
NT = 400                                    # 每周期时间步
OUT = os.environ.get("OUT", os.path.dirname(os.path.abspath(__file__)))

# ------------------------------------------------------------------ 机构几何（FL，米）
O1 = np.array([0.0847, 0.0260]); O2 = np.array([0.0634, 0.0040])
L_CRANK1, L_ROD = 0.028, 0.035
L_ROCK_B = 0.0301            # O2→B（摇臂长端）
L_ROCK_C = 0.0147            # O2→C（摇臂短端 = 平行四边形地杆）
L_PAR    = 0.039             # 曲柄2 = 从动杆（平行四边形长边，理想化取均值 39 mm）
L_TOP    = 0.0147            # 桨顶 E–D = 地杆
STEM, HEAD = 0.050, 0.042    # 桨：细杆长、桨头长
STEM_W, THK = 0.003, 0.0056
HEAD_W = CLOSED_W
PSI0   = np.deg2rad(E_("PSI_START", 136.0))   # 划水起点的摇臂角 O2→C（CAD 位姿 136°；摇臂竖直 = 90° 时桨板竖直、推力方向最正）
TH2_0  = np.deg2rad(-115.0)  # 曲柄2 O2→D 名义角
BLADE_OFF = np.pi            # 桨轴方向 = 摇臂 C→O2 方向（与 O2→C 反向）

def rot(a): return np.array([np.cos(a), np.sin(a)])

def rocker_from_crank1(th1):
    """四杆 O1–A–B–O2：给曲柄1角 th1，解摇臂角（O2→B 方向）"""
    A = O1 + L_CRANK1 * rot(th1)
    d = np.linalg.norm(A - O2);
    # B 在圆(O2, L_ROCK_B) 与 圆(A, L_ROD) 交点
    a = (L_ROCK_B**2 - L_ROD**2 + d**2) / (2 * d)
    h2 = L_ROCK_B**2 - a**2
    if h2 < 0: return None, A
    h = np.sqrt(h2); u = (A - O2) / d; p = O2 + a * u
    Bs = [p + h * np.array([-u[1], u[0]]), p - h * np.array([-u[1], u[0]])]
    B = min(Bs, key=lambda b: b[1])          # 取下方解（CAD 位姿 B 在 O2 右下）
    return np.arctan2(*(B - O2)[::-1]), A, B

def crank1_from_rocker(psi_B):
    """逆解：给定摇臂角（O2→B），求曲柄1角（用于报告舵机指令）"""
    B = O2 + L_ROCK_B * rot(psi_B)
    d = np.linalg.norm(B - O1)
    a = (L_CRANK1**2 - L_ROD**2 + d**2) / (2 * d); h2 = L_CRANK1**2 - a**2
    if h2 < 0: return np.nan
    h = np.sqrt(h2); u = (B - O1) / d; p = O1 + a * u
    As = [p + h * np.array([-u[1], u[0]]), p - h * np.array([-u[1], u[0]])]
    A = min(As, key=lambda a_: abs(np.arctan2(*(a_ - O1)[::-1]) - np.deg2rad(-35)))
    return np.arctan2(*(A - O1)[::-1])

def pose(psi_C, phi):
    """给摇臂角 psi_C（O2→C）与平行四边形相对角 phi（曲柄2 相对摇臂），返回关节点与桨轴方向。
       曲柄2 绝对角 th2 = psi_C + phi（servo_1 指令）。phi∈(0,π) 保证 D 在摇臂左侧。"""
    th2 = psi_C + phi
    C = O2 + L_ROCK_C * rot(psi_C)
    D = O2 + L_PAR * rot(th2)
    E = C + (D - O2)                       # 平行四边形闭合：CE ∥ O2D, ED ∥ O2C
    B = O2 + L_ROCK_B * rot(psi_C + np.pi) # 摇臂是过 O2 的直杆
    u = rot(psi_C + BLADE_OFF)             # 桨轴（从 E 指向桨尖），与摇臂平行
    n = np.array([-u[1], u[0]])            # 桨面法向（x–z 平面内，⟂ 桨轴）
    return dict(B=B, C=C, D=D, E=E, u=u, n=n, th2=th2)

# ------------------------------------------------------------------ 步态
def gait(t):
    """返回 (psi_C, th2, power) —— 摇臂角、曲柄2角、是否处于划水相。
       划水相：桨向后扫（psi 减小）。回收相：向前扫并抬起。正弦 + 占比 DUTY 的非对称时间映射。"""
    tau = (t / T_PER) % 1.0
    # 相位映射：前 DUTY 为划水（0→π），其余为回收（π→2π）
    ph = np.where(tau < DUTY, np.pi * tau / DUTY, np.pi + np.pi * (tau - DUTY) / (1 - DUTY))
    sw = np.deg2rad(SWEEP) / 2
    if VEL_LAW == "trap":                                 # 梯形速度律：相内归一化位移 s(u)，u = 相内进度
        r = min(max(RAMP, 0.02), 0.5)
        def s_of(u):
            u = np.clip(u, 0, 1)
            return np.where(u < r, u ** 2 / (2 * r), np.where(u <= 1 - r, u - r / 2, 1 - r - (1 - u) ** 2 / (2 * r))) / (1 - r)
        up = tau / DUTY; ur = (tau - DUTY) / (1 - DUTY)
        psi = np.where(tau < DUTY, PSI0 - 2 * sw * s_of(up), PSI0 - 2 * sw * (1 - s_of(ur)))
    else:
        psi = PSI0 - sw + sw * np.cos(ph)                # ph=0: 最前(PSI0)，ph=π: 最后(PSI0-SWEEP)
    # 平行四边形张角：划水相保持伸展 PHI_EXT；回收相蜷缩到 PHI_RET（回收开始收、回收结束前伸开，各占回收相 30%）
    pe, pr = np.deg2rad(PHI_EXT), np.deg2rad(PHI_RET)
    rec = (tau - DUTY) / (1 - DUTY)                       # 回收相进度 0..1（划水相为负）
    sf = min(max(SWITCH_FRAC, 0.05), 0.5)
    blend = np.where(rec < 0, 0.0, np.where(rec < sf, rec / sf, np.where(rec < 1 - sf, 1.0, (1 - rec) / sf)))
    blend = 0.5 - 0.5 * np.cos(np.pi * np.clip(blend, 0, 1))   # 平滑
    phi = pe + (pr - pe) * blend
    if np.any(phi < np.deg2rad(PHI_MIN)) or np.any(phi > np.deg2rad(PHI_MAX)):
        sys.exit(f"φ 超出 ({PHI_MIN}°,{PHI_MAX}°)：平行四边形不能越过摇臂")
    power = tau < DUTY
    return psi, phi, power

# ------------------------------------------------------------------ 桨的叶素离散
def blade_strips(fan_open):
    """返回沿桨轴的叶素：(s 起点, ds, 宽度 w, 法向力系数, 是否附加质量)；s 从 E 起"""
    strips = []
    ns = 10
    for i in range(ns):                                         # 细杆
        s0 = STEM * i / ns; strips.append((s0, STEM / ns, STEM_W, CD_STEM))
    nh = 20
    for i in range(nh):                                         # 桨头 / 扇
        s0 = STEM + HEAD * i / nh; sm = s0 + HEAD / nh / 2; r = sm - STEM
        if fan_open:
            w = 2 * r * np.sin(np.deg2rad(FAN_DEG) / 2) + STEM_W
            strips.append((s0, HEAD / nh, w, CN_FAN))
        else:
            # CAD 桨头：肩部 ~25 mm 内圆弧过渡到 30 mm
            w = HEAD_W * min(1.0, np.sqrt(max(r, 1e-6) / 0.025)) if r < 0.025 else HEAD_W
            strips.append((s0, HEAD / nh, w, CN_HEAD))
    return strips

def fan_fraction(t, power):
    """脚蹼展开度 0..1：划水开始由线拉开（过渡 8% 周期），在 FAN_CLOSE 相位松线自动收拢（默认划水结束）"""
    tau = (t / T_PER) % 1.0; tr = 0.08
    tc = DUTY if FAN_CLOSE < 0 else FAN_CLOSE
    f = np.zeros_like(tau)
    f = np.where(tau < tr, tau / tr, f)                                  # 划水开始：拉开
    f = np.where((tau >= tr) & (tau < tc), 1.0, f)                       # 全开
    f = np.where((tau >= tc) & (tau < tc + tr), 1 - (tau - tc) / tr, f)  # 松线：收拢
    return np.clip(f, 0, 1)

# ------------------------------------------------------------------ 主计算
def simulate(U, ncycle=2):
    dt = T_PER / NT; t = np.arange(0, ncycle * T_PER, dt)
    psi, phi, power = gait(t); ff = fan_fraction(t, power)
    P = [pose(p, q) for p, q in zip(psi, phi)]
    th2 = np.array([p["th2"] for p in P])
    E = np.array([p["E"] for p in P]); u = np.array([p["u"] for p in P]); n = np.array([p["n"] for p in P])
    tip = E + (STEM + HEAD) * u
    F = np.zeros((len(t), 2)); Fam = np.zeros((len(t), 2)); M_O2 = np.zeros(len(t))
    strips_o, strips_c = blade_strips(True), blade_strips(False)
    # 位置历史 → 速度/加速度（中心差分，周期性）
    def pos_of(s):  return E + s * u
    for (s0, ds, w_o, cn_o), (_, _, w_c, cn_c) in zip(strips_o, strips_c):
        sm = s0 + ds / 2
        X = pos_of(sm)
        V = (np.roll(X, -1, 0) - np.roll(X, 1, 0)) / (2 * dt)
        A_ = (np.roll(X, -1, 0) - 2 * X + np.roll(X, 1, 0)) / dt**2
        w = ff * w_o + (1 - ff) * w_c; cn = ff * cn_o + (1 - ff) * cn_c
        wrel = np.stack([-U - V[:, 0], -V[:, 1]], 1)            # 水相对叶素的速度
        wn = np.einsum("ij,ij->i", wrel, n)
        dF = (0.5 * RHO * cn * np.abs(wn) * wn * (w * ds))[:, None] * n           # 法向压差力
        an = np.einsum("ij,ij->i", A_, n)
        dFa = (-RHO * np.pi * w**2 / 4 * ds * an)[:, None] * n                    # 附加质量（相对静水）
        F += dF + dFa; Fam += dFa
        r = X - O2; M_O2 += r[:, 0] * (dF + dFa)[:, 1] - r[:, 1] * (dF + dFa)[:, 0]
    # 只取最后一个周期做统计
    m = t >= (ncycle - 1) * T_PER
    return dict(t=t, psi=psi, phi=phi, th2=th2, power=power, ff=ff, E=E, tip=tip, u=u, F=F, Fam=Fam, M=M_O2, m=m, P=P)

def main():
    os.makedirs(OUT, exist_ok=True)
    res = {U: simulate(U) for U in U_LIST}
    r0 = res[0.0]; m = r0["m"]; t = r0["t"][m] - r0["t"][m][0]
    th1 = np.array([crank1_from_rocker(p + np.pi) for p in r0["psi"][m]])   # 摇臂 O2→B 角 = psi_C + π

    # ---- 报告
    print(f"步态：T = {T_PER} s ({1/T_PER:.2f} Hz)，髋摆幅 {SWEEP:.0f}°，平行四边形 φ 伸展 {PHI_EXT:.0f}° / 蜷缩 {PHI_RET:.0f}°，划水占比 {DUTY:.2f}，脚蹼展开角 {FAN_DEG:.0f}°（半径 {FAN_R*1000:.0f} mm）")
    tipv = np.linalg.norm(np.gradient(r0["tip"][m], t, axis=0), axis=1)
    print(f"桨尖速度：最大 {tipv.max():.2f} m/s，划水相平均 {tipv[r0['power'][m]].mean():.2f} m/s；桨尖轨迹 x∈[{r0['tip'][m][:,0].min():.3f},{r0['tip'][m][:,0].max():.3f}] z∈[{r0['tip'][m][:,1].min():.3f},{r0['tip'][m][:,1].max():.3f}]")
    Re = tipv.max() * HEAD_W / NU; print(f"桨头 Re（基于桨尖速度与 30 mm 宽度）≈ {Re:.1e}")
    print(f"舵机指令范围：servo_2(曲柄1) {np.degrees(np.nanmin(th1)):.0f}°…{np.degrees(np.nanmax(th1)):.0f}°   servo_1(曲柄2) {np.degrees(r0['th2'][m].min()):.0f}°…{np.degrees(r0['th2'][m].max()):.0f}°")
    print("\n单腿周期平均力（+x = 推力 thrust，机体坐标）：")
    print(f"{'U [m/s]':>8} {'净推力 mN':>10} {'净推力 mN':>13}")
    print(f"{'':>8} {'(压差项)':>10} {'(含附加质量)':>13}")
    for U, r in res.items():
        mm = r["m"]; Fx = r["F"][mm, 0]; Fd = Fx - r["Fam"][mm, 0]; pw = r["power"][mm]
        print(f"{U:8.2f} {Fd.mean()*1e3:10.1f} {Fx.mean()*1e3:13.1f} | 划水相压差峰 {Fd[pw].max()*1e3:6.1f}  回收相压差峰 {Fd[~pw].min()*1e3:6.1f}  附加质量周期均值 {r['Fam'][mm,0].mean()*1e3:6.1f}  峰值|M_O2| {np.abs(r['M'][mm]).max()*1e3:5.1f} mN·m")

    # ---- 图 1：轨迹 + 力历程
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    ax = axs[0]
    for k in np.linspace(0, m.sum() - 1, 12).astype(int):
        Pk = r0["P"][np.where(m)[0][k]]; col = "#c0392b" if r0["power"][m][k] else "#2980b9"
        ax.plot([O2[0], Pk["C"][0], Pk["E"][0], Pk["D"][0], O2[0]], [O2[1], Pk["C"][1], Pk["E"][1], Pk["D"][1], O2[1]], color=col, alpha=.35, lw=1)
        tp = Pk["E"] + (STEM + HEAD) * Pk["u"]; ax.plot([Pk["E"][0], tp[0]], [Pk["E"][1], tp[1]], color=col, alpha=.6, lw=2)
    ax.plot(r0["tip"][m][:, 0], r0["tip"][m][:, 1], "k-", lw=1.2, label="桨尖轨迹")
    ax.plot(*O1, "ko"); ax.plot(*O2, "ko"); ax.text(O1[0] + .003, O1[1], "O1"); ax.text(O2[0] + .003, O2[1], "O2")
    ax.axhline(-0.010, color="tab:blue", ls="--", lw=1); ax.text(0.10, -0.006, "水线", color="tab:blue")
    ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_xlabel("x [m] → 前"); ax.set_ylabel("z [m]")
    ax.set_title("一个周期的腿姿态（红=划水相·伸展·脚蹼张开，蓝=回收相·蜷缩·收拢）"); ax.legend(loc="lower right")
    ax = axs[1]
    for U, r in res.items():
        mm = r["m"]; ax.plot(t / T_PER, r["F"][mm, 0] * 1e3, lw=2, label=f"U = {U:.1f} m/s")
    ax.plot(t / T_PER, r0["Fam"][m, 0] * 1e3, "k:", lw=1, label="附加质量项 (U=0)")
    ax.axhline(0, color="k", lw=.8); ax.axvspan(0, DUTY, color="#c0392b", alpha=.08); ax.text(DUTY / 2, ax.get_ylim()[1] * .9, "划水相", ha="center", color="#c0392b")
    ax.set_xlabel("t / T"); ax.set_ylabel("F_x [mN]  (+ 推力)"); ax.grid(alpha=.3); ax.legend(); ax.set_title("单腿 x 向力时间历程（准定常叶素）")
    ax = axs[2]
    ax.plot(t / T_PER, np.degrees(r0["psi"][m] - PSI0), label="髋摆动 Δψ [°]")
    ax.plot(t / T_PER, np.degrees(r0["phi"][m]), label="平行四边形张角 φ [°]（小=蜷缩）")
    ax.plot(t / T_PER, np.degrees(r0["th2"][m]), "-.", label="servo_1 绝对角 θ2 [°]")
    ax.plot(t / T_PER, r0["ff"][m] * 100, "--", label="脚蹼展开度 [%]")
    ax.plot(t / T_PER, tipv * 100, ":", label="桨尖速度 [cm/s]")
    ax.set_xlabel("t / T"); ax.grid(alpha=.3); ax.legend(); ax.set_title("步态输入与脚蹼状态")
    fig.suptitle(f"BODY2 单腿动力模型 · 第 1 级（准定常）  T={T_PER}s  摆幅{SWEEP:.0f}°  φ {PHI_EXT:.0f}°→{PHI_RET:.0f}°  扇{FAN_DEG:.0f}°", fontsize=13)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "leg_quasi_steady.png"), dpi=95); plt.close()

    # ---- 图 2：净推力对关键不对称参数的敏感性（U = 0，单腿）
    fig, axs = plt.subplots(1, 4, figsize=(20, 4.6))
    g = globals()
    def sweep_param(name, vals, ax, label, xt=None):
        keep = g[name]; yd, yt = [], []
        for v in vals:
            g[name] = v
            if name == "CLOSED_W": g["HEAD_W"] = v
            r = simulate(0.0); mm = r["m"]
            yd.append((r["F"][mm, 0] - r["Fam"][mm, 0]).mean() * 1e3); yt.append(r["F"][mm, 0].mean() * 1e3)
        g[name] = keep
        if name == "CLOSED_W": g["HEAD_W"] = keep
        xs = xt if xt is not None else vals
        ax.plot(xs, yd, "o-", label="压差项 (drag-based)"); ax.plot(xs, yt, "s--", label="含附加质量")
        ax.axhline(0, color="k", lw=.8); ax.set_xlabel(label); ax.set_ylabel("单腿周期平均净推力 [mN]，U=0"); ax.grid(alpha=.3); ax.legend(fontsize=8)
    sweep_param("DUTY", [0.3, 0.35, 0.4, 0.45, 0.5, 0.6], axs[0], "划水相占周期比例 DUTY（<0.5 = 划得快、收得慢）")
    sweep_param("CLOSED_W", [0.008, 0.012, 0.018, 0.024, 0.030], axs[1], "收拢态桨头宽度 [mm]", xt=[8, 12, 18, 24, 30])
    sweep_param("FAN_CLOSE", [0.2, 0.25, 0.3, 0.35, 0.4, 0.5], axs[2], "松线收拢相位 tau（0.5 = 划水结束才松）")
    sweep_param("PHI_RET", [20, 40, 60, 80, 105], axs[3], "回收相蜷缩张角 φ_ret [°]（105 = 不蜷缩）")
    fig.suptitle(f"敏感性（准定常叶素，U = 0）：基准 DUTY={DUTY} 收拢宽={CLOSED_W*1000:.0f}mm 松线相位={'划水结束' if FAN_CLOSE<0 else FAN_CLOSE} 扇={FAN_DEG:.0f}° φ_ret={PHI_RET:.0f}°", fontsize=12)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, "leg_sensitivity.png"), dpi=95); plt.close()

    # ---- 动画帧
    fr_dir = os.path.join(OUT, "anim"); os.makedirs(fr_dir, exist_ok=True)
    idx = np.where(m)[0][::5]
    for j, k in enumerate(idx):
        fig, ax = plt.subplots(figsize=(7, 6.5)); Pk = r0["P"][k]; pw = r0["power"][k]; f = r0["ff"][k]
        col = "#c0392b" if pw else "#2980b9"
        # 机构
        A_ = O1 + L_CRANK1 * rot(crank1_from_rocker(r0["psi"][k] + np.pi))
        ax.plot([O1[0], A_[0], Pk["B"][0]], [O1[1], A_[1], Pk["B"][1]], "-", color="#fd8d3c", lw=4, solid_capstyle="round")
        ax.plot([Pk["B"][0], Pk["C"][0]], [Pk["B"][1], Pk["C"][1]], "-", color="#31a354", lw=5, solid_capstyle="round")
        ax.plot([Pk["C"][0], Pk["E"][0]], [Pk["C"][1], Pk["E"][1]], "-", color="#636363", lw=4, solid_capstyle="round")
        ax.plot([O2[0], Pk["D"][0]], [O2[1], Pk["D"][1]], "-", color="#e6550d", lw=4, solid_capstyle="round")
        # 桨：细杆 + 桨头（宽度用线宽示意：收拢 30 mm，展开按扇角）
        sh = Pk["E"] + STEM * Pk["u"]; tp = Pk["E"] + (STEM + HEAD) * Pk["u"]
        ax.plot([Pk["E"][0], sh[0]], [Pk["E"][1], sh[1]], "-", color="#a63603", lw=3)
        wtip = (f * (2 * HEAD * np.sin(np.deg2rad(FAN_DEG) / 2)) + (1 - f) * HEAD_W)
        nrm = Pk["n"]
        # 用多边形画桨头在 x-z 面的“厚度”与在 y 向展开的示意（透视：宽度画成沿法向的阴影）
        poly = np.array([sh - nrm * 0.003, tp - nrm * 0.003, tp + nrm * 0.003, sh + nrm * 0.003])
        ax.fill(poly[:, 0], poly[:, 1], color="#a63603", alpha=.9)
        ax.text(tp[0] + .004, tp[1] - .004, f"脚蹼{'张开' if f > .5 else '收拢'}  宽 {wtip*1000:.0f} mm", color=col, fontsize=9)
        # 力
        Fk = r0["F"][k]; sc = 0.25
        ax.arrow(sh[0], sh[1], Fk[0] * sc, Fk[1] * sc, width=.0012, color=col, length_includes_head=True)
        ax.text(0.02, -0.088, f"t/T = {((r0['t'][k]/T_PER)%1):.2f}   {'划水相 (power)' if pw else '回收相 (recovery)'}   φ = {np.degrees(r0['phi'][k]):.0f}°\nF_x = {Fk[0]*1e3:+.1f} mN   F_z = {Fk[1]*1e3:+.1f} mN", fontsize=10,
                bbox=dict(fc="white", ec=col))
        ax.axhline(-0.010, color="tab:blue", ls="--", lw=1)
        ax.plot(r0["tip"][m][:, 0], r0["tip"][m][:, 1], "k-", lw=.6, alpha=.4)
        ax.set_xlim(-0.02, 0.14); ax.set_ylim(-0.10, 0.05); ax.set_aspect("equal"); ax.grid(alpha=.25)
        ax.set_title(f"BODY2 单腿划水 · 运动学 + 准定常力（U = 0）"); ax.set_xlabel("x [m] → 前"); ax.set_ylabel("z [m]")
        plt.tight_layout(); plt.savefig(os.path.join(fr_dir, f"f_{j:03d}.png"), dpi=80); plt.close()
    print(f"\n动画帧 {len(idx)} 张 → {fr_dir}/   (ffmpeg -framerate 16 -i f_%03d.png -pix_fmt yuv420p leg_anim.mp4)")

    # ---- 展开 / 收拢桨 STL（桨自身坐标：轴沿 -z，从 E 向下；面法向 x；y 为展开方向）
    def write_stl(tri, fn):
        nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
        rec = np.zeros(len(tri), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])); rec["n"] = nrm; rec["v"] = tri
        open(fn, "wb").write(b"\0" * 80 + struct.pack("<I", len(tri)) + rec.tobytes())
    def extrude(outline, thk):
        """outline: (N,2) 闭合多边形 (y, s)；挤出厚度 thk 沿 x；返回三角形 (M,3,3) 坐标 (x, y, z=-s)"""
        P = np.array(outline); N = len(P); tris = []
        c = P.mean(0)
        for i in range(N):                       # 两个端面（扇形/多边形为星形可从质心三角化）
            a, b = P[i], P[(i + 1) % N]
            for x, flip in ((thk / 2, False), (-thk / 2, True)):
                t3 = [[x, c[0], -c[1]], [x, a[0], -a[1]], [x, b[0], -b[1]]]
                tris.append(t3[::-1] if flip else t3)
            # 侧面
            p1, p2 = [thk / 2, a[0], -a[1]], [thk / 2, b[0], -b[1]]; q1, q2 = [-thk / 2, a[0], -a[1]], [-thk / 2, b[0], -b[1]]
            tris.append([p1, q1, q2]); tris.append([p1, q2, p2])
        return np.array(tris)
    # 收拢态：细杆 + CAD 桨头（简化为肩部圆弧 + 矩形）
    s = np.linspace(0, 1, 30)
    shoulder = [(HEAD_W / 2 * np.sin(np.pi / 2 * k), STEM + 0.025 * (1 - np.cos(np.pi / 2 * k))) for k in s]
    outl = [(-STEM_W / 2, 0), (STEM_W / 2, 0), (STEM_W / 2, STEM)] + [(w, z) for w, z in shoulder if w > STEM_W / 2] + \
           [(HEAD_W / 2, STEM + HEAD), (-HEAD_W / 2, STEM + HEAD)] + [(-w, z) for w, z in shoulder[::-1] if w > STEM_W / 2] + [(-STEM_W / 2, STEM)]
    write_stl(extrude(outl, THK), os.path.join(OUT, "paddle_closed.stl"))
    # 展开态：细杆 + 扇形（扇心在肩部 s=STEM）
    ang = np.deg2rad(FAN_DEG); th = np.linspace(ang / 2, -ang / 2, 40)   # 从右角到左角，避免多边形自交
    fan = [(FAN_R * np.sin(a), STEM + FAN_R * np.cos(a)) for a in th]
    outl_o = [(-STEM_W / 2, 0), (STEM_W / 2, 0), (STEM_W / 2, STEM)] + fan + [(-STEM_W / 2, STEM)]
    write_stl(extrude(outl_o, THK), os.path.join(OUT, "paddle_open.stl"))
    A_c = STEM * STEM_W + HEAD_W * HEAD - (1 - np.pi / 4) * HEAD_W * 0.025 / 2
    A_o = STEM * STEM_W + 0.5 * FAN_R**2 * ang
    print(f"桨 STL：paddle_closed.stl（面积≈{A_c*1e4:.1f} cm²）  paddle_open.stl（扇 {FAN_DEG:.0f}°，尖端宽 {2*FAN_R*np.sin(ang/2)*1000:.0f} mm，面积≈{A_o*1e4:.1f} cm²，×{A_o/A_c:.2f}）")

if __name__ == "__main__":
    main()
