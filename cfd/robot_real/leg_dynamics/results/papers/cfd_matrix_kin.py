#!/usr/bin/env python3
"""CFD 算例矩阵的共用运动学（沉浮 heave + 俯仰 pitch）。
只算运动学，不算受力 —— 受力交给 CFD。
"""
import numpy as np

RHO = 1000.0
C   = 0.042          # 弦长 chord [m]（展开态折扇半径）
THK = 0.0056
U   = 0.25           # 来流 [m/s]
SERVO_RATE = 820.0   # 7465W @8.4 V 带载 [°/s]

wrap = lambda x: (x + np.pi) % (2*np.pi) - np.pi


def _alpha_series(tt, om, h0, th0, psi, u=U):
    h  = h0*np.sin(om*tt)
    hd = h0*om*np.cos(om*tt)
    gW = np.arctan2(-hd, -u*np.ones_like(tt))
    th = th0*np.sin(om*tt + psi)
    return h, hd, th, wrap(gW - th - np.pi), gW


def solve_th0(om, h0, psi, al_max, u=U):
    """给定目标最大攻角，反解俯仰幅值 θ0（羽化）。"""
    tt = np.linspace(0, 2*np.pi/om, 721)
    ind0 = np.arctan(h0*om/u)
    lo, hi = 0.0, ind0*1.2
    f = lambda x: np.abs(_alpha_series(tt, om, h0, x, psi, u)[3]).max() - al_max
    if f(hi) > 0:   # 极端情况兜底
        return hi
    for _ in range(60):
        mid = 0.5*(lo + hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    return 0.5*(lo + hi)


def shape_plateau(x, k=2.6):
    """平台型（近方波）α 时程 —— Hover 2004 / Read 2003 的整形波形。"""
    return np.tanh(k*np.sin(x))/np.tanh(k)


class Case:
    """一个算例：给定 St、α_max、相位 ψ、俯仰轴、波形。"""
    def __init__(self, tag, label, St=0.35, al_max=20.0, psi=90.0, pivot=0.25,
                 wave="sin", le="sq", c=C, h0c=0.83, u=U, color="#2a78d6", note=""):
        self.tag, self.label, self.note, self.color = tag, label, note, color
        self.St, self.psi_d, self.pivot, self.wave, self.le = St, psi, pivot, wave, le
        self.c, self.u = c, u
        self.h0 = h0c*c
        self.f  = St*u/(2*self.h0)
        self.T  = 1/self.f
        self.om = 2*np.pi*self.f
        self.psi = np.radians(psi)
        self.al_max = np.radians(al_max)
        self.ind0 = np.arctan(self.h0*self.om/u)

        tt = np.linspace(0, self.T, 1441)
        self.th0 = solve_th0(self.om, self.h0, self.psi, self.al_max, u)
        h, hd, th, al, gW = _alpha_series(tt, self.om, self.h0, self.th0, self.psi, u)

        if wave in ("plateau", "alpha_sin"):
            t_pk = tt[np.argmax(al)]
            x = self.om*(tt - t_pk) + np.pi/2
            al = self.al_max*(shape_plateau(x) if wave == "plateau" else np.sin(x))
            th = np.unwrap(wrap(gW - al - np.pi))

        self.tt, self.h, self.hd, self.th, self.al = tt, h, hd, np.unwrap(th), al
        self.gW = gW
        self.chi = self.th0/self.ind0
        self.rate_pk = np.degrees(np.abs(np.gradient(self.th, tt)).max())

    def at(self, tau):
        """tau ∈ [0,1) → 插值状态"""
        i = tau*(len(self.tt) - 1)
        g = lambda a: np.interp(i, np.arange(len(a)), a)
        return dict(h=g(self.h), th=g(self.th), al=g(self.al), hd=g(self.hd))


def foil_poly(P, th, c, pivot, thk=THK):
    d = np.array([np.cos(th), np.sin(th)])
    n = np.array([-d[1], d[0]])*thk/2
    LE = P + pivot*c*d
    TE = P - (1 - pivot)*c*d
    return np.array([LE + n, TE + n, TE - n, LE - n]), LE, TE


# ---------------- 七个算例（按 08 号文献笔记定）
BLUE, RED, GOLD, GREEN, PUR, TEAL = "#2a78d6", "#c2185b", "#d4a017", "#16a085", "#7b3f9d", "#0f7b8a"

SCENES = [
    dict(no=1, title="① 验证算例 —— 复现 Anderson 1998",
         why="目标 η = 87 %、C_T ≈ 0.67。这一关过不了，后面所有数都没有基准。",
         cases=[Case("V", "NACA 0012 · AR 6 · 俯仰轴 1/3 弦",
                     St=0.30, al_max=20.0, psi=75.0, pivot=1/3, h0c=0.75, color=RED,
                     note="h₀/c = 0.75，ψ = 75°")]),
    dict(no=2, title="② St 扫描（α_max 固定 20°，靠 χ 调）",
         why="效率最优 0.25–0.35 vs 纯推力最优 0.60。两端差一倍多，先决定优化目标。",
         cases=[Case(f"St{s:.2f}", f"St = {s:.2f}", St=s, color=c)
                for s, c in zip((0.25, 0.35, 0.50, 0.60), (BLUE, GREEN, GOLD, RED))]),
    dict(no=3, title="③ α_max 扫描（St 固定 0.35）",
         why="文献窗口 15–25°，中心 20°。低于 10° 推力太小，高于 35° 前缘涡分离、效率崩。",
         cases=[Case(f"a{a:.0f}", f"α_max = {a:.0f}°", St=0.35, al_max=a, color=c)
                for a, c in zip((12, 20, 30), (BLUE, GREEN, RED))]),
    dict(no=4, title="④ 沉浮–俯仰相位 ψ 扫描",
         why="海豚的实测相位文献里查不到。90° 是从翼型实验搬来的假设，必须当变量扫。",
         cases=[Case(f"p{p:.0f}", f"ψ = {p:.0f}°", St=0.35, psi=p, color=c)
                for p, c in zip((75, 90, 105), (BLUE, GREEN, PUR))]),
    dict(no=5, title="⑤ 俯仰轴位置",
         why="Anderson 与 Read 全套数据都锚在 1/3 弦，但没人做过系统扫描。1/4 弦更靠前、气动更稳。",
         cases=[Case(f"x{int(x*100)}", f"轴 {x:.2f} c", St=0.35, pivot=x, color=c)
                for x, c in zip((0.25, 1/3), (BLUE, GOLD))]),
    dict(no=6, title="⑥ α 波形三选一（在 St = 0.50 上比）",
         why="正弦俯仰在高 χ 下 α 已自然接近平台，真正的对照是强制正弦 α。",
         cases=[Case("wth", "正弦俯仰 θ（现状）", St=0.50, wave="sin", color=BLUE),
                Case("wal", "强制正弦 α", St=0.50, wave="alpha_sin", color=GREEN),
                Case("wsq", "方波型 α", St=0.50, wave="plateau", color=RED)]),
    dict(no=7, title="⑦ 前缘：方边平板 vs 倒圆（同一套运动学）",
         why="海豚尾鳍前缘半径 0.032–0.046 c → 折算到 42 mm 弦是 1.3–2.8 mm。方边会过早分离。",
         cases=[Case("le_sq", "方边平板", St=0.35, le="sq", color=BLUE),
                Case("le_rd", "倒圆 R 1.7 mm", St=0.35, le="round", color=TEAL)]),
]

if __name__ == "__main__":
    for sc in SCENES:
        print(f"\n[{sc['no']}] {sc['title']}")
        for c in sc["cases"]:
            print(f"   {c.label:<32s} f = {c.f:5.2f} Hz  θ0 = {np.degrees(c.th0):5.1f}°  "
                  f"χ = {c.chi:.2f}  α_max = {np.degrees(np.abs(c.al).max()):4.1f}°  "
                  f"俯仰峰速 {c.rate_pk:5.0f} °/s  余量 {SERVO_RATE/c.rate_pk:4.1f}×")
