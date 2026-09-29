#!/usr/bin/env python3
"""
M3 平滑船体 —— 数值设计部分（与 CAD 无关，可单独运行画线图）
坐标：零件坐标系（密封身体0826.STEP）：x 沿船长（尖端 −x、钝端 +x），y 向上（对接面/水线 y = −10，顶缘 y = +15，原平底 y = −35），z 横向。
      与 CFD 坐标的关系：x_CFD = −x，z_CFD = y，y_CFD = z。
水下型线（M3）：水线轮廓 b_w(x) 直接取自 0826 外表面（B 样条面的 y = −10 截面，本身光滑）；
      龙骨深度 D(x) = T3·(1 − ξ)^p_b·(1 + ξ)^p_p·(1 + a·ξ)，ξ = (x − x_c)/L_h（T3 由排水量 385.1 cm³ 定、a 由浮心纵向位置与原船相同定；
      p_b = 0.8（钝端）、p_p = 0.85（尖端）使首尾龙骨线在水线处的倾角与上部外飘的首尾轮廓相接、无折角）；
      横剖面 b(t) = b_w·(1 − k·t)·√(1 − t²)，t = 深度分数（0 水线、1 龙骨）；k = s·D/b_w 使水线处斜率与上部外飘 s(x) 连续（C¹ 无折角）。
"""
import numpy as np, json, os
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
HERE = os.path.dirname(os.path.abspath(__file__))
P_B, P_P, KMAX, WL, DECK = 0.8, 0.85, 0.5, -10.0, 15.0     # 钝端 / 尖端指数
V_TARGET, LCB_TARGET = 385.1e3, 11.34          # mm³, mm（原船 y = −10 以下）

def _half(y, O):
    P = np.array(O[str(float(y))]); P = P[P[:, 1] > 1e-6]; P = P[np.argsort(P[:, 0])]
    xs, idx = np.unique(np.round(P[:, 0], 4), return_index=True); return xs, P[idx, 1]

class Design:
    def __init__(self, outlines_json=os.path.join(HERE, "orig_outlines.json")):
        O = json.load(open(outlines_json)); self.O = O
        xw, bw = _half(WL, O); x9, b9 = _half(WL + 1, O)
        self.xmin, self.xmax = float(xw.min()), float(xw.max())
        self.fb = UnivariateSpline(xw, bw, k=3, s=0); self.f9 = UnivariateSpline(x9, b9, k=3, s=0)
        self.xc = 0.5 * (self.xmin + self.xmax); self.Lh = 0.5 * (self.xmax - self.xmin)
        self.T3, self.a = 35.0, 0.0; self.solve()
    def bw(self, x): return np.clip(self.fb(np.clip(x, self.xmin, self.xmax)), 0, None)
    def slope(self, x):                       # 上部外飘在水线处的 db/dy（1 mm 差分）
        return np.clip(self.f9(np.clip(x, self.xmin, self.xmax)) - self.bw(x), 0, None)
    def xi(self, x): return (np.asarray(x) - self.xc) / self.Lh
    def D(self, x, T3=None, a=None):
        T3 = self.T3 if T3 is None else T3; a = self.a if a is None else a; xi = self.xi(x)
        return T3 * np.clip(1 - xi, 0, None) ** P_B * np.clip(1 + xi, 0, None) ** P_P * (1 + a * xi)
    def k(self, x, T3=None, a=None):
        b = self.bw(x); d = self.D(x, T3, a)
        return np.where(b > 1e-6, np.minimum(self.slope(x) * d / np.maximum(b, 1e-6), KMAX), 0.0)
    def area(self, x, T3=None, a=None):       # 剖面面积 = 2 ∫ b dy = 2 b_w D (π/4 − k/3)
        return 2 * self.bw(x) * self.D(x, T3, a) * (np.pi / 4 - self.k(x, T3, a) / 3)
    def volume_lcb(self, T3=None, a=None):
        x = np.linspace(self.xmin, self.xmax, 2001); A = self.area(x, T3, a)
        V = np.trapezoid(A, x); return V, np.trapezoid(A * x, x) / V
    def solve(self):
        def eq(p):
            V, L = self.volume_lcb(p[0], p[1]); return [V / V_TARGET - 1, (L - LCB_TARGET) / 10]
        sol = fsolve(eq, [35.0, 0.0]); self.T3, self.a = float(sol[0]), float(sol[1])
        self.V, self.LCB = self.volume_lcb()
    def section(self, x, n=41):
        """横剖面：从 (z=−b_w, y=−10) 经龙骨到 (z=+b_w, y=−10)；返回 (z, y) 数组"""
        t = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, n))          # 水线与龙骨附近加密
        b = self.bw(x) * (1 - self.k(x) * t) * np.sqrt(np.clip(1 - t ** 2, 0, None)); y = WL - self.D(x) * t
        z = np.concatenate([-b, b[::-1][1:]]); yy = np.concatenate([y, y[::-1][1:]]); return z, yy      # 左舷水线 → 龙骨 → 右舷水线
    def keel(self, x): return WL - self.D(x)
    def waterline(self, y, n=400):
        """y < −10 的水线：返回 (x, ±b) 半宽（无解处 b = 0）"""
        x = np.linspace(self.xmin, self.xmax, n); d = self.D(x); t = np.where(d > 1e-9, (WL - y) / np.maximum(d, 1e-9), 2.0)
        b = np.where(t <= 1, self.bw(x) * (1 - self.k(x) * np.clip(t, 0, 1)) * np.sqrt(np.clip(1 - np.clip(t, 0, 1) ** 2, 0, None)), 0.0)
        return x, b
    def deepest(self):
        x = np.linspace(self.xmin, self.xmax, 4001); d = self.D(x); i = int(np.argmax(d)); return float(x[i]), float(WL - d[i])
    def summary(self):
        xd, yd = self.deepest(); x = np.linspace(self.xmin, self.xmax, 2001)
        Aw = 2 * np.trapezoid(self.bw(x), x); Bwl = 2 * self.bw(x).max(); L = self.xmax - self.xmin; T = self.T3 * (1 + abs(self.a))
        return dict(LWL=L, BWL=Bwl, T3=self.T3, a=self.a, deepest_x=xd, deepest_y=yd, draft=WL - yd, V_cm3=self.V / 1e3, LCB=self.LCB,
                    Awp_cm2=Aw / 100, Cb=self.V / (L * Bwl * (WL - yd)), xmin=self.xmin, xmax=self.xmax)

if __name__ == "__main__":
    d = Design(); s = d.summary()
    for k, v in s.items(): print(f"{k:10s} {v:10.3f}")
    print("k(x) max", float(np.max(d.k(np.linspace(d.xmin + 0.5, d.xmax - 0.5, 500)))))
