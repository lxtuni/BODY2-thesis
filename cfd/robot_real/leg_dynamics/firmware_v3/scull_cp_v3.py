"""新固件（7 控制点接口）下的摇橹式尾鳍轨迹：大腿竖直不动 γ = −90°（E 固定在 (0,−30)），杆角 q_r 正弦摆动
   M(φ) = E − 45·dir(q_r)，q_r = QR_MEAN + QR_AMP·sin(2πφ)。7 点 = 同一段圆弧上下往返的折线：cp0 顶端 → cp3 底端 → (cp4 = cp3, 3 % 停顿) → 回到 cp0。
   直线段 + swim_smooth_enable = 0 （段内匀速），段时长按正弦分配。"""
import math, json, csv, numpy as np
from leg_v3 import ik, d, O2E, EM
QR_MEAN, QR_AMP = 40.0, 25.0        # q_r 15 … 65°：顶端杆朝船尾方向下倾 15°（尾鳍内侧尖离船底 ≥ 8 mm），底端下倾 65°
GAMMA = -85.0                       # 大腿近竖直；−90 会让 LF θ2 = 6°、RF θ2 = 170°，贴着 PWM 0–180 边界，取 −85 留 5° 余量
CAL = dict(LF=(53.0, 1.0, 96.0, 1.0), RF=(130.0, -1.0, 80.0, -1.0), LH=(59.0, 1.0, 90.0, 1.0), RH=(127.0, -1.0, 90.0, -1.0))
E = O2E * d(GAMMA)
def M_of(qr): return E - EM * d(qr)
# 7 点：按相位 φ_s（正弦相位）取样，cp0 = 顶端（sin = −1）
phis = [0.75, 0.75 + 1/6, 0.75 + 2/6, 1.25, 1.25, 1.25 + 1/6, 1.25 + 2/6]
qrs = [QR_MEAN + QR_AMP * math.sin(2 * math.pi * p) for p in phis]
pts = [M_of(q) for q in qrs]
DWELL = 0.03
knots = [0.0, 1/6, 2/6, 0.5, 0.5 + DWELL, 0.5 + DWELL + (0.5 - DWELL) / 3, 0.5 + DWELL + 2 * (0.5 - DWELL) / 3]
print("控制点（xM, yM）与逆解：")
rows = []
for i, (p, q, ph) in enumerate(zip(pts, qrs, phis)):
    s = ik(p[0], p[1])
    rows.append(dict(cp=i, phase=round(knots[i], 4), xM=round(p[0], 1), yM=round(p[1], 1), qr_target=q, **{k: s[k] for k in ("q1", "gamma", "qr", "ACO2", "O2EM", "r")}))
    servo = {leg: (o1 + s1 * s["q1"], o2 + s2 * s["gamma"]) for leg, (o1, s1, o2, s2) in CAL.items()}
    print(f" cp{i} @{knots[i]:.3f}: ({p[0]:6.1f}, {p[1]:6.1f})  q_r={s['qr']:5.1f} q1={s['q1']:6.1f} γ={s['gamma']:7.1f} |O2M|={s['r']:4.1f} ∠ACO2={s['ACO2']:5.1f} ∠O2EM={s['O2EM']:5.1f}  LF θ1/θ2={servo['LF'][0]:.0f}/{servo['LF'][1]:.0f} RF={servo['RF'][0]:.0f}/{servo['RF'][1]:.0f}")
# 连续轨迹（100 点）核对：折线 + 段内匀速 vs 理想正弦
N = 100; err = []
for k in range(N):
    ph = k / N; seg = max(i for i in range(7) if ph >= knots[i]); nxt = knots[seg + 1] if seg < 6 else 1.0
    t = (ph - knots[seg]) / (nxt - knots[seg]); a, b = pts[seg], pts[(seg + 1) % 7]; m = a + (b - a) * t
    s = ik(m[0], m[1]); err.append((m, s))
q1s = np.array([s["q1"] for _, s in err]); gms = np.array([s["gamma"] for _, s in err])
print(f"\n全程：q1 {q1s.min():.1f}…{q1s.max():.1f}°  γ {gms.min():.1f}…{gms.max():.1f}°  ∠ACO2 min {min(s['ACO2'] for _, s in err):.1f}°  ∠O2EM min {min(s['O2EM'] for _, s in err):.1f}°")
for f in (1.0, 1.5):
    dt = 1 / (f * N); print(f" {f} Hz 峰值角速度：q1 {np.abs(np.gradient(q1s, dt)).max():.0f} °/s（折线匀速段的均值，正弦峰值≈{QR_AMP*2*math.pi*f:.0f}）")
R = 62.2 - 15.0
print(f"\n尾鳍枢轴（离 E 47 mm，若 E 是外侧孔则 62）：顶端 y = {E[1] + 47*math.sin(math.radians(QR_MEAN-QR_AMP-180)):.1f} / {E[1] + 62*math.sin(math.radians(QR_MEAN-QR_AMP-180)):.1f}（船底在 O2 系 −39，内侧尖处 −34）；底端 y = {E[1] + 47*math.sin(math.radians(QR_MEAN+QR_AMP-180)):.1f}；横向沉浮 ±{47*math.sin(math.radians(QR_AMP)):.0f} / ±{62*math.sin(math.radians(QR_AMP)):.0f} mm")
json.dump(dict(qr_mean=QR_MEAN, qr_amp=QR_AMP, gamma=GAMMA, knots=knots, points=[[round(p[0], 1), round(p[1], 1)] for p in pts], rows=rows), open("scull_cp_v3.json", "w"), indent=1, default=float)
with open("scull_cp_v3_100pts.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["phase", "xM", "yM", "q1", "gamma", "q_rocker", "ACO2", "O2EM"])
    for k, (m, s) in enumerate(err): w.writerow([k / N] + [round(v, 3) for v in (m[0], m[1], s["q1"], s["gamma"], s["qr"], s["ACO2"], s["O2EM"])])
