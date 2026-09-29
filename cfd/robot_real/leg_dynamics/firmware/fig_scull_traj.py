import csv, math, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f_ in [f_ for f_ in font_manager.findSystemFonts() if "NotoSansCJK" in f_ and "Regular" in f_]:
    try: font_manager.fontManager.addfont(f_); plt.rcParams["font.family"] = font_manager.FontProperties(fname=f_).get_name(); break
    except Exception: pass
plt.rcParams["axes.unicode_minus"] = False
import scull_trajectory as st
rows = list(csv.DictReader(open("scull_traj_100ppc.csv")))
ph = np.array([float(r["phase"]) for r in rows]); q1 = np.array([float(r["q1"]) for r in rows]); q2 = np.array([float(r["q2"]) for r in rows])
yJ = np.array([float(r["yJ"]) for r in rows]); zJ = np.array([float(r["zJ"]) for r in rows])
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 6))
a.plot(yJ, zJ, "b-", lw=2, label="足端 J 轨迹（圆弧）")
for i, c in ((25, "#c0392b"), (0, "#555"), (75, "#1f77b4")):
    H, E, I, J = st.fk(q1[i], q2[i]); B = (0, st.AB); C = B + st.BC * np.array([math.cos(math.radians(q1[i])), math.sin(math.radians(q1[i]))])
    a.plot([0, H[0]], [0, H[1]], color=c, lw=2); 
    a.plot([H[0], I[0]], [H[1], I[1]], color=c, lw=2); a.plot([I[0], J[0]], [I[1], J[1]], color=c, lw=3)
    a.text(J[0] + 3, J[1], f"θ = {q2[i]:.0f}°", color=c, fontsize=9)
a.axhline(-39, color="gray", ls=":"); a.text(-20, -37, "船底（髋轴下 39）", color="gray", fontsize=8)
a.axhline(7, color="#2a9d8f", ls="--"); a.text(-20, 9, "压载后水线", color="#2a9d8f", fontsize=8)
a.plot(0, 0, "ko"); a.text(2, 3, "A 髋轴", fontsize=8)
a.set_aspect("equal"); a.grid(alpha=.3); a.set_xlabel("y 向前 / mm"); a.set_ylabel("z 向上 / mm"); a.set_title("腿平面：杆 −18° ± 18°，长杆竖直（红 顶端 / 灰 中 / 蓝 底端）", fontsize=10); a.legend(loc="lower left", fontsize=8)
b.plot(ph, q1, label="q1（曲柄）"); b.plot(ph, q2, label="q2（杆 = 尾鳍中立方向）"); b.set_xlabel("相位"); b.set_ylabel("角度 / °"); b.grid(alpha=.3); b.legend(); b.set_title("关节角（1.5 Hz 峰值 q1 194 °/s，q2 170 °/s）", fontsize=10)
plt.tight_layout(); plt.savefig("fig_scull_traj.png", dpi=140); print("ok")
