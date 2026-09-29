#!/usr/bin/env python3
"""
最优轨迹动态图：V3trap 运动学（梯形律，T 0.72 s，DUTY 0.382，摆角 50°，收拢 12 mm）+ CFD 两态合成力（开 τ_o = −0.06，闭 τ_c = 0.33）
  python3 make_optimal_anim.py [--out DIR] [--cycles 2] [--slow 5] [--fps 30] [--to -0.06 --tc 0.33] [--ref V3]
输出：optimal_gait_cfd_anim.mp4（放慢 slow 倍）、optimal_gait_cfd_storyboard.png（8 帧分镜）、frames_opt/（帧 PNG）
左：机构侧视（曲柄 1 + 连杆 + 摇臂 + 平行四边形 + 桨），CFD 合力箭头（在桨头中心），桨尖轨迹；右上小图：桨头俯视（扇的张开度，过渡 60 ms 示意）
右：F_x(τ) 与 V3 参考、桨头速度 + 开度、O2 力矩 / 功率；底：事件时间线
"""
import os, sys, json, argparse, subprocess, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=HERE); ap.add_argument("--cycles", type=float, default=2.0); ap.add_argument("--slow", type=float, default=5.0); ap.add_argument("--fps", type=int, default=30)
ap.add_argument("--to", type=float, default=-0.06); ap.add_argument("--tc", type=float, default=0.33); ap.add_argument("--scheme", default="V3trap"); ap.add_argument("--ref", default="V3"); ap.add_argument("--tfan", type=float, default=0.060, help="开合过渡时间 [s]（示意）")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
from dd_common import *
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle, FancyArrowPatch
from PIL import Image

d = load(os.path.join(LD, "results_verify", a.scheme)); ld = d["ld"]; T = d["T"]; DUTY = d["DUTY"]; tau = d["tau"]
c = composite(d, a.to, a.tc); Fx, Fz, My, P = c["F"][:, 0] * 1e3, c["F"][:, 2] * 1e3, c["M"][:, 1] * 1e3, c["P"] * 1e3
Fxs, Fzs, Mys, Ps = smooth(Fx, tau, 0.01), smooth(Fz, tau, 0.01), smooth(My, tau, 0.01), smooth(P, tau, 0.01)
Fx_mean, P_mean, My_peak = Fx.mean(), P.mean(), np.abs(Mys).max()
dref = load(os.path.join(LD, "results_verify", a.ref)); cref = composite(dref, 0.0, dref["DUTY"]); Fx_ref = smooth(cref["F"][:, 0] * 1e3, dref["tau"], 0.01); Fx_ref_mean = cref["F"][:, 0].mean() * 1e3
vP = np.linalg.norm(d["vP"], axis=1)
def interp(y, x): return np.interp(x % 1.0, tau, y, period=1.0)
tf = a.tfan / T
def fan_frac(x):
    """张开度 0..1：τ_o 开始张（tf 内完成），τ_c 开始收（tf 内完成）；两态合成本身是瞬时切换，这里只是过渡示意"""
    x = (x - a.to) % 1.0; L = (a.tc - a.to) % 1.0
    if x < tf: return x / tf
    if x < L: return 1.0
    if x < L + tf: return 1 - (x - L) / tf
    return 0.0
def phase_name(x):
    x %= 1.0; r = ld.RAMP
    if x < DUTY * r: return "划水 · 加速", "#c2185b"
    if x < DUTY * (1 - r): return "划水 · 恒速", "#c2185b"
    if x < DUTY: return "划水 · 减速（已收扇）", "#c2185b"
    if x < DUTY + (1 - DUTY) * r: return "回收 · 加速 + 蜷缩", "#2a78d6"
    if x < 1 - (1 - DUTY) * r: return "回收 · 恒速", "#2a78d6"
    return "回收 · 减速（已开扇）", "#2a78d6"
ev = sorted([(a.to % 1.0, "开扇指令"), ((a.to + tf) % 1.0, "全开"), (0.0, "划水开始"), (DUTY * ld.RAMP, "达恒速"), (DUTY * (1 - ld.RAMP), "减速开始"), (a.tc, "收扇指令"), (DUTY, "回收开始"), ((a.tc + tf), "全闭"), (1 - (1 - DUTY) * ld.RAMP, "回收减速")])

# ---- 一个周期的运动学与轨迹（画背景轨迹用）
psi_c, phi_c, _ = ld.gait(tau * T); PP = [ld.pose(p, q) for p, q in zip(psi_c, phi_c)]
tip_path = np.array([p["E"] + (ld.STEM + ld.HEAD) * p["u"] for p in PP]); head_path = np.array([p["E"] + (ld.STEM + ld.HEAD / 2) * p["u"] for p in PP])
O1, O2 = np.array(ld.O1), np.array(ld.O2)

nfr = int(round(a.cycles * T * a.slow * a.fps)); taus = np.linspace(0, a.cycles, nfr, endpoint=False)
fr_dir = os.path.join(a.out, "frames_opt"); os.makedirs(fr_dir, exist_ok=True)
for f in glob.glob(os.path.join(fr_dir, "*.png")): os.remove(f)
RED, BLUE, GOLD = "#c2185b", "#2a78d6", "#d4a017"

def draw_frame(i, x):
    fig = plt.figure(figsize=(16, 9), facecolor=SURF)
    gs = fig.add_gridspec(4, 2, width_ratios=[1.15, 1], height_ratios=[1.25, 1, 1, 0.42], left=0.045, right=0.975, top=0.885, bottom=0.06, hspace=0.45, wspace=0.17)
    axL = fig.add_subplot(gs[:3, 0]); axF = fig.add_subplot(gs[0, 1]); axV = fig.add_subplot(gs[1, 1]); axM = fig.add_subplot(gs[2, 1]); axT = fig.add_subplot(gs[3, :])
    for ax in (axL, axF, axV, axM, axT): style(ax)
    xm = x % 1.0; t_ms = xm * T * 1e3; pw = xm < DUTY; col = RED if pw else BLUE; f = fan_frac(xm); pname, pcol = phase_name(xm)
    psi, phi, _ = ld.gait(np.array([xm * T])); psi, phi = float(psi[0]), float(phi[0]); Pk = ld.pose(psi, phi)
    A_ = O1 + ld.L_CRANK1 * ld.rot(ld.crank1_from_rocker(psi + np.pi))
    E, u, n = Pk["E"], Pk["u"], Pk["n"]; sh = E + ld.STEM * u; tp = E + (ld.STEM + ld.HEAD) * u; hc = E + (ld.STEM + ld.HEAD / 2) * u
    # ---------- 左：机构
    axL.axhspan(-0.12, -0.010, color="#dbe9f6", alpha=0.55, lw=0, zorder=0); axL.axhline(-0.010, color="#4f86c6", lw=1.2, ls=(0, (5, 3))); axL.text(-0.018, -0.0085, "水线 z = −10 mm", color="#4f86c6", fontsize=9)
    axL.plot(tip_path[:, 0], tip_path[:, 1], color=MUTED, lw=0.8, alpha=0.6, zorder=1)
    j = int(xm * len(tau)); trail = np.arange(j - int(0.12 * len(tau)), j + 1) % len(tau); axL.plot(tip_path[trail, 0], tip_path[trail, 1], color=col, lw=2.2, alpha=0.8, zorder=2)
    axL.plot([O1[0], A_[0], Pk["B"][0]], [O1[1], A_[1], Pk["B"][1]], "-", color="#fd8d3c", lw=5, solid_capstyle="round", zorder=3)
    axL.plot([Pk["B"][0], Pk["C"][0]], [Pk["B"][1], Pk["C"][1]], "-", color="#31a354", lw=6, solid_capstyle="round", zorder=3)
    axL.plot([Pk["C"][0], E[0]], [Pk["C"][1], E[1]], "-", color="#636363", lw=5, solid_capstyle="round", zorder=3)
    axL.plot([O2[0], Pk["D"][0]], [O2[1], Pk["D"][1]], "-", color="#e6550d", lw=5, solid_capstyle="round", zorder=3)
    axL.plot([E[0], Pk["D"][0]], [E[1], Pk["D"][1]], "-", color="#9e9e9e", lw=3, solid_capstyle="round", zorder=3)
    for p_, nm in ((O1, "O1"), (O2, "O2")): axL.plot(*p_, "o", color=INK, ms=7, zorder=5); axL.text(p_[0] + 0.003, p_[1] + 0.003, nm, fontsize=9, color=INK)
    axL.plot([E[0], sh[0]], [E[1], sh[1]], "-", color="#7b3f00", lw=3.5, zorder=4)
    thk = 0.0022 + 0.0033 * f                                      # 侧视里用“厚度”示意张开（真实张开在 y 向）
    poly = np.array([sh - n * thk, tp - n * thk, tp + n * thk, sh + n * thk]); axL.fill(poly[:, 0], poly[:, 1], color=col, alpha=0.9, zorder=4)
    Fk = np.array([interp(Fxs, xm), interp(Fzs, xm)]) * 1e-3; sc = 0.10
    if np.linalg.norm(Fk) > 5e-3: axL.add_patch(FancyArrowPatch(hc, hc + Fk * sc, arrowstyle="simple,head_width=9,head_length=9,tail_width=3", color=GOLD, zorder=6, mutation_scale=1))
    axL.text(hc[0] + Fk[0] * sc + 0.004, hc[1] + Fk[1] * sc, f"F = {np.linalg.norm(Fk)*1e3:.0f} mN", color=GOLD, fontsize=10, fontweight="bold", zorder=6)
    axL.set_xlim(-0.025, 0.145); axL.set_ylim(-0.115, 0.055); axL.set_aspect("equal"); axL.set_xlabel("x [m] → 前进方向", color=INK2); axL.set_ylabel("z [m]", color=INK2)
    axL.set_title(f"t = {t_ms:5.0f} ms   τ = {xm:.3f}   {pname}", loc="left", fontsize=12.5, color=pcol, pad=8)
    box = (f"F_x = {interp(Fxs, xm):+6.0f} mN   F_z = {interp(Fzs, xm):+5.0f} mN\nM_y(O2) = {interp(Mys, xm):+5.1f} mN·m   P = {interp(Ps, xm):5.0f} mW\n"
           f"桨头 |v| = {interp(vP, xm):.2f} m/s   ψ = {np.degrees(psi):.0f}°   φ = {np.degrees(phi):.0f}°\n脚蹼 {'张开' if f > 0.99 else ('收拢' if f < 0.01 else ('张开中' if fan_frac(xm-0.01) < f else '收拢中'))}  ({100*f:.0f} %)")
    axL.text(0.02, 0.03, box, transform=axL.transAxes, fontsize=10, color=INK, va="bottom", bbox=dict(fc="white", ec=GRID, boxstyle="round,pad=0.4"))
    # 俯视小图：扇的张开度
    ins = axL.inset_axes([0.66, 0.70, 0.32, 0.28]); ins.set_facecolor("white"); ins.set_xlim(-45, 45); ins.set_ylim(-8, 48); ins.set_aspect("equal"); ins.axis("off")
    ins.add_patch(Rectangle((-1.5, -8), 3, 8, color="#7b3f00"))
    ang = 12 + (120 - 12) * f
    if f > 0.02: ins.add_patch(Wedge((0, 0), 42, 90 - ang / 2, 90 + ang / 2, color=col, alpha=0.9))
    ins.add_patch(Rectangle((-6, 0), 12, 42, color=col if f <= 0.02 else "#00000000", ec=col, lw=1.2))
    ins.text(0, 45, f"桨头俯视：{'扇 120°' if f > 0.99 else ('收拢 12 mm' if f < 0.01 else f'过渡 {100*f:.0f} %')}", ha="center", fontsize=8.5, color=INK2)
    # ---------- 右上：F_x
    axF.axvspan(0, DUTY, color=RED, alpha=0.05, lw=0); axF.plot(dref["tau"], Fx_ref, color=MUTED, lw=1.1, ls=(0, (4, 3)), label=f"参考 {a.ref}（余弦 · 30 mm）：F̄x {Fx_ref_mean:.1f} mN")
    axF.plot(tau, Fxs, color="#e8c6d3", lw=1.2); k = tau <= xm + 1e-9; axF.plot(tau[k], Fxs[k], color=RED, lw=2.2, label=f"最优（梯形 · 12 mm · 开 {a.to:+.2f} / 闭 {a.tc:.2f}）：F̄x {Fx_mean:.1f} mN")
    axF.axhline(Fx_mean, color=RED, lw=1, ls=(0, (3, 3))); axF.axhline(0, color=GRID, lw=1); axF.plot(xm, interp(Fxs, xm), "o", color=GOLD, ms=9, mec=INK, zorder=6)
    axF.set_ylabel("推力 F_x [mN]", color=INK2); axF.set_xlim(0, 1); axF.set_ylim(-260, 520); axF.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK2)
    axF.set_title("CFD 合成推力（第 3 周期，open/closed 两态按开合时刻拼接）", loc="left", fontsize=10.5, color=INK, pad=6)
    # ---------- 右中：速度 + 开度
    axV.axvspan(0, DUTY, color=RED, alpha=0.05, lw=0); axV.plot(tau, vP, color=INK2, lw=1.8, label="桨头中心速度 |v_P|（梯形律：41 ms 斜坡 + 恒速）"); axV.plot(xm, interp(vP, xm), "o", color=GOLD, ms=8, mec=INK, zorder=6)
    ff = np.array([fan_frac(t_) for t_ in tau]); ax2 = axV.twinx(); ax2.fill_between(tau, 0, ff, color=RED, alpha=0.12, lw=0); ax2.plot(tau, ff, color=RED, lw=1.2); ax2.set_ylim(0, 3.2); ax2.set_yticks([0, 1]); ax2.set_yticklabels(["闭", "开"], color=RED, fontsize=9); ax2.tick_params(colors=RED)
    for s_ in ("top", "right"): ax2.spines[s_].set_visible(False)
    axV.set_ylabel("|v_P| [m/s]", color=INK2); axV.set_xlim(0, 1); axV.set_ylim(0, 0.42); axV.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK2)
    axV.set_title("桨头速度与脚蹼开度（开度过渡 60 ms 为示意）", loc="left", fontsize=10.5, color=INK, pad=6)
    # ---------- 右下：力矩 + 功率
    axM.axvspan(0, DUTY, color=RED, alpha=0.05, lw=0); axM.plot(tau, Mys, color="#7a5cc0", lw=1.6, label=f"M_y(O2)，峰 {My_peak:.1f} mN·m"); axM.plot(xm, interp(Mys, xm), "o", color=GOLD, ms=8, mec=INK, zorder=6)
    axM.set_ylabel("M_y [mN·m]", color="#7a5cc0"); axM.set_ylim(-45, 25); axM.set_xlim(0, 1)
    ax3 = axM.twinx(); ax3.plot(tau, Ps, color="#16a085", lw=1.4, label=f"功率 P，均值 {P_mean:.1f} mW"); ax3.set_ylabel("P [mW]", color="#16a085"); ax3.set_ylim(-45, 140); ax3.tick_params(colors="#16a085")
    for s_ in ("top", "right"): ax3.spines[s_].set_visible(False)
    h1, l1 = axM.get_legend_handles_labels(); h2, l2 = ax3.get_legend_handles_labels(); axM.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK2)
    axM.set_title("O2 力矩与机构输入功率", loc="left", fontsize=10.5, color=INK, pad=6); axM.set_xlabel("周期相位 τ = t / T", color=INK2)
    for ax in (axF, axV, axM): ax.axvline(xm, color=GOLD, lw=1, alpha=0.8)
    # ---------- 底：时间线
    axT.set_xlim(0, 1); axT.set_ylim(0, 1); axT.set_yticks([]); axT.spines["left"].set_visible(False)
    axT.axvspan(0, DUTY, color=RED, alpha=0.18, lw=0); axT.axvspan(DUTY, 1, color=BLUE, alpha=0.14, lw=0)
    axT.text(DUTY / 2, 0.78, f"划水 {DUTY*T*1e3:.0f} ms", ha="center", color=RED, fontsize=9.5); axT.text((1 + DUTY) / 2, 0.78, f"回收 {(1-DUTY)*T*1e3:.0f} ms", ha="center", color=BLUE, fontsize=9.5)
    for k_, (te, nm) in enumerate(ev):
        axT.axvline(te, color=INK2, lw=0.8); axT.text(te, 0.03 + 0.25 * (k_ % 3), f"{nm} {te*T*1e3:.0f} ms", ha="left" if te < 0.9 else "right", va="bottom", fontsize=7.6, color=INK2)
    axT.axvline(xm, color=GOLD, lw=2.5); axT.set_xlabel(f"周期 T = {T:.2f} s（{1/T:.2f} Hz）   动画放慢 {a.slow:.0f}×", color=INK2, fontsize=9.5)
    fig.suptitle(f"BODY2 单腿最优轨迹 · 运动学 + CFD 力   梯形速度律 · 摆角 50°（ψ 120°→70°）· 收拢 12 mm · 开扇 τ = {a.to:+.2f} / 收扇 τ = {a.tc:.2f}   →  F̄x = {Fx_mean:.0f} mN（V3 的 {Fx_mean/Fx_ref_mean:.2f} 倍），P̄ = {P_mean:.0f} mW",
                 fontsize=13, x=0.045, ha="left", y=0.975, color=INK)
    fig.text(0.045, 0.945, "合成：τ ∈ [τ_o, τ_c) 取 open 算例的力，其余取 closed；切换假定瞬时（折扇过渡本身未模拟）。箭头 = 桨上的 CFD 合力（x–z），比例 0.1 m/N。U = 0 系泊。", fontsize=9, color=MUTED)
    fig.savefig(os.path.join(fr_dir, f"f_{i:04d}.png"), dpi=100); plt.close(fig)

for i, x in enumerate(taus):
    draw_frame(i, x)
    if i % 30 == 0: print(f"  帧 {i}/{nfr}", flush=True)
# 末尾停 1.5 s
last = os.path.join(fr_dir, f"f_{nfr-1:04d}.png")
for k in range(int(1.5 * a.fps)):
    os.link(last, os.path.join(fr_dir, f"f_{nfr+k:04d}.png"))
mp4 = os.path.join(a.out, "optimal_gait_cfd_anim.mp4")
subprocess.run(["ffmpeg", "-y", "-v", "error", "-framerate", str(a.fps), "-i", os.path.join(fr_dir, "f_%04d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart", mp4], check=True)
print("→", mp4)
# 分镜：一个周期内 8 帧，只裁左侧机构面板（整幅仪表盘缩小后看不清）
sel = [int(round(t_ / a.cycles * nfr)) % nfr for t_ in np.linspace(0, 1, 9)[:-1]]
crop = (40, 62, 790, 768); cw, chh = crop[2] - crop[0], crop[3] - crop[1]
grid = Image.new("RGB", (4 * cw, 2 * chh), "white")
for k, i in enumerate(sel):
    im = Image.open(os.path.join(fr_dir, f"f_{i:04d}.png")).crop(crop); grid.paste(im, ((k % 4) * cw, (k // 4) * chh))
sb = os.path.join(a.out, "optimal_gait_cfd_storyboard.png"); grid.save(sb); print("→", sb)
# 另存一张完整仪表盘（τ = 0.20，恒速划水中）作报告插图
full = os.path.join(a.out, "optimal_gait_cfd_frame.png"); Image.open(os.path.join(fr_dir, f"f_{int(round(0.20 / a.cycles * nfr)):04d}.png")).save(full); print("→", full)
