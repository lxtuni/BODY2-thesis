#!/usr/bin/env python3
"""
深挖 C6 + C1 数据 → 最优轨迹（2026‑09‑16）
  python3 deep_dive_optimal.py [--out DIR]
输入：results_verify/{V3,V3c12,V3c06,V3trap}（open/closed 力与力矩时程 + gait.json），results_leg（现步态）
输出（--out，默认本目录）：
  fig_dd_timeseries.png   一个周期内 速度 / F_x / M_y / 功率 的时程：V3（余弦 30 mm）vs V3c12（余弦 12 mm）vs V3trap（梯形 12 mm）
  fig_dd_unsteady.png     划水相 CFD vs 准定常：非定常放大倍数随行程（弦数）衰减；累积冲量 余弦 vs 梯形
  fig_dd_recovery.png     回收阻力 vs 收拢宽度（CFD 三点 + 准定常预测）；有效阻力比 a 与 DUTY*
  fig_dd_timing_map.png   (τ_o, τ_c) 二维开合时刻图：V3c12 与 V3trap
  fig_dd_waterfall.png    从现步态到最优的逐项增益
  deep_dive_optimal.json  全部数字（报告与动画共用）
"""
import os, sys, json, argparse, importlib.util, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); LD = os.path.abspath(os.path.join(HERE, "..", ".."))
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=HERE); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
from dd_common import *
import matplotlib.pyplot as plt
runs = {k: os.path.join(LD, "results_verify", k) for k in ("V1", "V2", "V3", "V3c12", "V3c06", "V3trap")}; runs["current"] = os.path.join(LD, "results_leg")
D = {k: load(R) for k, R in runs.items() if os.path.exists(os.path.join(R, "force_open.dat"))}
OUTJ = dict(schemes={})
def summarize(k, to=None, tc=None):
    d = D[k]; to = 0.0 if to is None else to; tc = d["DUTY"] if tc is None else tc; c = composite(d, to, tc); tau = d["tau"]; T = d["T"]
    Fx, Fz, My, P = c["F"][:, 0], c["F"][:, 2], c["M"][:, 1], c["P"]
    st = tau < d["DUTY"]
    return dict(T=T, DUTY=d["DUTY"], to=to, tc=tc, Fx=Fx.mean() * 1e3, Fz=Fz.mean() * 1e3, P=P.mean() * 1e3, Fx_stroke=Fx[st].mean() * 1e3, Fx_rec=Fx[~st].mean() * 1e3,
                impulse=Fx.mean() * T * 1e3, impulse_stroke=np.trapezoid(Fx[st], tau[st] * T) * 1e3, My_peak=np.abs(smooth(My, tau)).max() * 1e3, P_peak=smooth(P, tau).max() * 1e3,
                Fx_peak=smooth(Fx, tau).max() * 1e3, Fx_min=smooth(Fx, tau).min() * 1e3, vP_max=np.linalg.norm(d["vP"], axis=1).max(), tip_max=np.linalg.norm(np.gradient(d["tip"], d["tt"], axis=0), axis=1).max(),
                open_only=d["F"]["open"][:, 0].mean() * 1e3, closed_only=d["F"]["closed"][:, 0].mean() * 1e3, closed_rec=d["F"]["closed"][~st, 0].mean() * 1e3, closed_stroke=d["F"]["closed"][st, 0].mean() * 1e3,
                open_stroke=d["F"]["open"][st, 0].mean() * 1e3, open_rec=d["F"]["open"][~st, 0].mean() * 1e3)
for k in D: OUTJ["schemes"][k] = summarize(k)
for k, s in OUTJ["schemes"].items(): print(f"{k:8s} Fx {s['Fx']:5.1f} (划水 {s['Fx_stroke']:6.1f} / 回收 {s['Fx_rec']:6.1f})  冲量 {s['impulse']:.1f} mN·s  My峰 {s['My_peak']:.1f}  P {s['P']:.1f} mW  closed-only 回收段均值 {s['closed_rec']:.1f}  open 划水段 {s['open_stroke']:.1f}")

# ------------------------------------------------------------------ 图 1：时程
fig, axs = plt.subplots(4, 1, figsize=(12.5, 12), sharex=True, facecolor=SURF, gridspec_kw=dict(hspace=0.12, left=0.08, right=0.98, top=0.93, bottom=0.06))
keys = ["V3", "V3c12", "V3trap"]
for ax in axs: style(ax)
for k in keys:
    d = D[k]; c = composite(d, 0.0, d["DUTY"]); tau = d["tau"]; col = COL[k]; lw = 2.0 if k == "V3trap" else 1.6
    v = np.linalg.norm(d["vP"], axis=1)
    axs[0].plot(tau, v, color=col, lw=lw, label=NAME[k] + f"   T = {d['T']:.2f} s, DUTY = {d['DUTY']:.3f}")
    axs[1].plot(tau, smooth(c["F"][:, 0], tau, 0.01) * 1e3, color=col, lw=lw, label=f"{NAME[k]}：F̄x = {c['F'][:,0].mean()*1e3:.1f} mN")
    axs[2].plot(tau, smooth(c["M"][:, 1], tau, 0.01) * 1e3, color=col, lw=lw, label=f"{NAME[k]}：|M_y| 峰 {np.abs(smooth(c['M'][:,1], tau)).max()*1e3:.1f} mN·m")
    axs[3].plot(tau, smooth(c["P"], tau, 0.01) * 1e3, color=col, lw=lw, label=f"{NAME[k]}：P̄ = {c['P'].mean()*1e3:.1f} mW，峰 {smooth(c['P'], tau).max()*1e3:.0f} mW")
    axs[1].axhline(c["F"][:, 0].mean() * 1e3, color=col, lw=0.9, ls=(0, (3, 3)), alpha=0.8)
for ax in axs:
    for k, ls in (("V3", (0, (5, 3))), ("V3trap", "-")):
        ax.axvline(D[k]["DUTY"], color=COL[k], lw=0.9, ls=ls, alpha=0.6)
    ax.axvspan(0, D["V3trap"]["DUTY"], color="#c2185b", alpha=0.035, lw=0)
axs[0].set_ylabel("桨头中心速度 |v_P| [m/s]", color=INK2); axs[1].set_ylabel("推力 F_x [mN]（合成：τ < DUTY 开，其余闭）", color=INK2)
axs[2].set_ylabel("O2 处力矩 M_y [mN·m]", color=INK2); axs[3].set_ylabel("机构输入功率 P [mW]", color=INK2); axs[3].set_xlabel("周期相位 τ = t / T（0 = 划水开始；竖线 = 各方案的 DUTY）", color=INK2)
axs[0].text(0.395, 0.30, "梯形律：41 ms 斜坡后桨头恒速 0.285 m/s（桨尖 0.40）\n余弦律：只在 τ ≈ 0.22 达到同样峰值", color="#c2185b", fontsize=9.5, va="bottom")
axs[1].text(0.47, -215, "梯形律减速段（τ 0.32–0.38，桨仍张开）：附加质量把推力拉到 −230 mN\n→ 这就是「减速一开始就收扇」值 +7.5 mN 的原因", color="#c2185b", fontsize=9.3, va="bottom")
axs[1].text(0.60, -120, "回收相：30 mm 收拢（蓝）比 12 mm（绿/红）多出的负推力 = 回收阻力节省", color=INK2, fontsize=9.3)
for ax in axs: ax.legend(frameon=False, fontsize=8.8, loc="upper right", labelcolor=INK2)
axs[0].set_xlim(0, 1)
fig.suptitle("一个周期的 CFD 时程（第 3 周期，两态合成于 τ = DUTY）：余弦 → 梯形速度律，收拢 30 → 12 mm", fontsize=13, x=0.08, ha="left", y=0.975, color=INK)
fig.text(0.08, 0.945, "F_x：+x = 推力；功率 P = −(F·v_E + M_E·ω)（对机构做功）；0.01 T 滑动平均", fontsize=9.5, color=MUTED)
plt.savefig(os.path.join(a.out, "fig_dd_timeseries.png"), dpi=110); plt.close()

# ------------------------------------------------------------------ 图 2：非定常放大倍数 + 累积冲量
fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.4), facecolor=SURF, gridspec_kw=dict(wspace=0.28, left=0.05, right=0.985, top=0.80, bottom=0.14))
for ax in axs: style(ax)
uns = {}
for k in ("V3", "V3trap"):
    d = D[k]; tau = d["tau"]; st = tau < d["DUTY"]; Fo = d["F"]["open"][:, 0]; Fq = qs_force(d, "open")[:, 0]
    v = np.linalg.norm(d["vP"], axis=1); s = np.concatenate([[0], np.cumsum(v[:-1] * np.diff(tau) * d["T"])])     # 桨头中心走过的路程
    c_chord = d["ld"].HEAD                                                                                    # 弦 = 桨头长 42 mm
    axs[0].plot(tau[st] / d["DUTY"], smooth(Fo, tau, 0.01)[st] * 1e3, color=COL[k], lw=2, label=f"{NAME[k].split('  ')[0]} CFD（open）")
    axs[0].plot(tau[st] / d["DUTY"], Fq[st] * 1e3, color=COL[k], lw=1.4, ls=(0, (4, 3)), label=f"{NAME[k].split('  ')[0]} 准定常 ½ρC_N A v_n²（C_N 1.3）")
    ok = st & (v > 0.12)
    ratio = smooth(Fo, tau, 0.02) / np.where(np.abs(Fq) > 1e-4, Fq, np.nan)          # CFD / 准定常(C_N 1.3) → 有效 C_N = 1.3 × ratio
    axs[1].plot(s[ok] / c_chord, 1.3 * ratio[ok], color=COL[k], lw=2, label=NAME[k].split("  ")[0] + f"：均值 C_N,eff = {1.3*np.nanmean(ratio[ok]):.1f}")
    Io = np.concatenate([[0], np.cumsum(Fo[st][:-1] * np.diff(tau[st]) * d["T"])]) * 1e3; Iq = np.concatenate([[0], np.cumsum(Fq[st][:-1] * np.diff(tau[st]) * d["T"])]) * 1e3
    axs[2].plot(tau[st] / d["DUTY"], Io, color=COL[k], lw=2, label=f"{NAME[k].split('  ')[0]} CFD：{Io[-1]:.1f} mN·s")
    axs[2].plot(tau[st] / d["DUTY"], Iq, color=COL[k], lw=1.4, ls=(0, (4, 3)), label=f"{NAME[k].split('  ')[0]} 准定常：{Iq[-1]:.1f} mN·s")
    uns[k] = dict(I_cfd=float(Io[-1]), I_qs=float(Iq[-1]), stroke_chords=float(s[st][-1] / c_chord), ratio_mean=float(np.nanmean(ratio[ok])),
                  ratio_first_chord=float(np.nanmean(ratio[ok & (s < c_chord)])), ratio_last=float(np.nanmean(ratio[ok & (s > s[st][-1] - c_chord)])))
axs[1].axhline(1.3, color=MUTED, lw=1, ls=(0, (3, 3))); axs[1].text(0.05, 1.38, "定常平板 C_N ≈ 1.3", color=MUTED, fontsize=8.5)
axs[1].axvspan(0, 2, color="#2a78d6", alpha=0.05, lw=0); axs[1].set_ylim(0, 12); axs[1].set_xlim(0, 2)
axs[1].text(1.0, 11.2, "涡形成阶段（≲ 2 弦，Kim & Gharib 2011）", color="#2a78d6", fontsize=9, ha="center")
axs[0].set_xlabel("划水相进度 τ / DUTY", color=INK2); axs[0].set_ylabel("F_x [mN]", color=INK2); axs[0].set_title("A  划水相推力：CFD（open 算例）vs 准定常 ½ρC_N A v_n²\n     两种速度律的 CFD 都是准定常的 ≈ 3.5–4 倍", loc="left", fontsize=11, color=INK, pad=9)
axs[1].set_xlabel("桨头中心走过的路程 s / c（c = 42 mm；|v_P| > 0.12 m/s 段）", color=INK2); axs[1].set_ylabel("有效法向力系数 C_N,eff = F / (½ρ A v_n²)", color=INK2)
axs[1].set_title("B  有效 C_N 随行程衰减：起动后 ≈ 5–6 → 1.5 弦后 ≈ 4（恒速）/ 2.5（余弦已在减速）\n     前 1 弦两种速度律重合 → 衰减由路程决定，不由速度律决定", loc="left", fontsize=11, color=INK, pad=9)
axs[2].set_xlabel("划水相进度 τ / DUTY", color=INK2); axs[2].set_ylabel("累积冲量 ∫F_x dt [mN·s]", color=INK2); axs[2].set_title(f"C  一次划水的冲量：梯形比余弦 +{100*(uns['V3trap']['I_cfd']/uns['V3']['I_cfd']-1):.0f} %（CFD）\n     准定常（同峰值速率、同 50° 摆角）预测 +{100*(uns['V3trap']['I_qs']/uns['V3']['I_qs']-1):.0f} % —— 两者一致", loc="left", fontsize=11, color=INK, pad=9)
axs[0].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower left"); axs[1].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="center right"); axs[2].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
fig.suptitle("梯形律为什么只多 18 %：同峰值速率下 ∫v² 本来就只多 20 %（摆角没变、划水时间缩短）；理论的 +60 % 要靠同时加大摆角 → C2", fontsize=13, x=0.05, ha="left", y=0.965, color=INK)
fig.text(0.05, 0.905, "准定常参考：桨头中心法向速度 v_n，A = 18.5 cm²（扇）+ 杆，C_N 1.3 / C_D 1.0，不含附加质量；CFD 力为 open 算例第 3 周期，0.01–0.02 T 滑动平均", fontsize=9.5, color=MUTED)
plt.savefig(os.path.join(a.out, "fig_dd_unsteady.png"), dpi=110); plt.close()
OUTJ["unsteady"] = uns; print("unsteady:", json.dumps(uns, indent=None))

# ------------------------------------------------------------------ 图 3：回收阻力 vs 收拢宽度；a 与 DUTY*
fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.2), facecolor=SURF, gridspec_kw=dict(wspace=0.3, left=0.05, right=0.985, top=0.80, bottom=0.14))
for ax in axs: style(ax)
W = {"V3": 30.0, "V3c12": 12.0, "V3c06": 5.6}; rec = {}; S_ = {k: OUTJ["schemes"][k]["Fx"] for k in W}
for k, w in W.items():
    d = D[k]; st = d["tau"] < d["DUTY"]; Fc = d["F"]["closed"][:, 0]
    Fq = qs_force(d, "closed", w * 1e-3)[:, 0]; v = np.linalg.norm(d["vP"], axis=1); vr = np.sqrt((v[~st] ** 2).mean()); vp = np.sqrt((v[st] ** 2).mean())
    rec[k] = dict(w=w, rec_mean=float(Fc[~st].mean() * 1e3), cyc_mean=float(Fc.mean() * 1e3), qs_rec=float(Fq[~st].mean() * 1e3), stroke_closed=float(Fc[st].mean() * 1e3),
                  a_raw=float(-Fc[~st].mean() / d["F"]["open"][st, 0].mean()), a=float(-Fc[~st].mean() / d["F"]["open"][st, 0].mean() * (vp / vr) ** 2), v_ratio=float(vr / vp))   # 有效阻力面积比 a = (F_rec/F_pw)·(v_p/v_r)²
ws = np.array([rec[k]["w"] for k in W]); fr = np.array([-rec[k]["rec_mean"] for k in W]); fq = np.array([-rec[k]["qs_rec"] for k in W])
kfit, bfit = np.polyfit(ws, fr, 1); wl = np.linspace(0, 32, 50)
axs[0].plot(wl, kfit * wl + bfit, color=INK2, lw=1.2, ls=(0, (4, 3)), label=f"线性拟合：{kfit:.2f} mN/mm × w + {bfit:.1f} mN")
axs[0].plot(ws, fq, "s-", color=MUTED, lw=1.2, ms=6, label="准定常（C_N 1.2 × w × 42 mm + 杆，无附加质量）")
for k in W: axs[0].plot(rec[k]["w"], -rec[k]["rec_mean"], "o", color=COL[k], ms=10, zorder=5, label=f"{k}  {rec[k]['w']:.0f} mm：{-rec[k]['rec_mean']:.1f} mN")
axs[0].set_xlim(0, 32); axs[0].set_ylim(0, None); axs[0].set_xlabel("收拢态桨头宽度 w [mm]", color=INK2); axs[0].set_ylabel("回收相平均阻力 −F̄x [mN]（closed 算例，V3 运动学）", color=INK2)
axs[0].set_title(f"A  回收阻力 ∝ 收拢宽度（{kfit:.2f} mN/mm，截距 ≈ 0：杆与附加质量不是地板）\n     12 mm 时回收拖累只剩周期平均的 {100*0.55*(-rec['V3c12']['rec_mean'])/S_['V3c12']:.0f} %，羽化到 5.6 mm 只再拿回 {(-rec['V3c12']['rec_mean'])-(-rec['V3c06']['rec_mean']):.1f} mN（相均值）→ 不值", loc="left", fontsize=11, color=INK, pad=9)
axs[0].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
# a 与 DUTY*
aa = np.linspace(0.02, 1.0, 300); xs = np.sqrt(1 + 1 / aa) - 1; duty_star = xs / (1 + xs)
axs[1].plot(aa, duty_star, color=INK2, lw=1.6, label="准定常最优 DUTY* = x*/(1+x*)，x* = √(1+1/a) − 1（不限舵机速率）")
axs[1].axhline(0.5, color=MUTED, lw=1, ls=(0, (3, 3))); axs[1].text(0.62, 0.51, "舵机速率上限 → 回收不能比划水快 → DUTY ≤ 0.5（同路程）", color=MUTED, fontsize=8.5)
for k in W: axs[1].plot(rec[k]["a"], np.interp(rec[k]["a"], aa, duty_star), "o", color=COL[k], ms=10, zorder=5, label=f"{k}：a = {rec[k]['a']:.2f} → DUTY* = {np.interp(rec[k]['a'], aa, duty_star):.2f}")
axs[1].set_xscale("log"); axs[1].set_xlabel("有效阻力面积比 a = (F̄_rec / F̄_pw)·(v_p / v_r)²（同一运动学的 closed / open 算例）", color=INK2); axs[1].set_ylabel("DUTY*（划水占比）", color=INK2)
axs[1].set_title("B  阻力比越小，理论上回收应越快（DUTY* 越大）\n     实际被舵机速率卡住：两相都用最大速率，DUTY 由路程比定 ≈ 0.4–0.5", loc="left", fontsize=11, color=INK, pad=9)
axs[1].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower right")
# 回收段时程
for k in W:
    d = D[k]; tau = d["tau"]; rc = tau >= d["DUTY"]
    axs[2].plot(tau[rc], smooth(d["F"]["closed"][:, 0], tau, 0.01)[rc] * 1e3, color=COL[k], lw=1.8, label=f"{k}  {rec[k]['w']:.0f} mm  回收段均值 {rec[k]['rec_mean']:.1f} mN")
axs[2].axhline(0, color=MUTED, lw=1); axs[2].set_xlabel("周期相位 τ（回收相 0.45 → 1）", color=INK2); axs[2].set_ylabel("F_x [mN]（closed 算例）", color=INK2)
axs[2].set_title("C  回收相阻力时程：整段随宽度等比例缩小\n     τ ≈ 0.6 的负峰（蜷缩 + 前扫加速）和 τ ≈ 0.85 的正峰（减速）都 ∝ 宽度", loc="left", fontsize=11, color=INK, pad=9)
axs[2].legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower right")
fig.suptitle("C6：收拢宽度的作用 —— 回收阻力与宽度成正比；12 mm 已把回收拖累压到 6 %，再收窄的绝对收益很小", fontsize=13, x=0.05, ha="left", y=0.965, color=INK)
plt.savefig(os.path.join(a.out, "fig_dd_recovery.png"), dpi=110); plt.close()
OUTJ["recovery"] = dict(points=rec, fit_slope_mN_per_mm=float(kfit), fit_intercept_mN=float(bfit)); print("recovery:", json.dumps(rec, indent=None))

# ------------------------------------------------------------------ 图 4：(τ_o, τ_c) 二维图
to_g = np.round(np.arange(-0.16, 0.121, 0.01), 3); tc_g = np.round(np.arange(0.20, 0.601, 0.01), 3)
fig, axs = plt.subplots(1, 2, figsize=(14.5, 6.2), facecolor=SURF, gridspec_kw=dict(wspace=0.22, left=0.06, right=0.97, top=0.80, bottom=0.12))
maps = {}
for ax, k in zip(axs, ("V3c12", "V3trap")):
    d = D[k]; Z = np.array([[composite(d, to, tc)["F"][:, 0].mean() * 1e3 for to in to_g] for tc in tc_g]); PZ = np.array([[composite(d, to, tc)["P"].mean() * 1e3 for to in to_g] for tc in tc_g])
    i, j = np.unravel_index(np.argmax(Z), Z.shape); style(ax)
    cf = ax.contourf(to_g, tc_g, Z, levels=24, cmap="viridis"); cs = ax.contour(to_g, tc_g, Z, levels=8, colors="white", linewidths=0.6, alpha=0.7); ax.clabel(cs, fontsize=7.5, fmt="%.0f")
    ax.plot(0, d["DUTY"], "o", color="white", ms=9, mec=INK, mew=1.2); ax.annotate(f"现用切换：开 τ=0 / 闭 τ=DUTY\n{Z[np.argmin(abs(tc_g-d['DUTY'])), np.argmin(abs(to_g))]:.1f} mN", (0, d["DUTY"]), xytext=(0.02, d["DUTY"] + 0.09), color="white", fontsize=9, arrowprops=dict(arrowstyle="->", color="white"), bbox=dict(boxstyle="round,pad=0.3", fc="#1b1b1b", ec="none", alpha=0.6))
    ax.plot(to_g[j], tc_g[i], "D", color="#ffd166", ms=10, mec=INK, mew=1.2); ax.annotate(f"最优：开 {to_g[j]:+.2f} / 闭 {tc_g[i]:.2f}（= 划水 {100*tc_g[i]/d['DUTY']:.0f} %）\n{Z[i,j]:.1f} mN（+{100*(Z[i,j]/Z[np.argmin(abs(tc_g-d['DUTY'])), np.argmin(abs(to_g))]-1):.0f} %），P̄ {PZ[i,j]:.1f} mW", (to_g[j], tc_g[i]), xytext=(-0.155, tc_g[i] - 0.12), color="#ffd166", fontsize=9, arrowprops=dict(arrowstyle="->", color="#ffd166"), bbox=dict(boxstyle="round,pad=0.3", fc="#1b1b1b", ec="none", alpha=0.75))
    ax.axhline(d["DUTY"], color="white", lw=0.8, ls=(0, (3, 3)), alpha=0.7); ax.axvline(0, color="white", lw=0.8, ls=(0, (3, 3)), alpha=0.7)
    ax.set_xlabel("开扇时刻 τ_o（负 = 回收末段提前张开）", color=INK2); ax.set_ylabel("收扇时刻 τ_c", color=INK2)
    ax.set_title(f"{NAME[k]}   T = {d['T']:.2f} s, DUTY = {d['DUTY']:.3f}", loc="left", fontsize=11, color=INK, pad=8)
    plt.colorbar(cf, ax=ax, pad=0.02, shrink=0.9).set_label("周期平均推力 F̄x [mN]", color=INK2)
    maps[k] = dict(to=to_g.tolist(), tc=tc_g.tolist(), Fx=Z.tolist(), P=PZ.tolist(), best_to=float(to_g[j]), best_tc=float(tc_g[i]), best_Fx=float(Z[i, j]), best_P=float(PZ[i, j]),
                   base_Fx=float(Z[np.argmin(abs(tc_g - d["DUTY"])), np.argmin(abs(to_g))]), best_tc_frac=float(tc_g[i] / d["DUTY"]))
fig.suptitle("开合时刻二维扫描（两态合成，第 3 周期）：最优都在「划水 82–86 % 处收、回收结束前 0.06 T 开」，且最优区很平（±0.05 内损失 < 3 %）", fontsize=13, x=0.06, ha="left", y=0.965, color=INK)
fig.text(0.06, 0.895, "限制：合成假定开/闭切换瞬时完成，未模拟折扇过渡本身（C7）；真实过渡 ≈ 60 ms ≈ 0.08 T，指令要比这里的时刻再提前约半个过渡时间", fontsize=9.5, color=MUTED)
plt.savefig(os.path.join(a.out, "fig_dd_timing_map.png"), dpi=110); plt.close()
OUTJ["timing_maps"] = {k: {kk: vv for kk, vv in v.items() if kk not in ("Fx", "P")} for k, v in maps.items()}; OUTJ["timing_maps_full"] = maps
print("timing:", json.dumps(OUTJ["timing_maps"], indent=None))

# ------------------------------------------------------------------ 图 5：瀑布
S = OUTJ["schemes"]; mt = maps["V3trap"]; m12 = maps["V3c12"]
steps = [("现步态\n(T 1.25, 60°, DUTY 0.5)", S["current"]["Fx"], MUTED, ""),
         ("V3 步态\n(T 0.81, 50°, 120→70°)", S["V3"]["Fx"], COL["V3"], "周期 ×0.65、行程居中竖直"),
         ("收拢 30 → 12 mm", S["V3c12"]["Fx"], COL["V3c12"], f"回收阻力 {S['V3']['Fx_rec']:.1f} → {S['V3c12']['Fx_rec']:.1f} mN"),
         ("梯形速度律\n(T 0.72, DUTY 0.38)", S["V3trap"]["Fx"], COL["V3trap"], f"冲量 +{100*(uns['V3trap']['I_cfd']/uns['V3']['I_cfd']-1):.0f} % × 频率 +{100*(S['V3']['T']/S['V3trap']['T']-1):.0f} %；|M_y| 峰 ×{S['V3trap']['My_peak']/S['V3']['My_peak']:.2f}"),
         (f"提前收扇 τ_c {mt['best_tc']:.2f}\n(= 划水 {100*mt['best_tc_frac']:.0f} %)", composite(D["V3trap"], 0.0, mt["best_tc"])["F"][:, 0].mean() * 1e3, "#8e1750", "减速段就收，不再拖水"),
         (f"提前开扇 τ_o {mt['best_to']:+.2f}", mt["best_Fx"], "#5e0f36", "回收减速段就开"),
         ("(可选) 羽化 5.6 mm", mt["best_Fx"] + (S["V3c06"]["Fx"] - S["V3c12"]["Fx"]), "#7fb800", f"只多 {S['V3c06']['Fx']-S['V3c12']['Fx']:.1f} mN → 不做")]
fig, ax = plt.subplots(figsize=(13, 6.6), facecolor=SURF, gridspec_kw=dict(left=0.06, right=0.98, top=0.88, bottom=0.25)); style(ax)
prev = 0
for i, (lab, val, col, note) in enumerate(steps):
    if i == 0: ax.bar(i, val, color=col, width=0.62, zorder=3)
    else:
        ax.bar(i, val - prev, bottom=prev, color=col, width=0.62, zorder=3); ax.plot([i - 1 + 0.31, i - 0.31], [prev, prev], color=INK2, lw=0.8, ls=(0, (2, 2)))
        ax.text(i, max(val, prev) + 1.5, f"{val:.1f} mN\n({'+' if val >= prev else ''}{100*(val/prev-1):.0f} %)", ha="center", va="bottom", fontsize=9.5, color=INK)
    if i == 0: ax.text(i, val + 1.5, f"{val:.1f} mN", ha="center", va="bottom", fontsize=9.5, color=INK)
    if note: ax.text(i, -3.5, note.replace("；", "\n"), ha="center", va="top", fontsize=7.8, color=MUTED)
    prev = val
ax.set_xticks(range(len(steps))); ax.set_xticklabels([s[0] for s in steps], fontsize=9.5, color=INK2); ax.tick_params(axis="x", pad=34); ax.set_ylabel("周期平均推力 F̄x [mN]（U = 0 系泊，单腿）", color=INK2); ax.set_ylim(0, prev * 1.18)
ax.set_title(f"从现步态到最优轨迹：{S['current']['Fx']:.1f} → {mt['best_Fx']:.0f} mN（× {mt['best_Fx']/S['current']['Fx']:.1f}）；相对 V3 × {mt['best_Fx']/S['V3']['Fx']:.2f}", loc="left", fontsize=12.5, color=INK, pad=10)
fig.text(0.06, 0.03, "前四项是独立 CFD 算例；后三项是两态合成的开合时刻扫描（切换假定瞬时）。代价：梯形律使力矩峰值 19 → 36 mN·m、功率 14 → 21 mW；单位功率推力 4.1 → 3.6 → 约 4.0 mN/mW（早收扇把一部分功率省回来）", fontsize=9, color=MUTED)
plt.savefig(os.path.join(a.out, "fig_dd_waterfall.png"), dpi=110); plt.close()
OUTJ["waterfall"] = [dict(label=s[0].replace("\n", " "), Fx=float(s[1]), note=s[3]) for s in steps]

# ------------------------------------------------------------------ 最优方案：V3trap + 最优开合 → 完整统计 + 事件表
d = D["V3trap"]; best = summarize("V3trap", mt["best_to"], mt["best_tc"]); OUTJ["optimal"] = best
T = d["T"]; DUTY = d["DUTY"]; ramp = d["gait"]["RAMP"]
ev = dict(T_s=T, f_Hz=1 / T, DUTY=DUTY, stroke_ms=DUTY * T * 1e3, recovery_ms=(1 - DUTY) * T * 1e3, ramp_ms=ramp * DUTY * T * 1e3, ramp_rec_ms=ramp * (1 - DUTY) * T * 1e3,
          close_cmd_tau=mt["best_tc"], close_cmd_ms=mt["best_tc"] * T * 1e3, open_cmd_tau=mt["best_to"], open_cmd_ms=(mt["best_to"] % 1.0) * T * 1e3,
          decel_start_tau=DUTY * (1 - ramp), decel_start_ms=DUTY * (1 - ramp) * T * 1e3, rec_decel_start_tau=1 - ramp * (1 - DUTY), rec_decel_start_ms=(1 - ramp * (1 - DUTY)) * T * 1e3,
          vP_max=best["vP_max"], tip_max=best["tip_max"], psi_rate_max_deg_s=float(np.degrees(np.abs(d["om"]).max())),
          sweep=d["gait"]["SWEEP"], psi0=d["gait"]["psi0_deg"], phi_ext=d["gait"]["PHI_EXT"], phi_ret=d["gait"]["PHI_RET"], switch_frac=d["gait"]["SWITCH_FRAC"], closed_w_mm=d["gait"]["CLOSED_W"] * 1e3)
OUTJ["events"] = ev
print("optimal:", json.dumps(best, indent=None)); print("events:", json.dumps(ev, indent=None))
json.dump(OUTJ, open(os.path.join(a.out, "deep_dive_optimal.json"), "w"), indent=1, ensure_ascii=False)
print("→", a.out)
