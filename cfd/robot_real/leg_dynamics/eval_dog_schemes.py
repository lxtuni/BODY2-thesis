import sys, os, json, numpy as np
os.chdir("/home/claude/leg_dynamics"); sys.argv = ["analyze_stroke.py", "results_leg"]
src = open("analyze_stroke.py").read(); cut = src.index("# ---------------- 扫描")
g = {"__file__": os.path.abspath("analyze_stroke.py"), "__name__": "__main__"}
import matplotlib; matplotlib.use("Agg")
exec(compile(src[:cut], "analyze_stroke.py", "exec"), g)
cand = g["candidate"]; T = g["T"]; DUTY = g["DUTY"]; SWEEP = g["SWEEP"]; PSI0 = g["PSI0"]; DPHI0 = g["DPHI0"]
J_SW0, V_SW0, T_SW0, J_back, V_BACK0, P_rec0, J_rec0 = g["J_SW0"], g["V_SW0"], g["T_SW0"], g["J_back"], g["V_BACK0"], g["P_rec0"], g["J_rec0"]
print("consts", J_SW0, V_SW0, T_SW0, J_back, V_BACK0, P_rec0, J_rec0)
def with_recovery(c, t_rec, dphi=DPHI0, U=0.0):
    """把 candidate 的回收相时长强制为 t_rec（慢蜷缩摊到整个回收相），重算损失与周期平均"""
    t_sw_eff = t_rec / 2; v_sw = V_SW0 * (dphi / DPHI0) * (T_SW0 / t_sw_eff)
    J_sw = J_SW0 * (dphi / DPHI0) ** 2 * (T_SW0 / t_sw_eff) * ((v_sw + U) / v_sw) ** 2
    t_back = t_rec; v_back = V_BACK0 * (c["A"] / SWEEP) * (DUTY * T / t_back)
    J_bk = J_back * (c["A"] / SWEEP) * (v_back / V_BACK0) * ((v_back + U) / v_back) ** 2
    J_rec = J_sw + J_bk; T_ = c["t_p"] + t_rec; P_rec = P_rec0 * T * abs(J_rec) / abs(J_rec0)
    P_pow = c["Pmean"] * c["T"] - P_rec0 * T * abs(c["J_rec"]) / abs(J_rec0)
    o = dict(c); o.update(t_rec=t_rec, T=T_, f=1 / T_, duty=c["t_p"] / T_, J_sw=J_sw, J_bk=J_bk, J_rec=J_rec, Fmean=(c["J_pow"] + J_rec) / T_, Pmean=(P_pow + P_rec) / T_)
    o["eff"] = o["Fmean"] / o["Pmean"]; return o
out = {}
cur = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "tuck"); out["current"] = cur
out["full_slowtuck"] = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck")
# A：狗式占空比 0.35，周期不变 1.25 s → 划水 0.4375 s（速率 ×1.43），回收 0.8125 s 慢蜷缩
rfA = 0.625 / 0.4375
cA = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", rf=rfA); out["dog_T1.25"] = with_recovery(cA, 0.8125)
# B：狗式占空比 0.35，划水相时长不变 0.625 s（速率不变），回收拉长到 1.161 s → T = 1.786 s
cB = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, 0.0, "rate", "slowtuck", rf=1.0); out["dog_T1.79"] = with_recovery(cB, 0.625 / 0.35 * 0.65)
# C：狗式 + 动力相在 90° 结束（136→90°，A=46°），周期不变：t_p = 46/60×0.4375? 用同一速率 ×1.43 → t_p=0.335，回收 0.915
cC = cand(46, 113, DPHI0, 0.0, "rate", "slowtuck", rf=rfA); out["dog_T1.25_end90"] = with_recovery(cC, 1.25 - cC["t_p"])
for U in (0.05, 0.1):
    cAu = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, U, "rate", "slowtuck", rf=rfA); out[f"dog_T1.25_U{U}"] = with_recovery(cAu, 0.8125, U=U)
    out[f"current_U{U}"] = cand(SWEEP, PSI0 - SWEEP / 2, DPHI0, U, "rate", "tuck")
keys = ["A", "psi_c", "rf", "t_p", "t_rec", "T", "f", "duty", "J_pow", "J_rec", "Fmean", "Fz_mean", "Pmean", "eff", "peak_Fx", "M_peak", "v_tip_max"]
for k, v in out.items():
    print(f"{k:<18}", {kk: (round(v[kk], 3) if isinstance(v[kk], float) else v[kk]) for kk in keys if kk in v})
json.dump({k: {kk: v[kk] for kk in keys if kk in v} for k, v in out.items()}, open("results_leg/dog_schemes.json", "w"), indent=1)
