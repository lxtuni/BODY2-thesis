import json, os
HERE="/home/claude/leg_dynamics"; OUT=os.path.join(HERE,"results_verify")
K=["current","V1","V2","V3"]; LAB={"current":"现步态（基准）","V1":"V1 前沿点","V2":"V2 狗式占空比","V3":"V3 低中心角"}
G={k:json.load(open(os.path.join(OUT,f"gait_params_{k}.json"))) for k in K}
S=json.load(open(os.path.join(OUT,"verify_analysis.json")))
MODEL={"current":dict(F0=21.0,Fz=12.1,M=9.0,P=7.3),"V1":dict(F0=47.7,Fz=23.3,M=18.5,P=17.3),"V2":dict(F0=36.5,Fz=16.7,M=18.5,P=11.6),"V3":dict(F0=42.6,Fz=13.0,M=17.9,P=15.4)}
CLR={"current":(14.4,17.8),"V1":(13.5,18.0),"V2":(14.5,18.0),"V3":(9.0,17.9)}
RUN={k:{f:json.load(open(os.path.join(OUT,k,f"run_{f}.json"))) for f in ("open","closed")} for k in ("V1","V2","V3")}
RUN["current"]=None
mn=lambda s:s.replace("-","−")
def f(x,n=1): return mn(f"{x:.{n}f}")
rows=[]
def R(name,fn,unit=""):
    rows.append([name+(f" [{unit}]" if unit else "")]+[fn(k) for k in K])
def sec(t): rows.append(["**"+t+"**"]+[""]*len(K))

sec("步态输入参数（leg_dynamics.py 环境变量）")
R("T　周期",lambda k:f(G[k]["T"],2),"s")
R("DUTY　划水占空比",lambda k:f(G[k]["DUTY"],2))
R("SWEEP　髋摆幅",lambda k:f(G[k]["SWEEP"],0),"°")
R("PSI_START　划水起点摇臂角",lambda k:f(G[k]["PSI_START"],0),"°")
R("PHI_EXT　伸展张角",lambda k:f(G[k]["PHI_EXT"],0),"°")
R("PHI_RET　蜷缩张角",lambda k:f(G[k]["PHI_RET"],0),"°")
R("SWITCH_FRAC　蜷缩/伸展占回收比",lambda k:f(G[k]["SWITCH_FRAC"],2))
R("FAN_DEG　脚蹼展开角",lambda k:f(G[k]["FAN_DEG"],0),"°")
R("FAN_CLOSE　松线相位",lambda k:"−1（= DUTY）")
R("CLOSED_W　收拢桨头宽",lambda k:f(G[k]["CLOSED_W_mm"],0),"mm")

sec("派生时序与角度")
R("f　频率",lambda k:f(G[k]["f"],2),"Hz")
R("ψ 范围（起 → 止）",lambda k:f"{G[k]['PSI_START']:.0f} → {G[k]['PSI_END']:.0f}°")
R("ψ 中心角",lambda k:f(G[k]["PSI_CENTER"],0),"°")
R("划水时长 DUTY·T",lambda k:f(G[k]["t_power"],3),"s")
R("回收时长 (1−DUTY)·T",lambda k:f(G[k]["t_recovery"],3),"s")
R("　└ 蜷缩 / 保持 / 伸展",lambda k:f"{G[k]['t_tuck']:.3f} / {G[k]['t_hold_tucked']:.3f} / {G[k]['t_extend']:.3f}","s")

sec("运动学（模型，U = 0）")
R("髋角速度 划水峰值",lambda k:f(G[k]["w_power_max_deg_s"],0),"°/s")
R("髋角速度 划水平均",lambda k:f(G[k]["w_power_mean_deg_s"],0),"°/s")
R("髋角速度 回收峰值",lambda k:f(G[k]["w_rec_max_deg_s"],0),"°/s")
R("桨尖速度 峰值",lambda k:f(G[k]["vtip_max"],2),"m/s")
R("桨尖速度 划水相平均",lambda k:f(G[k]["vtip_power_mean"],2),"m/s")
R("划水相桨尖水下深度",lambda k:f"{G[k]['tip_depth_power_mm'][0]:.0f}–{G[k]['tip_depth_power_mm'][1]:.0f}","mm")
R("回收相 E 点抬升",lambda k:f(G[k]["E_lift_mm"],0),"mm")
R("桨顶最高 z（水线 −10）",lambda k:f({"current":49.6,"V1":43.1,"V2":42.6,"V3":47.3}[k],1),"mm")
R("与船体最小间隙 划水 / 回收",lambda k:f"{CLR[k][0]:.1f} / {CLR[k][1]:.1f}","mm")
R("servo_2（曲柄1）指令范围",lambda k:f"{G[k]['servo2_crank1_deg'][0]:.0f}…{G[k]['servo2_crank1_deg'][1]:.0f}°".replace("-","−"))
R("servo_1（曲柄2）指令范围",lambda k:f"{G[k]['servo1_crank2_deg'][0]:.0f}…{G[k]['servo1_crank2_deg'][1]:.0f}°")

sec("CFD 结果（overPimpleDyMFoam，U = 0，两态合成，第 3 周期）")
R("F̄_x　周期平均推力",lambda k:"**"+f(S[k]["Fx_mean"])+"**","mN")
R("　└ 行程模型预测 / 偏差",lambda k:f"{MODEL[k]['F0']:.1f} / {100*(S[k]['Fx_mean']/MODEL[k]['F0']-1):+.0f} %".replace("-","−"),"mN")
R("　└ 相对现步态",lambda k:"—" if k=="current" else f"{100*(S[k]['Fx_mean']/S['current']['Fx_mean']-1):+.0f} %")
R("F̄_z　周期平均垂向力",lambda k:f(S[k]["Fz_mean"]),"mN")
R("F̄_z / F̄_x",lambda k:f(S[k]["Fz_mean"]/S[k]["Fx_mean"],2))
R("F_x 峰 / 谷",lambda k:f"{S[k]['Fx_peak']:.0f} / {S[k]['Fx_min']:.0f}".replace("-","−"),"mN")
R("F_z 峰 / 谷",lambda k:f"{S[k]['Fz_peak']:.0f} / {S[k]['Fz_min']:.0f}".replace("-","−"),"mN")
R("划水相平均力矢量 大小 / 倾角",lambda k:f"{S[k]['Fmag_pw']:.0f} mN / {S[k]['tilt_pw_deg']:.0f}°")
R("划水相 / 回收相平均 F_x",lambda k:f"{S[k]['Fx_pw']:.0f} / {S[k]['Fx_rec']:.0f}".replace("-","−"),"mN")
R("|M_y| 峰值（髋轴）",lambda k:f(S[k]["My_peak"]),"mN·m")
R("P̄　水动力功率",lambda k:f(S[k]["P_mean"]),"mW")
R("P 峰值",lambda k:f(S[k]["P_peak"],0),"mW")
R("F̄_x / P̄　效率",lambda k:f(S[k]["Fx_mean"]/S[k]["P_mean"],2),"mN/mW")
R("两对腿叠加 F_z 峰峰（升沉激励）",lambda k:f(S[k]["Fz_pair_pp"],0),"mN")
R("逐周期 F̄_x：1 / 2 / 3",lambda k:" / ".join(f"{p['Fx']:.1f}" for p in S[k]["per_cycle"]),"mN")
R("　└ 2→3 变化",lambda k:(lambda c:f"{100*(c[-1]['Fx']/c[-2]['Fx']-1):+.1f} %".replace("-","−"))(S[k]["per_cycle"]))

sec("算例开销（8 核，每方案 open + closed）")
R("时间步数 open / closed",lambda k:"—" if RUN[k] is None else f"{RUN[k]['open']['steps']} / {RUN[k]['closed']['steps']}")
R("机时 open / closed",lambda k:"—" if RUN[k] is None else f"{RUN[k]['open']['clock_s']/3600:.2f} / {RUN[k]['closed']['clock_s']/3600:.2f}","h")

hdr=["参数"]+[LAB[k] for k in K]
esc=lambda x: x.replace("|","\\|") if x.startswith("|M_y|") else x
md=["| "+" | ".join(hdr)+" |","|"+"---|"*(len(K)+1)]+["| "+" | ".join([esc(r[0])]+r[1:])+" |" for r in rows]
txt="# BODY2 单腿步态参数表（现步态 + V1 / V2 / V3）\n\n"+"\n".join(md)+"\n\n复现命令（在 leg_dynamics/）：\n\n```\n"
for k in K[1:]:
    g=G[k]; txt+=f"# {LAB[k]}\nT={g['T']} DUTY={g['DUTY']} SWEEP={g['SWEEP']:.0f} PSI_START={g['PSI_START']:.0f} SWITCH_FRAC={g['SWITCH_FRAC']} python3 make_leg_cfd.py --case ~/run/leg_{k}_open --fan open --U 0 --cycles 3 --np 8\n"
txt+="```\n"
open(os.path.join(OUT,"gait_params_table.md"),"w").write(txt)
json.dump(dict(header=hdr, rows=rows, cmd=[dict(name=LAB[k], env=f"T={G[k]['T']} DUTY={G[k]['DUTY']} SWEEP={G[k]['SWEEP']:.0f} PSI_START={G[k]['PSI_START']:.0f} SWITCH_FRAC={G[k]['SWITCH_FRAC']}") for k in K[1:]]),
          open(os.path.join(OUT,"gait_params_table.json"),"w"), indent=1, ensure_ascii=False)
print(txt)
