#!/usr/bin/env python3
"""
把 open / closed 两个算例渲染好的帧按相位合成一套"两态合成"帧：τ < DUTY 取 open 帧，其余取 closed 帧（与力的合成规则一致）。
    python3 compose_verify_frames.py <open_case> <closed_case> <out_dir> [--cycles-keep 2] [--frames-name frames_side]
帧序号 ↔ 时间：paddle_video.py 按算例根目录的时间目录顺序逐帧输出（含 t=0），这里用同样的顺序对应；两算例时间步按 1e-6 容差匹配。
默认只保留最后 --cycles-keep 个周期（第一周期是起动瞬态）。输出 out_dir/frame.%04d.png（复制）。
"""
import os, sys, json, shutil, argparse
ap = argparse.ArgumentParser()
ap.add_argument("open_case"); ap.add_argument("closed_case"); ap.add_argument("out")
ap.add_argument("--cycles-keep", type=float, default=2.0); ap.add_argument("--frames-name", default="frames_side")
a = ap.parse_args()

def times_of(case):
    ts = []
    for d in os.listdir(case):
        try:
            v = float(d)
            if os.path.isdir(os.path.join(case, d)): ts.append(v)
        except ValueError: pass
    return sorted(ts)
def frames_of(case):
    fr = os.path.join(case, a.frames_name)
    fs = sorted(f for f in os.listdir(fr) if f.startswith("frame.") and f.endswith(".png")) if os.path.isdir(fr) else []
    return [os.path.join(fr, f) for f in fs]
g = json.load(open(os.path.join(a.open_case, "gait.json"))); T = g["T"]; DUTY = g["DUTY"]
to, tc = times_of(a.open_case), times_of(a.closed_case); fo, fc = frames_of(a.open_case), frames_of(a.closed_case)
if len(fo) != len(to) or len(fc) != len(tc):
    print(f"!! 帧数与时间目录数不一致：open {len(fo)} 帧 / {len(to)} 时间步，closed {len(fc)} / {len(tc)}（paddle_video.py 用了 SKIP？或算例没 reconstruct？）")
    n = min(len(fo), len(to)); to, fo = to[:n], fo[:n]; n = min(len(fc), len(tc)); tc, fc = tc[:n], fc[:n]
t_end = max(to[-1], tc[-1]); t_start = max(0.0, t_end - a.cycles_keep * T)
os.makedirs(a.out, exist_ok=True)
for f in os.listdir(a.out):
    if f.startswith("frame."): os.remove(os.path.join(a.out, f))
n = 0; picked = []
for i, t in enumerate(to):
    if t < t_start - 1e-6 or t > t_end - T / 40: continue                    # 去掉最后一个与周期起点重复的 t = 3T 帧
    tau = (t / T) % 1.0
    if tau < DUTY - 1e-9 or abs(tau - 1.0) < 1e-9: src = fo[i]; tag = "open"
    else:
        j = min(range(len(tc)), key=lambda k: abs(tc[k] - t))
        if abs(tc[j] - t) > 1e-6 * max(1.0, t) + 1e-4: src = fo[i]; tag = "open(缺closed帧)"
        else: src = fc[j]; tag = "closed"
    shutil.copyfile(src, os.path.join(a.out, f"frame.{n:04d}.png")); picked.append((t, tag)); n += 1
print(f"{n} 帧 → {a.out}   t = {picked[0][0]:.3f} … {picked[-1][0]:.3f} s（T = {T}, DUTY = {DUTY}）  open {sum(1 for _, x in picked if x == 'open')} / closed {sum(1 for _, x in picked if x == 'closed')}")
json.dump(dict(T=T, DUTY=DUTY, frames=[dict(i=i, t=t, src=tag) for i, (t, tag) in enumerate(picked)]), open(os.path.join(a.out, "frames.json"), "w"), indent=1)
