# -*- coding: utf-8 -*-
"""
可视化 — 在 MuJoCo 交互窗口里播放步态, 亲眼看机器人怎么游.

用法:
  python view.py config.json                 # 播放 results/best.json 里的最优步态
  python view.py config.json results/best.json
  python view.py config.json --demo          # 播放一组手工示例步态(还没优化时用)
鼠标: 左键拖=旋转视角, 右键拖=平移, 滚轮=缩放, 空格=暂停, Esc=退出
"""
import sys, json, time
import numpy as np
import mujoco, mujoco.viewer
from simulate import Swimmer, load_cfg


def demo_params(gait):
    """手工示例: 对角步态 + 脚蹼二次谐波顺桨(能产生推力)."""
    g, P = gait, np.pi
    b = g.blocks()
    x = np.full(g.dim, 0.5)
    x[b["freq"][0]] = (1.2 - g.freq_range[0])/(g.freq_range[1]-g.freq_range[0])
    if "duty" in b: x[b["duty"][0]] = 0.5
    k = b.get("amp1", (2, 2))[0]
    # h1 幅度(按组: 髋/膝/脚蹼)
    amps1 = [0.6, 0.3, 0.15][:g.n_amp] + [0.4]*max(0, g.n_amp-3)
    for j in range(g.n_amp): x[k+j] = amps1[j]/g.amp_range[1]
    k += g.n_amp
    # h1 相位: 对角步态(按名字判断左右前后)
    kn = np.radians(60)
    for i, nm in enumerate(g.names):
        leg = "FL" if "FL" in nm else "FR" if "FR" in nm else "BL" if "BL" in nm else "BR"
        base = 0.0 if leg in ("FL", "BR") else P
        j = 0 if "1." in nm else (1 if "2." in nm else 2)
        ph = base + (kn if j == 1 else 0.0)
        x[k+i] = (ph % (2*P))/(2*P)
    k += g.n_phase
    if g.H > 1:
        amps2 = [0.0, 0.0, 0.55][:g.n_amp] + [0.0]*max(0, g.n_amp-3)
        for j in range(g.n_amp): x[k+j] = amps2[j]/g.amp_range[1]
        k += g.n_amp
        for i in range(g.n_phase): x[k+i] = (np.radians(225) % (2*P))/(2*P)
        k += g.n_phase
    x[k:k+g.n_off] = 0.5
    return np.clip(x, 0, 1)


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    cfg = load_cfg(cfg_path)
    sw = Swimmer(cfg["model"], cfg)
    print(sw.info())

    if "--demo" in sys.argv:
        x = demo_params(sw.gait); tag = "手工示例步态"
    else:
        path = next((a for a in sys.argv[2:] if a.endswith(".json")),
                    cfg.get("outdir", "results") + "/best.json")
        try:
            x = np.array(json.load(open(path))["x"]); tag = f"最优步态 ({path})"
        except Exception:
            print(f"[!] 读不到 {path}, 改用手工示例步态"); x = demo_params(sw.gait); tag = "手工示例步态"

    print(f"\n播放: {tag}")
    print(sw.gait.describe(sw.gait.expand(x)))
    r = sw.rollout(x)
    print(f"\n该步态: 速度={r['speed']:+.4f} m/s  净位移={r['dist']:+.3f} m  偏航={r['yaw']:.1f}°\n")

    m, d = sw.model, sw.data
    p = sw.gait.decode(x)
    sw.reset_posture(p)
    dt = m.opt.timestep
    settle = int(sw.settle/dt)
    with mujoco.viewer.launch_passive(m, d) as v:
        t0 = time.time(); k = 0
        while v.is_running():
            t = max(0.0, k*dt - sw.settle)
            d.ctrl[:] = sw.gait.ctrl(p, 0.0) if k < settle else sw.gait.ctrl(p, t)
            sw.hydro.apply(d)
            sw.foil.apply(d)
            mujoco.mj_step(m, d)
            if k % 10 == 0:
                v.sync()
                lag = k*dt - (time.time()-t0)
                if lag > 0: time.sleep(lag)
            k += 1
            if k*dt > 60:            # 循环播放
                sw.reset_posture(p); k = 0; t0 = time.time()


if __name__ == "__main__":
    main()
