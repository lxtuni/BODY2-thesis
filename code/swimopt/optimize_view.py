# -*- coding: utf-8 -*-
"""
边训练边看 —— CMA-ES 寻优过程可视化: 每试一组参数就在窗口里播放一次, 跑完自动换下一组.

用法:
  python optimize_view.py config.json 200            # 每组都看(默认4倍速)
  python optimize_view.py config.json 200 --every 5  # 每5组看1组(其余后台快跑)
  python optimize_view.py config.json 200 --best     # 只看"刷新纪录"的那些
  python optimize_view.py config.json 200 --rt       # 实时速度播放(慢, 看得清)
窗口: 左键拖=转视角, 滚轮=缩放, 关掉窗口=提前结束并保存当前最优
"""
import sys, os, json, time, csv
import numpy as np
import mujoco, mujoco.viewer
import cma
from simulate import Swimmer, load_cfg


def rollout_view(sw, x, viewer=None, realtime=False, render_every=8):
    """Always score with Swimmer.rollout; optionally replay the same gait in the viewer."""
    result = sw.rollout(x)
    if viewer is None or not viewer.is_running() or not result["ok"]:
        return result
    m, d = sw.model, sw.data
    p = sw.gait.decode(sw.gait.expand(x))
    sw.reset_posture(p)
    dt = m.opt.timestep
    n_settle, n_steps = int(sw.settle/dt), int(sw.T/dt)
    d.ctrl[:] = sw.gait.ctrl(p, 0.0)
    t_wall = time.time()
    for k in range(n_settle):
        sw.hydro.apply(d); sw.foil.apply(d); mujoco.mj_step(m, d)
        if k % render_every == 0:
            viewer.sync()
            if not viewer.is_running(): return None
    for k in range(n_steps):
        t = k*dt
        d.ctrl[:] = sw.gait.ctrl(p, t)
        sw.hydro.apply(d); sw.foil.apply(d); mujoco.mj_step(m, d)
        if k % render_every == 0:
            viewer.sync()
            if not viewer.is_running(): return None
            if realtime:
                lag = (n_settle+k)*dt - (time.time()-t_wall)
                if lag > 0: time.sleep(lag)
    return result


def main():
    args = sys.argv[1:]
    cfg_path = args[0] if args and args[0].endswith(".json") else "config.json"
    budget = next((int(a) for a in args if a.isdigit()), None)
    every = 1
    if "--every" in args: every = int(args[args.index("--every")+1])
    best_only = "--best" in args
    realtime = "--rt" in args

    cfg = load_cfg(cfg_path)
    sw = Swimmer(cfg["model"], cfg)
    print(sw.info())
    budget = budget or cfg.get("budget", 250)
    outdir = cfg.get("outdir", "results"); os.makedirs(outdir, exist_ok=True)
    es = cma.CMAEvolutionStrategy(sw.gait.x0_opt(), cfg.get("sigma0", 0.25),
                                  {"bounds": [0, 1], "popsize": cfg.get("popsize", 10),
                                   "maxfevals": budget, "verbose": -9, "seed": cfg.get("seed", 1)})
    log = open(os.path.join(outdir, "log.csv"), "w", newline="", encoding="utf-8")
    wr = csv.writer(log); wr.writerow(["eval", "gen", "fitness", "speed", "yaw", "ok"] +
                                      [f"x{i}" for i in range(sw.gait.dim_opt)])
    n_eval, gen, t0 = 0, 0, time.time()
    best = dict(fitness=-1e9, x=list(sw.gait.x0_opt()))
    mode = "只看刷新纪录的" if best_only else (f"每{every}组看1组" if every > 1 else "每组都看")
    print(f"\n开始寻优 (预算 {budget} 次, {mode}, {'实时' if realtime else '快进'}播放)")
    print("关掉窗口 = 提前结束并保存当前最优\n")

    with mujoco.viewer.launch_passive(sw.model, sw.data) as v:
        stop = False
        while not es.stop() and n_eval < budget and not stop:
            X = es.ask(); F = []
            for x in X:
                n_eval += 1
                show = v.is_running() and (not best_only) and (n_eval % every == 0)
                r = rollout_view(sw, x, v if show else None, realtime)
                if r is None: stop = True; break          # 窗口被关
                F.append(-r["fitness"])
                wr.writerow([n_eval, gen, f"{r['fitness']:.5f}", f"{r['speed']:.5f}",
                             f"{r['yaw']:.2f}", int(r["ok"])] + [f"{q:.4f}" for q in x])
                flag = ""
                if r["fitness"] > best["fitness"]:
                    best = dict(fitness=r["fitness"], speed=r["speed"], yaw=r["yaw"],
                                x=list(map(float, x)))
                    json.dump(best, open(os.path.join(outdir, "best.json"), "w"), indent=1)
                    flag = "  ★ 新纪录!"
                    if best_only and v.is_running():       # 只看纪录模式: 立刻回放这一组
                        rollout_view(sw, x, v, realtime)
                print(f"  #{n_eval:4d}  速度 {r['speed']:+.4f} m/s  偏航 {r['yaw']:5.1f}°  "
                      f"fitness {r['fitness']:+.4f}{flag}")
            if F and not stop:
                es.tell(X[:len(F)], F); gen += 1
                print(f"--- 第 {gen} 代 | 已评估 {n_eval} | 当前最优 {best['fitness']:+.4f} "
                      f"(速度 {best.get('speed',0):.4f} m/s) | sigma {es.sigma:.4f} | 用时 {time.time()-t0:.0f}s ---")
    log.close()
    print("\n=== 最优步态 ===")
    print(sw.gait.describe(sw.gait.expand(best["x"])))
    print(f"速度 {best.get('speed',0):.4f} m/s | 偏航 {best.get('yaw',0):.1f}° | "
          f"共 {n_eval} 次 | 耗时 {time.time()-t0:.0f}s")
    print(f"已存 {outdir}/best.json —— 双击 3_看最优结果.bat 可重播")


if __name__ == "__main__":
    main()
