# -*- coding: utf-8 -*-
"""
④ 优化层 — CMA-ES 步态寻优 (Hansen & Ostermeier 2001).
与机器人无关: 只调用 Swimmer.evaluate(x).

用法:  python optimize.py config.json [评估次数]
输出:  results/log.csv (每次评估) + results/best.json (最优参数) + 收敛曲线数据
"""
import sys, os, json, time, csv
import numpy as np
import cma
from simulate import Swimmer, load_cfg


def main(cfg_path, budget=None):
    cfg = load_cfg(cfg_path)
    sw = Swimmer(cfg["model"], cfg)
    print(sw.info())
    dim = sw.gait.dim_opt
    budget = int(budget or cfg.get("budget", 500))
    outdir = cfg.get("outdir", "results")
    os.makedirs(outdir, exist_ok=True)

    es = cma.CMAEvolutionStrategy(
        sw.gait.x0_opt(), cfg.get("sigma0", 0.25),
        {"bounds": [0, 1], "popsize": cfg.get("popsize", 10),
         "maxfevals": budget, "verbose": -9, "seed": cfg.get("seed", 1)})

    log = open(os.path.join(outdir, "log.csv"), "w", newline="", encoding="utf-8")
    wr = csv.writer(log); wr.writerow(["eval", "gen", "fitness", "speed", "yaw", "ok",
                                      "power", "pitch_rms", "heave_rms", "servo_rate"] +
                                      [f"x{i}" for i in range(dim)])
    n_eval, gen, t0 = 0, 0, time.time()
    best = dict(fitness=-1e9)
    hist = []
    while not es.stop() and n_eval < budget:
        X = es.ask()
        F = []
        for x in X:
            r = sw.rollout(x)
            F.append(-r["fitness"])
            n_eval += 1
            wr.writerow([n_eval, gen, f"{r['fitness']:.5f}", f"{r['speed']:.5f}",
                         f"{r['yaw']:.2f}", int(r["ok"]), f"{r['power']:.6f}",
                         f"{r['pitch_rms']:.6f}", f"{r['heave_rms']:.6f}",
                         f"{r['servo_rate']:.4f}"] + [f"{v:.4f}" for v in x])
            if r["fitness"] > best["fitness"]:
                best = dict(fitness=r["fitness"], speed=r["speed"], yaw=r["yaw"],
                            power=r["power"], pitch_rms=r["pitch_rms"],
                            heave_rms=r["heave_rms"], servo_rate=r["servo_rate"],
                            x=list(map(float, x)))
                json.dump(best, open(os.path.join(outdir, "best.json"), "w"), indent=1)
        es.tell(X, F)
        gen += 1
        hist.append(best["fitness"])
        print(f"gen {gen:3d} | evals {n_eval:4d} | best fitness {best['fitness']:+.4f} "
              f"(速度 {best['speed']:.4f} m/s, 偏航 {best['yaw']:.1f}°) | sigma {es.sigma:.4f} | {time.time()-t0:.0f}s")
    log.close()
    json.dump({"history": hist}, open(os.path.join(outdir, "convergence.json"), "w"))
    print("\n=== 最优步态 ===")
    print(sw.gait.describe(sw.gait.expand(best["x"])))
    print(f"速度 {best['speed']:.4f} m/s | 偏航 {best['yaw']:.1f}° | "
          f"共 {n_eval} 次仿真 | 耗时 {time.time()-t0:.0f}s")
    print(f"结果已存: {outdir}/best.json, log.csv, convergence.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.json",
         sys.argv[2] if len(sys.argv) > 2 else None)
