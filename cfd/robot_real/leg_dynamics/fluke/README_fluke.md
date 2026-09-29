# 升力型尾鳍单腿 CFD（+ 阻力型航速对照）

日期 2026-09-23。脚本在 `leg_dynamics\` 根目录，数据在 `leg_dynamics\fluke\`。

## 一条命令跑完

```bash
cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics
nohup bash run_fluke.sh > ~/run/run_fluke.out 2>&1 &
tail -f ~/run/run_fluke.out              # 总进度
bash progress.sh ~/run/fluke_L0          # 单个算例进度
```

只跑其中几个：`ORDER="L0 L2U000" nohup bash run_fluke.sh > ~/run/run_fluke.out 2>&1 &`。
脚本可以重复执行：跑完的跳过，被杀的从 latestTime 续算。

## 算例（默认顺序）

| 名字 | 内容 | U m/s | 次数 | 预计 |
|---|---|---|---|---|
| L0 | A2 尾鳍，运动 = Flare 单翼基准的实测关节角（相位平均），实际攻角 75 分位 19.9° | 0.223 | 1 | 3–4 h |
| L2U000 | 同 L0 的运动（固定步态，只改航速），系泊 | 0 | 1 | 3–4 h |
| D1U008 / D1U015 / D1U022 | 阻力型 V3trap（open + closed 两态），与 results_verify/V3trap 同参数 | 0.08 / 0.15 / 0.223 | 6 | 各 2.5 h |
| L1a10 / L1a28 | Flare 攻角律 α 10° / 28°（实际 α75 13° / 30°） | 0.223 | 2 | 各 3–4 h |
| L2U008 / L2U015 | 同 L0 的运动 | 0.08 / 0.15 | 2 | 各 3–4 h |
| V1 | 翼型基准 Schouveiler 2005：NACA0012 c 0.1 × 展 0.6 m，h0/c 0.75，轴 c/3，ψ 90°，St 0.25，α_max 15° | 0.4 | 1 | 8–12 h（最大，放最后） |

合计 13 次计算。

## 设置要点

- **几何（2026-09-24 改）**：默认用 `fluke/A2_blade_cfd.stl`（CFD 专用干净叶片：与打印件同平面形/截面，去掉 Ø1.9 销孔和限位 T 臂，尾缘在 95 % 弦处截平；原点 = 销轴，弦线水平，`--ref-deg 0`；俯视面积 11.2 cm²，与 Flare STL 表面偏差中位 0.12 mm）。原因：销孔和 0.1 mm 尖尾缘在网格里产生碎单元，把 Δt 压到 1e-4 s 以下。
- 旧几何 `fluke/A2_fin.stl` = flare-live `body2_fin3/meshes/LF_fin.stl`（A2 叶片 + 指节，原点在销轴，零位弦角 −25°，已用 PCA 核对）。俯视面积 11.8 cm²；C_T 用 **11.5 cm²**（A2 设计平面形）。Flare `bench_foil.py` 里 C_T 用的是 34 × 60 mm 外包矩形（20.4 cm²），所以 Flare 的 0.19 换成同一面积是 **0.34**。
- **只算尾鳍，不含脚蹼杆**。Flare 的 9.7 mN 是整腿（尾鳍 + 脚蹼杆）净推力；要严格对比，Flare 那边补一次输出"尾鳍单独的世界系受力"（`LF_fin` 的 Joint Force Torque 是尾鳍局部系下的约束力，均值不能直接用）。
- **运动**：铰链绕膝（大腿竖直 40 mm）按 7 点轨迹的时间律摆 ±24.9°，2 Hz，竖向行程 39.2 mm；尾鳍绝对角按运动表给定。L2 系列固定 L0 的运动只改 U（螺旋桨"固定转速、改进速"的做法），所以 U 越低实际攻角越大。
- **网格**：部件网格 1.2 mm（表面 0.6 mm，≈ 弦长 / 50），背景 2.5 mm，洞切体素 0.7 mm（叶片梢部最薄约 2.7 mm）。沙箱 v1912 实测：部件 17.5 万 + 背景 53.9 万单元，两级 checkMesh 通过；洞切正常（背景里约 150 个洞单元）；运动方向与运动表对照误差 0.008 mm。
- **单相、无自由液面**，顶面是离尾鳍 12 cm 的滑移刚盖（等价深水，与 Flare 条件相同）。层流。
- 力矩参考点 = 髋轴 O2；铰链力矩、功率在后处理里换算。

## 洞切检查与起步平滑（2026-09-24 加）

- **v2606 上 L0 第一次的真因**：背景网格 194824 单元 = blockMesh 原样，snappyHexMesh 没把运动区加密到 2.5 mm → 10 mm 粗单元装不下尾鳍，背景 99.9 % 被判成洞。`run_fluke.sh` 现在网格做完先查"blockMesh → snappy 加密后"的单元数，没加密就停下并打印 log.snappyHexMesh 末尾。
- 进程数默认 = 物理核数（原来 min(nproc,16) 在超线程机器上开了 16 个 = 超订）。

- `run_fluke.sh` 对尾鳍算例先只做网格（`MESH_ONLY=1 ./Allrun`），再启动求解并读第一步的 Overset analysis：**洞单元 > 3 % 总单元，或 < 50 = 洞切失败**（v2606 上 L0 第一次出现 19.5 万洞，正常应 ≈ 700），自动停掉、换下一种洞切设置重启，网格不重做。顺序：`id`（inverseDistance，体素 0.68 mm）→ `id1mm`（≈ 1.3 mm = 120³；v2606 上对这片薄叶片洞 = 0，物体内部没挖掉）→ `track2`（trackingInverseDistance，每区一组体素）→ `track` → `cvw`（cellVolumeWeight）。选中的写在算例目录 `overset_mode.txt`；失败的 log 留作 `log.overPimpleDyMFoam.fail_<模式>`。只想用某几种：`OVERSET_MODES="track2 cvw"`。
- `--ramp 0.5`：前 0.5 个周期从 t=0 位姿、速度 0 平滑起步（smoothstep），避免突然起动的压力尖峰和 Δt 塌缩；平均取最后 2 个周期，不受影响。周期数 4 → 3。
- `--limitU <m/s>`（默认关）：fvOptions limitVelocity，只在 Δt 仍被伪速度卡住时再开（整机算例用过 1.2）。

## 结果

每个算例跑完自动收集到 `results_fluke\<名字>\` 并打印汇总；随时手动汇总：

```bash
python3 fluke_forces.py summary          # → results_fluke/summary.json
```

输出：周期平均推力 T、F̄z、C_T、η、输入功率 P、T/P、铰链力矩峰值、实际攻角；阻力型给 V3trap 原最优时机（τ_c 0.33 / τ_o −0.06）和该航速下重扫的最优时机两个值；最后一张"推力–航速"对照表（升力型 L0/L2 vs 阻力型 V3trap）。

## 还没包含的（下一步再做）

- S：自由液面（尾鳍在行程顶端可能出水，需要 overInterDyMFoam，风险最高）
- P：被动尾鳍（尾鳍角换成 MuJoCo 算出的被动转角 θ(t)，用同一个 make_fluke_cfd.py）
- W：带脚蹼杆的整腿（第二个 overset 区域）
- 前后两腿串列（尾流影响，Flare 相位扫描选出 2–3 个相位后再做）

## 文件

| 文件 | 作用 |
|---|---|
| `fluke_traj.py` | 运动表生成：`--bench`（Flare 实测）/ `--alpha`（Flare 攻角律，逐行移植 run_float.py）/ `--heave-pitch`（翼型基准） |
| `make_fluke_cfd.py` | 算例生成（与 make_leg_cfd.py 同一流程；补丁仍叫 paddle，progress.sh / resume_case.sh 照用） |
| `fluke_forces.py` | 收集与汇总 |
| `run_fluke.sh` | 批量运行 |
| `fluke/A2_fin.stl`、`fluke/bench_LF_active.csv`、`fluke/motion_*.csv/.json` | 几何、Flare 基准数据、已生成的运动表 |
