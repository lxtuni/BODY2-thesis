# 单腿非定常 CFD（第 2 级）—— 重叠网格 overset + overPimpleDyMFoam

## 思路
- 桨板 = 刚体，按 `leg_dynamics.py` 的运动学（髋摆动 ψ(t) + 平行四边形 φ(t)，φ 始终在摇臂左侧）运动；
  `make_leg_cfd.py` 把轨迹导出为 `constant/6DoF.dat`（tabulated6DoFMotion：平动 E(t)−E0，绕 y 转角 −(ψ−ψ0)），已用 moveDynamicMesh 校核方向。
- 网格：部件网格（与桨轴对齐的旋转盒，snappyHexMesh 贴体，外边界 `overset`）+ 背景网格（blockMesh + 两级区域加密）→ mergeMeshes → topoSet 分区 → setFields 写 zoneID。
- **洞切 (hole cutting) 关键设置**：默认 30³ 体素铺满整个背景域，分辩不了 5.6 mm 厚的桨板，会把整个域切成洞（U 恒为 0）。
  生成器在 `fvSchemes/oversetInterpolation` 里加了 `searchBox`（运动包络）+ `searchBoxDivisions`（体素 ≈1.2 mm）。
- 脚蹼两态：`--fan open`（全程扇形张开）与 `--fan closed`（全程收拢）各算一次；`paddle_forces.py` 按相位合成（张开相位取 open，其余取 closed）。
- 单相水（桨头在水下 2–8 cm，先不算自由液面；域顶为滑移面近似刚盖）；层流（桨头 Re≈10⁴）；`--turb kOmegaSST` 可切换。
- 输出：`postProcessing/forces/0/force.dat`（每步，力矩参考点 = 髋轴 O2），`vorticity` 场随写出帧保存。

## 跑（Ubuntu/WSL；先 `cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics`）
    # 基准步态（视频估计：T=1.25 s，摆幅 60°，划水结束松线），U=0（水箱静止）
    python3 make_leg_cfd.py --case ~/run/leg_open   --fan open   --U 0
    python3 make_leg_cfd.py --case ~/run/leg_closed --fan closed --U 0
    cd ~/run/leg_open   && ./Allrun        # 8 核约 30–45 min
    cd ~/run/leg_closed && ./Allrun
    # 后处理（两态合成 + 与准定常模型对比）
    python3 /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/paddle_forces.py --open ~/run/leg_open --closed ~/run/leg_closed
    # 视频（side = y=0 切面；3d = 压力着色桨板 + 涡量等值面）
    bash /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/paddle_video.sh ~/run/leg_open both

    # 改进步态（准定常模型推荐）：占比 0.3、峰速松线、摆幅 90°、蜷缩 60°、扇 150°
    DUTY=0.3 FAN_CLOSE=0.15 SWEEP=90 PHI_RET=60 FAN_DEG=150 CLOSED_W=0.012 \
      python3 make_leg_cfd.py --case ~/run/leg2_open --fan open --U 0
    （closed 同理；前进速度用 --U 0.1）

## 看进度
    tail -f ~/run/leg_open/log.overPimpleDyMFoam
    awk '/^Time = /{t=$3} /ClockTime =/{c=$7} END{printf "t=%.3f/3.75 s  用时 %.1f 分  约剩 %.1f 分\n",t,c/60,(3.75-t)*c/t/60}' ~/run/leg_open/log.overPimpleDyMFoam

## 网格与代价（沙箱 v1912 实测同设置）
- 部件 41.5k + 背景 220k ≈ 262k 单元，checkMesh OK（非正交 35°，偏斜 0.84）；洞单元约 650（桨内部）。
- 自适应 Δt ≈ 1.2–3 ms（maxCo 0.8），3 个周期约 1500–2000 步；8 核估 30–45 min/算例。
- 加密：`--bg 0.003 --comp 0.002`（约 2.5 倍单元数，时间 ×3）。

## 并行分区
默认 `hierarchical`（np=8 → n (2 2 2)）。scotch 分区在 closed 桨算例上触发过求解器第一步崩溃
`Attempt to cast type patch to type lduInterface`（overset 并行的已知边角情况）；若 hierarchical 也报同样错误，减到 4 块（n (2 2 1)）。

## v2606 兼容性备忘（沙箱只有 v1912，这几处若报错按此改）
1. `fvSchemes` → `searchBoxDivisions (nx ny nz)` 若报类型错误，改成每个 zone 一组：`((nx ny nz) (nx ny nz))`。
2. `constant/dynamicMeshDict` → 若报 `solver` 关键字，改为 `motionSolver solidBody;` 并把 `solidBodyCoeffs {...}` 的内容提到顶层（`cellZone`, `solidBodyMotionFunction`, `tabulated6DoFMotionCoeffs`）。
3. `fluxRequired` 已弃用但仍可读，忽略警告。
4. 层流下若 PIMPLE 发散：`nCorrectors 3→4`，`maxCo 0.8→0.5`，或 `--turb kOmegaSST`。

## 结果怎么读
- `paddle_forces.png`：最后一个周期 F_x(t)（+x 推力），open / closed / composite 三条 + 准定常两条虚线。
- 周期平均净推力：composite 是两态脚蹼的估计；open−closed 的差就是"脚蹼张开值多少"。
- 和准定常最该对比的两点：① 附加质量项（准定常 −10 mN 那一坨）CFD 里到底有多大；② 划水相峰值形状。

## 报告与视频（两个算例都跑完之后，按顺序）
    cd ~/run/leg_closed && mpirun -np 8 redistributePar -reconstruct -parallel > log.redistributePar.reconstruct 2>&1   # Allrun 若已重构可跳过
    cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics
    python3 paddle_forces.py --open ~/run/leg_open --closed ~/run/leg_closed   # 合成 + 拷 force_*.dat / gait_*.json 到 results/
    bash paddle_video.sh ~/run/leg_open both        # 速度云图 + 涡量 + 整腿骨架叠画（LEG=0 可关掉骨架），自动拷到 results/
    bash paddle_video.sh ~/run/leg_closed both
    python3 analyze_leg.py results/                 # → results/report_data.json + fig_forces/phases/cycles/paddles/cfd_side/cfd_3d.png
    node make_leg_report.js results/ .              # → results/BODY2_单腿动态水动力CFD报告.docx（要 node + `npm i docx`；没有 node 就把 results/ 交给 Claude 出 docx/pdf）
- 只有 open 时，analyze_leg.py / make_leg_report.js 会自动出"预览版"（文件名带 _预览），补上 closed 后重跑即为 v1.0。
- 视频文件：results/paddle_side_leg_open.mp4、paddle_3d_leg_open.mp4、…_leg_closed.mp4；旧的 paddle_side.mp4 / paddle_3d.mp4 是没有骨架叠画的版本，可删。
- 报告里的静态图（leg_kinematics_FL.png、stroke_direction.png）从 leg_dynamics/ 或上一级 robot_real/ 找。

## 连杆的水动力（"只算了脚蹼，连杆呢"）
- `python3 linkage_drag.py results/`：整套连杆（曲柄1、链接杆、摇臂、曲柄2、从动杆，截面按 parts/*.stl 实测）在水线以下部分的条带法阻力，
  输出 results/fig_linkage.png + linkage.json，并与 results/force_open.dat 的脚蹼 CFD 力对比；`--wl -0.005` 改水线看敏感性。analyze_leg.py 末尾会自动调用。
- `python3 make_leg_cfd.py --case ~/run/leg_open_leg --fan open --U 0 --leg`：把划水相里入水的曲柄2、从动杆段（零件实测截面，最大湿长 + 2 mm）
  作为刚体附在桨上一起算（φ 不变时它们本来就和桨绕 O2 同转，运动精确），壁面 patch `links`，forces = 桨 + 连杆总力，forcesLinks = 只有连杆段。
  paddle_forces.py 会把它拷成 results/force_links_open.dat，analyze_leg.py / 报告自动引用。网格约 +50%（部件盒扩大），时间 ×1.3。
  只对 open 有意义（closed 用于回收相，那时连杆在水面以上）；`--wl` 决定湿长。
- 不建议把每根杆做成独立 overset 部件：杆件在销孔处相互接触、部件网格相互重叠，洞切/供体查找最容易失败。
