# 整机动态步态 CFD（船体 + 四腿对角步态 + 自由液面）—— 准备说明

脚本都在 `robot_real/leg_dynamics/`，船体几何 `robot_real/hull_body_only.stl`（密封身体 + 盖板，封闭，17 628 三角形）。
算例放 `~/run/`，结果自动拷到 `robot_real/leg_dynamics/results_full/`。

## 1. 模型里有什么
- **船体**固定在背景网格里（贴体 snappy，表面 2 mm），水线 z = −0.010（船体中部），来流从 +x 吹向 −x = 机器人以 U 前进（拖曳水池等效）。
- **四块脚蹼**各是一个 overset 部件，运动 = `leg_dynamics.py` 的连杆运动学（你确认过的轨迹），后腿 = 前腿机构沿 x 平移 −143.5 mm。
  对角步态：FL+RR 相位 0，FR+RL 相位 0.5（`--gait diag`，另有 sync / bound / walk）。
- **自由液面**：overInterDyMFoam（VOF），k-ω SST，两相（水 / 空气），顶面 atmosphere。
- **两态脚蹼**仍然是 open / closed 两个算例，`full_forces.py` 按每条腿自己的相位拼接（划水相取 open，其余取 closed）。
- **受力输出**：`forces_FL / FR / RL / RR / hull / all` 六个函数对象，各给 **x、y、z 三个方向**的压差力、黏性力和力矩（参考点 = 船体中点、水线高度）。
  +x 推力，+y 左，+z 上；`forces_all` 就是整机合力。

## 1b. 几何精确到什么程度
- **收拢态桨 = STEP 零件 扇形足 的精确平面轮廓与板厚**（`paddle_closed_cad.json`，从 parts/扇形足_1.stl 抽出）：桨头从 s = 50 mm 处 3 mm 宽平滑展到 75 mm 处 30 mm 宽、末端 90 mm 圆角，
  投影面积 11.7 cm²（之前的理想化轮廓 13.3 cm²）；桨头板厚 2.7 mm、偏在桨杆一侧（与 CAD 一致）。桨杆用 3 × 4.5 mm 矩形杆代替 CAD 的变截面杆
  （CAD 杆在 s ≈ 45 mm 处只有 0.5 mm 厚，网格划不出来；杆的水动力可忽略），销孔按被销轴填实处理。`--paddle model` 可切回理想化轮廓。
- **张开态 = 120° 扇形模型**（STEP 里没有张开的形状），板厚取与桨头板相同的 2.7 mm（`--thk-open` 可改），内侧按间隙裁剪（见 §2）。
- **连杆（曲柄 1、链接杆、摇臂、曲柄 2、从动杆）不作为壁面进 CFD**，但它们在水下的阻力按零件实测截面用条带法逐时刻算出（含前进速度 U 的相对速度），
  `full_forces.py` 默认把它加进每条腿的 F_x、F_z 里，并在表里单独列出"其中连杆"（`--no-linkage` 关闭）。单腿算例的 `--leg` 选项可用 CFD 核对这个估算。
- 桨面网格：部件网格 2.5 mm、桨面加密一级 1.25 mm（2.7 mm 板 ≈ 2 层单元，沙箱 checkMesh OK、两面法向各 740 面）；若你的 checkMesh 报错，`--paddle-level 2`。

## 2. 准备阶段发现的关键问题：桨与船体的间隙（见 fig_clearance.png）
overset 洞切要求运动物体与固定壁面之间至少留出约 **10 mm**（≈ 2 层部件单元 + 2 层背景单元 + 体素）；
间隙不够时前沿单元找不到供体（donor），会连锁把整个背景判成洞——沙箱里复现过，整个域 U ≡ 0。
按 CAD 位置（前腿桨中心 y = 68.8 mm）：
- 收拢桨在回收相离船侧只有 **7.9 mm**（tau ≈ 0.63，桨头经过船侧靠水线处）；
- 120° 张开扇在划水相末（tau ≈ 0.4–0.5）直接**碰到船底棱**，回收相里更是与船侧相交（真机是弹性扇，会被船体挤住）。

生成器的处理（都可以用参数改）：
1. 前腿桨中心默认外移到 **y = 75 mm**（`--y-front 0.075`，比 CAD 多 6 mm）→ 收拢桨最小间隙 13.5 mm；后腿 y = 70 mm 不动（16.4 mm）。
2. `--fan open` 时逐腿检查整周期间隙，不够就**只裁扇的内侧**（靠船体一侧）直到满足 `--min-gap`：
   前腿内侧裁到 20 mm（外侧仍 36 mm，面积 ×0.78），后腿裁到 22 mm（×0.81）；间隙检查按真实轮廓与板厚逐点计算。裁掉的主要是"回收相里扇还张着"这个 open 算例本身的假象；
   划水相里前扇本来也会蹭船底棱，所以前腿的裁剪有一部分是物理的。
3. 船体表面加密到 2 mm、部件盒靠船一侧 15 mm 条带加密到 1.25 mm，保证间隙里有足够的供体单元。
4. `gait.json` 里记录了 `gaps`（每条腿的最小间隙）和 `fan_clip`（裁剪量），报告里会引用。

**需要你判断的两件事**：真机的扇在划水末尾是否真的会碰到船体？把前腿外移 6 mm 能否接受（或者你告诉我真机的实际 y）？
想坚持 CAD 位置：`--y-front 0.0688` 会因为收拢桨间隙 7.9 mm < 10 mm 直接退出；可以用 `--min-gap 0.008 --bg 0.003 --comp 0.002`
（网格 ×2、时间 ×3）试，但不保证不断链。

## 2b. 先看一眼进入仿真的模型（几何核对，不用跑 CFD，几十秒）
```bash
python3 full_geometry_check.py ~/run/full_closed_U01 --video     # 生成算例后就能跑；open 同理
```
输出到 `results_full/`：`geom_<fan>_views.png`（t=0 三视图 + 四块桨整周期离船体的最近距离曲线，含 5.6 mm 桨厚）、
`geom_<fan>_montage.png`（周期内 6 个时刻）、`geom_<fan>.mp4`（一个周期的俯视/侧视/艉视动画）、`geom_<fan>_clearance.json`。
曲线全程 ≥ 10 mm 就是不会碰船。当前默认几何：closed（CAD 桨）全程 ≥ 17.8 mm，open（裁扇后）≥ 11.0 mm。

## 3. 规模与代价（沙箱实测，1 核串行）
- closed 算例 1.09 M 单元（背景 85 万 + 4 × 4 万部件），open 约 1.3 M；网格生成串行 ~15 min。
- 求解：洞切 8 640 个洞、供体 36 553 个单元，正常；启动第 2 步界面 Courant 冲到 4，Δt 自动缩到 0.1 ms，十几步后回到 ~1–2 ms。
- 每步 1 核 43 s → 8 核约 6–8 s；2 个周期 (2.5 s) ≈ 1500–2500 步 → **每个算例 3–5 h**（磁盘：每周期 20 帧 × ~150 MB ≈ 6 GB/算例，`--writes 10` 可减半）。
- 建议顺序：`full_closed_U01` → `full_open_U01`（先拿到 U = 0.1 的整机结果）→ 再 U = 0 和 0.2 → `--sweep` 找自航点。

## 4. 命令（一键版）
```bash
cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics
nohup bash run_full.sh 0.1 > ~/run/run_full_U01.out 2>&1 &      # 生成 → 几何核对 → 跑 closed → 跑 open → 受力 → 四视角视频 + 四合一
tail -f ~/run/run_full_U01.out                                    # 看到哪一步了
```
`run_full.sh` 可以重复执行，已完成的步骤自动跳过（断电/中断后再跑一次即可接着做）。其它航速：`bash run_full.sh 0 ` / `bash run_full.sh 0.2`，
最后 `python3 full_forces.py --sweep "results_full/full_U*.json"` 得自航点。
视频：`full_video.sh <case> quad` = 3d（透视）| side（左腿平面切面）上排，top（自由面高程）| rear（艉视）下排，2×2 合成 `full_quad_<case>.mp4`（1920 宽）；
单个视角 `bash full_video.sh <case> 3d|side|top|rear`。

## 4b. 分步命令
```bash
cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics

# (0) 只生成、看一眼间隙表和网格规模（不跑）
python3 make_full_cfd.py --case ~/run/full_closed_U01 --fan closed --U 0.1
python3 make_full_cfd.py --case ~/run/full_open_U01   --fan open   --U 0.1

# (1) 跑（每个 3–5 h；Allrun = 4 个部件网格 → 背景网格 → 合并 → setFields → 分区 → 求解 → 重构）
cd ~/run/full_closed_U01 && nohup ./Allrun > allrun.out 2>&1 &
# 网格做完后先核对（合并网格 checkMesh 应 OK；洞的数量应是几千到一万，不是几十万）：
grep -a "cells:" ~/run/full_closed_U01/log.checkMesh.merged
grep -a -m1 -A3 "calculated" ~/run/full_closed_U01/log.overPimpleDyMFoam 2>/dev/null || grep -a -m1 -A3 "calculated" ~/run/full_closed_U01/log.overInterDyMFoam
# 进度
awk '/^Time = /{t=$3} /ClockTime =/{c=$7} END{printf "t=%.3f/2.50 s  用时 %.1f 分  约剩 %.1f 分\n",t,c/60,(2.5-t)*c/t/60}' ~/run/full_closed_U01/log.overInterDyMFoam

# (2) 受力：整机 x/y/z 合力、各腿（含连杆条带法阻力）、船体、净力、俯仰力矩；两态拼接（只有一个算例也能出）
python3 full_forces.py --open ~/run/full_open_U01 --closed ~/run/full_closed_U01
# (3) 视频：quad = 3d（船体 + 四桨 + 自由面 + 水下速度切面）| side（左腿平面切面）/ top（自由面高程）| rear（艉视）四合一；单视角把 quad 换成 3d/side/top/rear
bash full_video.sh ~/run/full_closed_U01 quad
# (4) 其它航速 → 净推力–航速曲线与自航点
python3 make_full_cfd.py --case ~/run/full_closed_U00 --fan closed --U 0   # 同理 U02 / open
python3 full_forces.py --sweep "results_full/full_U*.json"
```
可选：`--no-hull`（四腿、无船体、单相 `--no-fs`，Stage 3a，约 1 h，看腿间干扰）、`--legs FL,RR`、`--cycles 3`、`--writes 10`、`--np 8`。

## 4c. 跑得太慢时（第一次跑 closed 的实测：Δt 卡在 0.15 ms、每步 16.5 s → 75 h）
原因两条：① 桨面附近空气侧一个薄片单元里 2.8 m/s 的伪速度顶住了 maxCo，Δt 缩到 0.15 ms；② overset 每步重建子网格与供体搜索（mesh.update 8.6 s / 步，
与体素数无关——沙箱验证过体素少 10 倍时间不变）。对策已写进生成器默认值（`--limitU 1.2`、`--maxCo 1.0`、`--maxAlphaCo 1.0`、p_rgh 容差放宽、`maxU` 监控）；
正在跑的算例用 `resume_case.sh` 改完从最近写出时刻续跑：
```bash
bash resume_case.sh ~/run/full_closed_U01 --tracking      # 停 → 改设置（含 trackingInverseDistance）→ 续跑 → 跑完自动重构
bash progress.sh                                          # 2 分钟后看 Δt、mesh.update；--tracking 若报错或洞数不对，去掉它再执行一次
GEN_ARGS="--overset trackingInverseDistance" nohup bash run_full.sh 0.1 > ~/run/run_full_U01.out 2>&1 &   # 之后的 open 算例带同样设置
```

## 5. 跑起来后怎么判断正常
- 求解开始处 `calculated / interpolated / hole`：hole 应为几千（桨内部 + 部件盒伸进船体的部分）；若 hole 占了大半 → 间隙/供体问题，把 log 前 200 行发我。
- `Phase-1 volume fraction` 应基本不变（≈0.765），`Min(alpha.water)` 在 −1e-6 量级即可。
- 头 20 步 Δt 会很小，之后应回到 1 ms 以上；若一直 <0.3 ms，看 `Courant Number max` 出现在哪（多半是自由面与 overset 交界的伪速度），可把 `maxAlphaCo` 改 1.0。
- v2606 若报 `solver` 关键字：`constant/dynamicMeshDict` 里改成 `motionSolver multiSolidBodyMotionSolver;`；`searchBoxDivisions` 类型错误则按 fvSchemes 里的注释改成每个 zone 一组。

## 6. 之后
跑完把 `results_full/` 里的 json/png/mp4 给我（或直接说一声），我把整机报告（含 x/y/z 三向合力、各腿贡献、船体阻力、自航点、升沉/纵摇激励）做出来。
