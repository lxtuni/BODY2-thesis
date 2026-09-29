# 水动力保真度核对工具 / Hydro-fidelity validation toolkit

三个工具，对应从"算一算"到"在仿真里测"到"补上缺失物理"。

## 1. `validate_hydro.py` —— 纸面参照 (一秒出结果)
```
python validate_hydro.py
```
打印：解析叶素阻力预测的终端速度、Re/Fr/St、以及 Webots 完全没建模的项（附加质量等）。
改顶部参数即可对你自己的步态/桨重算。当前结论：解析 ~0.20 m/s vs 仿真 0.015 m/s（差 ~13×），
St≈5（远高于高效区 0.2–0.4），**附加质量 ~137 g（桨自重 4.6×），峰值力比阻力推力还大**。

## 2. `worlds/force_rig.wbt` —— 测力台 (在仿真里量)
一条桨装在带**力矩反馈**的电机上、固定底座、泡在水里。控制器
`controllers/force_rig/force_rig.py` 驱动桨摆动，读电机力矩反馈，扣掉惯量项得到
**Webots 实际产生的水动力力矩**，并与解析的"阻力"和"附加质量"逐拍对比，写出 `rig_log.csv`。

用法：
1. `File ▸ Open World ▸ worlds/force_rig.wbt`，点 ▶；看控制台。
2. 控制器顶部 `MODE`：
   - `"oscillate"`（默认）：正弦摆动。会看到**实测水动力力矩 ≈ 解析阻力**，
     而**附加质量那一项是缺的**（控制台直接给出"缺失附加质量占实测的百分比"）。
   - `"steady"`：匀速旋转，α=0，纯阻力 → 反推**有效 Cd**，和你设的 1.10 对比，验证阻力定律本身。
3. 用 Excel/Python 打开 `rig_log.csv` 画 `tau_hydro_meas / tau_drag / tau_added` 曲线最直观。

> 测力原理：`tau_motor = I_hinge·α + tau_hydro` ⇒ `tau_hydro = tau_motor − I_hinge·α`。
> 桨设成中性浮力，所以重力/浮力不产生力矩，测得的就是纯水动力。

## 3. `plugins/physics/added_mass/` —— 补上附加质量 (C 物理插件)
把内置流体缺的**附加质量力** `F = −m_add·a_n`（`m_add = ρ·π·(w/2)²·L`，沿桨宽面法向）
加回来。对应你文献总结里 Han 2025 的经验公式 `F_A = 2πρa³·dv/dt`。

编译并启用：
1. 确保设了环境变量 `WEBOTS_HOME`（指向 Webots 安装目录）。
2. 在 `plugins/physics/added_mass/` 里运行 `make`（Windows 用 Webots 自带的工具链；
   通常**直接在 Webots 里打开带该插件的世界，Webots 会自动编译**）。
3. 启用方式 = 在某个世界的 `WorldInfo` 里加一行 `physics "added_mass"`：
   - 现成的：`worlds/swim_addedmass.wbt`（游泳 demo + 插件）。
   - 测力台：在 `force_rig.wbt` 的 `WorldInfo` 里加 `physics "added_mass"`，
     再跑——这时**实测水动力力矩会从"只有阻力"变成"阻力 + 附加质量"**，
     和控制器里的 `tau_drag + tau_added` 对上，直接验证缺口被补上。

插件靠 DEF 名字找桨：`DEF FL_PADDLE/FR_PADDLE/RL_PADDLE/RR_PADDLE`（游泳）和 `DEF RIG_PADDLE`（测力台），
世界里都已加好；插件会跳过当前世界里不存在的 DEF。
> 注意：附加质量(137g) 远大于桨自重(30g)，显式积分偏刚；插件里已做加速度低通+`GAIN`，
> 若发抖就把 `added_mass.c` 顶部 `GAIN` 调到 0.5、或把世界 `basicTimeStep` 降到 1。

## 参照层级 / fidelity ladder
仿真(内置流体) → 解析叶素(Li 2019) → +附加质量(Han 2025) → CFD → 拖曳水池/真机实测。
每上一层差距来源更清楚；本工具覆盖前三层，并给出和 CFD/实验对比的量纲基准(Re/Fr/St)。

## 内置流体仍然没有的项 (做研究时要心里有数)
非定常涡脱/涡环推力、环量升力、自由液面兴波与溅水、Cd 随 Re 变化、腿间流场干扰
(Qu 2025 指出真实是非线性叠加)、以及为数值稳定加的 `defaultDamping` 这种非物理阻尼。
