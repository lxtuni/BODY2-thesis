# robot0526 四足机器人水动力仿真 (Webots 内置流体)

Hydrodynamic simulation of the robot0526 quadruped in Webots using the **built-in fluid model**
(`Fluid` node + per-link `ImmersionProperties`).

本项目把你上传的 `robot0526.urdf`（28 自由度、约 125 g、~13 cm）落到 Webots 里做**划水（dog-paddle）水动力仿真**。
采用 Webots **内置流体**方案：浮力由 Webots 自动计算，阻力用阻力系数 (Cd) 设置——无需自己写物理插件即可看到机器人浮在水面、划水前进。

---

## 1. 一个必须先知道的前提 / Key constraint

> **Webots 内置浸没（浮力）只支持基本几何体**：`Box / Cylinder / Capsule / Sphere`，
> **不支持 STL 网格 (IndexedFaceSet/Mesh)**。

你的 URDF 几何是 SolidWorks 导出的 STL，所以**不能直接对原网格加水动力**。
因此本项目用**基本体 (Box) 近似**了浸水连杆：

* 用 `urdf2webots` 把 URDF 转成 `protos/Robot0526.proto`（保留真实运动学/外形，仅供参考与陆地用）；
* 另外手写了一个**自包含、可直接运行**的演示世界 `worlds/swim_demo.wbt`，
  其中机身=长方体、每条腿简化为 **2 自由度桨**（髋摆 + 腕部顺桨），尺寸取自对你 STL 的实测。

> 想用你**精确的 28-DOF 模型**做水动力？见 [第 6 节](#6-用你的真实模型-robot0526proto)。

---

## 2. 运行步骤 / How to run

1. 安装 **Webots R2023b 或更新版**（本项目用 R2025a 头，向下兼容；官网免费：https://cyberbotics.com）。
2. 打开 Webots → `File ▸ Open World…` → 选择本项目的 `worlds/swim_demo.wbt`。
3. 点工具栏 **▶ 播放**。机器人会浮在蓝色水体表面、四条桨交替划水**向 +X 方向前进**。
4. 控制台（Console）每秒打印一次前进速度，例如：
   ```
   [paddle_swim] t=  5.0s  x=+0.041m  v=+0.030 m/s  yaw=+1.2 deg
   ```

> Webots 会在第一次运行时自动编译 Python 控制器，无需手动配置。
> 控制器用的是标准 `controller` 模块，确保 Webots 的 Python 已配好（菜单 `Tools ▸ Preferences ▸ Python command`）。

---

## 3. 水动力是怎么算的 / How the hydrodynamics work

Webots 的内置流体把力拆成两部分，**都由仿真器自动施加**，控制器只负责驱动关节：

**(a) 浮力 Buoyancy（自动）** —— 来自 `Fluid` 节点。对每个浸入水中的 `boundingObject`，按阿基米德原理：

```
F_buoy = ρ_fluid · g · V_immersed        (方向向上)
```

`V_immersed` 是该基本体浸没部分的体积，Webots 实时计算。**机器人浮还是沉，取决于「boundingObject 体积×水密度」与「连杆质量」的对比**（见第 5 节调浮力）。

**(b) 阻力 Drag（靠 `ImmersionProperties` 的 Cd）** —— 沿连杆局部坐标轴：

```
F_drag,i ≈ -½ · ρ_fluid · Cd_i · A_i · v_i · |v_i|   −   c_visc · v_i
```

* `Cd_i` = `dragForceCoefficients`（x,y,z 各一个）；
* `A_i` = 该基本体**垂直于第 i 轴**的浸没横截面积（Webots 按几何自动算）；
* `c_visc` = `viscousResistanceForceCoefficient`（低速线性阻尼）。

这正是你那份文献总结里 Li 2019 / Qu 2025 的**叶素-阻力**思想 `F=½ρCd·A·v²` 的引擎内置版。

**桨的设计**：桨叶 `Box (0.006 × 0.050 × 0.070)` 在局部 **X 方向最薄**，故垂直于 X 的面积最大、`Cd_x=1.10`（≈平板，对应 Li 2019 的 1.17）。
划水**发力相**桨叶宽面正对水流→大阻力→推力；**回收相**腕关节把桨**顺桨 (feather)** 转 ~75°→迎水面积骤减→低阻力。净推力向前。

---

## 4. 文件结构 / Project layout

```
webots_swim/
├─ worlds/
│   └─ swim_demo.wbt              ← 主场景：水池 + 机器人（直接打开这个）
├─ controllers/
│   └─ paddle_swim/
│       └─ paddle_swim.py         ← 划水步态控制器（含参数）
├─ protos/
│   ├─ Robot0526.proto           ← 你的真实 28-DOF 模型（urdf2webots 转换，供精确几何使用）
│   └─ robot0526/meshes/*.STL     ← 原始网格（Robot0526.proto 引用）
└─ README.md
```

---

## 5. 调参 / Tuning（最常用）

所有水动力参数都在 `worlds/swim_demo.wbt` 里，步态参数在 `controllers/paddle_swim/paddle_swim.py` 顶部。

| 目标 | 改哪里 | 怎么改 |
|---|---|---|
| **游得更快** | 控制器 `FREQ`、`HIP_AMP` | 调高频率(0.9→1.3 Hz)、加大摆幅 |
| **桨推力更大** | 世界文件 桨叶 `dragForceCoefficients 1.10 …` | 增大第一个数 (X 向 Cd)；或加大桨叶 Y/Z 尺寸 |
| **机器人下沉** | 减小连杆 `mass`，或增大 `boundingObject` 尺寸 | 浮力∝体积，重力∝质量 |
| **浮太高/不稳** | 增大机身 `mass 0.180`，或减小机身 `Box` 高度 | 让吃水更深 |
| **顺桨不彻底（回收阻力大）** | 控制器 `FEATHER` | 调大到 ~1.5 rad |
| **模拟水流** | `Fluid { streamVelocity 0.1 0 0 }` | 给水一个流速，测抗流能力 |
| **更稳但慢 / 更快但晃** | 控制器 `PHASE` | 对角(默认 TLPG)更稳；改成顺序相位(LSPG)更快——见 Qu 2025 |

> **浮力速算**：机身 Box `0.22×0.10×0.05 = 5.5e-4 m³`，满浸浮力 `≈0.55 kg`；
> 全机质量 `≈0.29 kg` → 机身约 1/2 吃水即平衡。想让它坐得更低就加机身质量或减小 Box 高度。

### 文献 Cd 参考值（来自你那份总结）

| 形状 | Cd | 用途 |
|---|---|---|
| 平板正对水流 (flat plate) | **1.1–1.17** | 桨叶发力面 (Li 2019: 1.17) |
| 半圆柱面 | 0.42 | 桨叶/腿恢复面 (Li 2019) |
| 长方体机身 | ~1.05 | 机身 (Li 2019 / Qu 2025) |

---

## 6. 用你的真实模型 Robot0526.proto

`protos/Robot0526.proto` 是用官方 `urdf2webots` 从你的 URDF 转出来的**精确 28-DOF 模型**，运动学、质量、外形都对。
但如上所述，它的 `boundingObject` 是 **STL 网格 → 不能产生浮力/浸没阻力**。要让它在水里游，需要两步：

1. **把每个浸水连杆的 `boundingObject` 从 Mesh 换成基本体**（Box/Cylinder），尺寸照搬本 README 第 2 节的实测值（已在 `swim_demo.wbt` 里给好范例）。
2. **给每个浸水连杆加 `immersionProperties [ ImmersionProperties { fluidName "water" … } ]`**，Cd 用上表。

> 因为该模型每条腿是 7 段细趾、关节轴是斜向单位向量，直接逐段加浸没既繁琐、又远超「看它游起来」的需求；
> 所以推荐：**先用 `swim_demo.wbt` 跑通并调好水动力参数，再按同样写法把基本体浸没属性搬到 Robot0526.proto 的对应连杆**。
> 也可只对「发力的远端桨/蹼连杆 (L_3 / L_4 / L_6)」加浸没，近端细杆忽略，工作量更小。

把 `Robot0526` 实例化进世界（替换 `swim_demo.wbt` 里的内联 `Robot{}`）的写法：

```vrml
Robot0526 {
  translation 0 0 0.02
  controller "paddle_swim"   # 注意：真实模型电机名不同，需相应改控制器里的电机名
}
```
（真实模型的电机名形如 `joint_rl`, `JFL` …；可在控制器里 `print(robot.getDevice(name))` 逐一确认。）

---

## 7. 更高精度（可选）：自定义物理插件

内置流体是**稳态阻力 + 自动浮力**，不含附加质量 (added mass)、不含多腿流场干扰。
若要复刻论文级精度（叶素积分、附加质量、非定常项），可写 **Webots physics plugin**（C/ODE，用 `dBodyAddForceAtPos` 对每段连杆按 `F=½ρCd·A·v²` 逐时步加力）。
这对应你那份《四足机器人水模型文献总结》里 Li 2019 / Qu 2025（叶素-阻力）与 Han 2025（附加质量项 `F_A=2πρa³·dv/dt`）的做法。需要的话我可以再给一份插件模板。

---

## 8. 常见问题 / Troubleshooting

| 现象 | 原因 / 解决 |
|---|---|
| 打不开世界、提示找不到 PROTO | `swim_demo.wbt` 是自包含的，不依赖外部 PROTO；确认 Webots ≥ R2023b |
| 机器人直接沉底 | 浮力不足：减小连杆质量或加大 boundingObject；确认连杆在水面 (z≈0) 以下 |
| 机器人乱飞/爆炸 | 时间步过大或初始穿模：把 `WorldInfo.basicTimeStep` 调到 4–8；机器人起始 z 略高于水面 |
| 划水但不前进 | 桨叶 Cd 各向同性了：确认桨叶 `dragForceCoefficients` 第一个数(X)远大于另两个；或顺桨方向反了，把 `FEATHER` 取反号 |
| 控制器不动 | Console 看报错；确认 `controller "paddle_swim"` 且 Python 路径已配 |

---

*生成日期 2026-06-13 · 基于你上传的 robot0526.urdf 与《四足机器人水模型文献总结》。*
