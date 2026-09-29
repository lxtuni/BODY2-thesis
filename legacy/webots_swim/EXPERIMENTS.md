# 三个验证实验 · 具体操作步骤 / Three validation experiments — step by step

全部用 `worlds/force_rig.wbt` + `controllers/force_rig/force_rig.py`（必要时加插件）。
每步改完都要 **Ctrl+Shift+R 重载世界**。控制台结果每 1~2 秒打印一次，详细数据在控制器目录的 `rig_log.csv`。

---
## 实验 1 · 阻力标定（反推有效 Cd）
目的：不让桨动，给水恒定来流，测稳态受力 → 反推 Cd，和你设的 1.10 / 文献值对比。

1. `File ▸ Open World ▸ worlds/force_rig.wbt`。
2. 左侧场景树展开 `Fluid "water"`，把 **`streamVelocity`** 从 `0 0 0` 改成 **`0.20 0 0`**（沿 X 来流冲桨宽面）。`Ctrl+S`。
3. 打开 `controllers/force_rig/force_rig.py`，顶部设 **`MODE = "drag"`**，并让 **`V_STREAM = 0.20`** 与上一步一致。
4. `Ctrl+Shift+R` 重载 ▸ ▶。控制台读 `Cd_eff`，应≈ 你设的 `1.10`：
   `[drag] ... F=0.0xx N -> Cd_eff=1.1x (set=1.10)`
5. 把 `streamVelocity` 和 `V_STREAM` 一起改到 0.10 / 0.30 重复。**真实流体 Cd 会随 Re 略变；Webots 给的是常数——这个"Cd 不随 Re 变"本身就是一处差距。**

---
## 实验 2 · 单桨动态测力（复刻 Qu 2025，验证阻力项）
目的：桨正弦摆动，在**速度最大**瞬间（纯阻力）比较"实测水动力 vs 解析叶素阻力"。

1. `worlds/force_rig.wbt`：把 `Fluid.streamVelocity` 改回 **`0 0 0`**。
2. 控制器：**`MODE = "oscillate"`**（可调 `F` 频率、`AMP` 摆幅）。
3. 重载 ▸ ▶。读控制台 `[osc] DRAG bin (v-peak)`：
   `measured ≈ blade-element，ratio≈1` → 说明 **Webots 的阻力和叶素模型一致**（阻力项可信）。
4. 想看曲线：用 Excel/Python 打开 `rig_log.csv`，画 `tau_hydro_meas` 与 `tau_drag_analytic` 随时间。

---
## 实验 3 · 附加质量缺口（量化 Webots 缺的项）
目的：在**速度≈0、加速度最大**瞬间，阻力≈0，剩下的应是附加质量。Webots 无插件时这里≈0 → 缺口；加插件后补上。

A) **先看缺口**（无插件）
1. 仍用实验 2 的设置（`oscillate`，`streamVelocity 0 0 0`）。
2. 读控制台 `ADDED-MASS bin (v~0)`：`measured ≈ 0`，而 `analytic added-mass = 某正值`。
   → 这个差值就是 **Webots 完全没建模的附加质量**。

B) **补上缺口**（启用 C 插件）
3. 给 `worlds/force_rig.wbt` 的 `WorldInfo` 加一行：`physics "added_mass"`（`Ctrl+S`）。
4. 确保设了环境变量 `WEBOTS_HOME`。`Ctrl+Shift+R` 重载——Webots 会**自动编译** `plugins/physics/added_mass/`，
   控制台出现 `[added_mass] init: m_add = 137.x g ...`。
5. 再读 `ADDED-MASS bin (v~0)`：`measured` 现在 ≈ `analytic added-mass` → **缺口被补上**。
   阻力 bin 基本不变，说明插件只补了加速度项，没动阻力。

> 编译失败常见原因：没设 `WEBOTS_HOME`、或没有 C 编译器（Windows 装 Webots 时一般自带 MinGW）。
> 也可在该插件目录手动 `make`。改了 `added_mass.c` 后要删掉生成的 `.dll` 再重载，让它重编。

---
## 把三层连起来看
- 实验 1 → 阻力定律本身对不对（有效 Cd）。
- 实验 2 → Webots 阻力项 = 叶素阻力（√）。
- 实验 3 → Webots 缺附加质量（×）→ 插件补上（√）。
再往上：把 `rig_log.csv` 的稳态/峰值力和 **CFD** 或 **拖曳水池/真机** 实测对比，即可定标到真实水模型。
