# swimopt × 升力型尾鳍腿：已完成的核心 + 交给下一个对话的任务清单

日期 2026-09-21。三个新文件放进 swimopt 根目录即可，**不改动原有文件**（原文件的改动在任务 1 里，很小）。

## 已完成（核心，已在 MuJoCo 3.13 里跑通）

| 文件 | 内容 |
|---|---|
| `hydro_foil.py` | 升力型翼面模型 `FoilModel`：¾ 弦点取相对来流 → α → 准定常 C_L/C_D（Helmbold 有限展弦比 + 失速后平板法向力）→ 力作用在 ¼ 弦点 → 被动铰力矩自然产生；附加质量走 hydro.py 同款"质量戏法"。名字含 `foil` 的 box geom 自动纳入 |
| `robots/body2_fluke.xml` | 新腿型 + A4 尾鳍：大腿固定 γ = −85°，每腿一个主动铰 `shank_*`（= 固件 q_r）+ 一个被动铰 `pitch_*`（叶簧 2.5 mN·m/rad，限位 ±40°）；船体 Morison 箱 + 三段碰撞轮廓（平底 / 舭部 / 上部）；压载块把质心压到浮心之下 |
| `config_fluke.json` | 对应配置：尾鳍从 Morison 排除、`foil` 参数、单谐波 7 维搜索空间、纵摇/垂荡罚项的权重字段 |

跑通结果（未优化，杆 40° ± 25°，对角步态）：1.0 Hz 1.6 cm/s，1.5 Hz 2.9 cm/s；尾鳍被动转角 ±31°；纵摇 RMS 0.8°。
同相位（现在真机的后腿）纵摇 RMS 10°——和视频里的点头一致，说明模型抓到了主要物理。

## 任务清单（按顺序，1–3 必做）

### 1. 把 FoilModel 接进 simulate.py（10 行）
```python
from hydro_foil import FoilModel                                   # 顶部
# Swimmer.__init__ 里 self.hydro = ... 之后：
self.foil = FoilModel(self.model, cfg.get("foil", {}), cfg["hydro"].get("rho", 1000.0), cfg["hydro"].get("water_z", 0.0))
# rollout 里两处 self.hydro.apply(d) 后面各加一行：
self.foil.apply(d)
```
注意 `self.data = mujoco.MjData(self.model)` 必须在 FoilModel 之后创建（它改了 body_mass）——原代码顺序已满足。

### 2. 代价函数加稳定性与舵机约束（simulate.py rollout，~15 行）
- 每步记录 trunk 的 pitch 和 z，循环结束后取后半段 RMS：`pitch_rms`, `heave_rms`
- `fit = speed − w_yaw·… − w_energy·power − cfg["w_pitch"]·pitch_rms − cfg["w_heave"]·heave_rms`
- 舵机角速度：`rate = max |d(ctrl)/dt|`，超过 `cfg["servo_rate_max"]`（14 rad/s）按 `−0.05·(rate/max − 1)` 罚
- 四杆可达性：对 ctrl 里 q_r 的极值调用 `leg_v3.ik()`（在 leg_dynamics/firmware_v3/），∠ACO2 < 10° 或不可达 → fitness = −1

### 3. 冒烟测试
```bash
python optimize.py config_fluke.json 60     # 60 次评估应在 5 min 内跑完，best 速度 > 2 cm/s
python view.py config_fluke.json            # 看：大腿不动、杆摆、尾鳍自己转、船不点头
```
优化输出里的 `phase1` 四个相位应收敛到对角（FL≈BR, FR≈BL, 两组差 0.5）。

### 4. 标定（把模型对到已有数据）
- **尾鳍**：用 `leg_dynamics/hydro/sculling_fsi.json`（UVLM）：同一 f / 摆幅 / k 下，比较平均推力和被动转角幅值。单桨测试可以在 xml 里把 trunk 改成 `<joint type="slide" axis="1 0 0"/>` + 大质量，读推力。误差 > 20 % 就调 `foil.CLa`、`alpha_stall_deg`
- **船体阻力**：`trunk` 的 `cd[0]` 用 CFD 船体阻力（`docs/notes/05_C6_C1结果`）在 0.1–0.3 m/s 反算
- **质量/水线**：把 `trunk_hull` + `vis_ballast` 的总质量改成真机实测（含压载），静置后 O2 高度应在 −0.011 m（水线到船顶）

### 5. 与固件对接
- 最优步态 → `scull_cp_v3.py`（`QR_MEAN`/`QR_AMP`/相位）→ 7 点 + 相位偏移，格式与 `README_scull_v3.md` 相同
- 校验舵机角 0–180、峰值角速度 < 820 °/s

### 6. 论文用的对比实验（都用同一套代价函数、同一预算）
| 组 | 模型 | 说明 |
|---|---|---|
| A | `robots/body2.xml`（扇形足，纯 Morison） | 阻力型基线 |
| B | `body2_fluke.xml`，`pitch` 刚度 1e6（尾鳍焊死） | 证明被动铰的作用 |
| C | `body2_fluke.xml` 现状 | 升力型 |
| D | C + 叶簧刚度作为优化变量（0.3–5 mN·m/rad） | 弹簧与步态协同最优 |
指标：速度、推力/功率、纵摇 RMS、垂荡 RMS。

### 7. 顺手要修的两个坑
- **`robots/toy_quad.xml` 缺 `<compiler angle="radian"/>`**：MJCF 默认按"度"解释 `range`，所以玩具模型的关节被限制在 ±1.3°，它现在的 0.056 m/s 是在关节限位上硬顶出来的。加这一行后重跑示例，README 里的数字要更新。`body2_fluke.xml` 已加。
- `hydro_foil.py` 第一步会警告 `invalid value`（mj_forward 之前 geom_xmat 全零）：在 `Swimmer.__init__` 末尾加一句 `mujoco.mj_forward(self.model, self.data)` 即可。

### 8. 可选（有时间再做）
- 尾鳍碰撞轮廓换成船体 STL 凸包（现在是三段箱体近似的舭部）
- 大腿改成可动关节（`thigh_*`），看 γ 参与摆动是否有收益
- A2（尾鳍不外移）与 A4 对比：把 `fluke_*` body 的 y 偏移 0.015 改 0
