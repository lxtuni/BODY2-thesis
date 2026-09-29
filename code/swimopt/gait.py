# -*- coding: utf-8 -*-
"""
③ 控制层 — 通用步态参数化, 与具体机器人无关.

自动扫描模型里的执行器(actuator), 为每个执行器生成正弦指令:
    q_i(t) = offset_i + Σ_{h=1..H} amp_i^h · sin(2π·h·f·t + phase_i^h)

H=2(双谐波)是必要的: 单谐波(H=1)时反相腿的推力精确反号→四腿抵消, 净推力恒为0;
二次谐波在相位平移 π 下不变, 正好对应真机被动脚蹼的"整流"顺桨行为.
参数向量 = [f] + 各执行器的 (amp^h, phase^h, offset)
可用 groups 把若干执行器的 amp/offset 绑在一起(共享), 降低维度:
    groups = [{"match":"1.1", "share":["amp","offset"]},
              {"match":"2.1", "share":["amp","offset"]}]
换机器人时: 执行器数量变了, 参数向量自动变长, 其余代码不用改.
"""
import numpy as np
import mujoco


class SineGait:
    def __init__(self, model, cfg=None):
        cfg = cfg or {}
        self.m = model
        self.names = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i) or f"act{i}"
                      for i in range(model.nu)]
        self.n = model.nu
        self.freq_range = cfg.get("freq_range", [0.2, 3.0])
        self.amp_range = cfg.get("amp_range", [0.0, 0.9])
        self.off_range = cfg.get("offset_range", [-0.5, 0.5])
        self.duty_range = cfg.get("duty_range", [0.25, 0.75])   # 发力相占空比(划水狗论文: 25%/33%/50%)
        self.ramp_t = cfg.get("ramp_t", 1.5)
        self.groups = cfg.get("groups", [])
        self.H = int(cfg.get("harmonics", 2))          # 谐波数
        self.amp2_scale = float(cfg.get("amp2_scale", 1.0))
        self._build_index()
        # ---- 参数冻结: 只优化其中一部分, 其余固定 ----
        self.base_x = np.full(self.dim, 0.5)           # 被冻结的参数取这里的值
        self.opt_flags = dict(freq=True, duty=True, amp=True, phase=True, offset=True)
        self.opt_flags.update({k: bool(v) for k, v in cfg.get("optimize", {}).items()})
        self.preset = cfg.get("preset_phase", "")      # 冻结相位时用哪个预设步态
        self._build_mask()

    # ---- 参数索引: 决定哪些执行器共享 amp/offset ----
    def _gid_of(self, name):
        for k, g in enumerate(self.groups):
            if g["match"].lower() in name.lower():
                return k
        return None

    def _build_index(self):
        self.amp_idx, self.off_idx = {}, {}
        slots_amp, slots_off = [], []
        for i, nm in enumerate(self.names):
            gk = self._gid_of(nm)
            share = self.groups[gk].get("share", []) if gk is not None else []
            # amp
            key = ("g", gk) if ("amp" in share) else ("i", i)
            if key not in self.amp_idx:
                self.amp_idx[key] = len(slots_amp); slots_amp.append(key)
            # offset
            key2 = ("g", gk) if ("offset" in share) else ("i", i)
            if key2 not in self.off_idx:
                self.off_idx[key2] = len(slots_off); slots_off.append(key2)
        self.n_amp = len(slots_amp)
        self.n_off = len(slots_off)
        self.n_phase = self.n            # 相位始终每执行器独立
        # 每个谐波各有一套 amp/phase; offset 只有一套
        self.dim = 2 + self.H*(self.n_amp + self.n_phase) + self.n_off   # freq+duty
        # 便于反查
        self._amp_slot = [self.amp_idx[("g", self._gid_of(nm))] if ("amp" in (self.groups[self._gid_of(nm)].get("share", []) if self._gid_of(nm) is not None else [])) else self.amp_idx[("i", i)]
                          for i, nm in enumerate(self.names)]
        self._off_slot = [self.off_idx[("g", self._gid_of(nm))] if ("offset" in (self.groups[self._gid_of(nm)].get("share", []) if self._gid_of(nm) is not None else [])) else self.off_idx[("i", i)]
                          for i, nm in enumerate(self.names)]

    # ---- 参数块索引 & 冻结 ----
    def blocks(self):
        """返回 {块名: (起, 止)} —— 便于按块冻结/显示."""
        b, k = {}, 0
        b["freq"] = (0, 1); k = 1
        b["duty"] = (k, k+1); k += 1
        for h in range(self.H):
            b[f"amp{h+1}"] = (k, k+self.n_amp); k += self.n_amp
            b[f"phase{h+1}"] = (k, k+self.n_phase); k += self.n_phase
        b["offset"] = (k, k+self.n_off)
        return b

    def preset_phases(self, name):
        """按腿名生成常见步态的相位(归一化到[0,1]). name: diag/fb/lr/wave/inphase"""
        P2 = 2*np.pi
        order = {"FL": 0, "FR": 1, "BL": 2, "BR": 3}
        table = {
            "diag":    {"FL": 0, "BR": 0, "FR": np.pi, "BL": np.pi},
            "fb":      {"FL": 0, "FR": 0, "BL": np.pi, "BR": np.pi},
            "lr":      {"FL": 0, "BL": 0, "FR": np.pi, "BR": np.pi},
            "wave":    {"FL": 0, "FR": np.pi/2, "BR": np.pi, "BL": 3*np.pi/2},
            "inphase": {"FL": 0, "FR": 0, "BL": 0, "BR": 0},
        }
        t = table.get(name, table["diag"])
        out = np.zeros(self.n_phase)
        for i, nm in enumerate(self.names):
            leg = next((L for L in order if L in nm), "FL")
            # 同一条腿内不同关节保留一个默认错相(髋0/膝60°/腕0)
            j = 1 if (".2" in nm or "2." in nm) else 0
            ph = t[leg] + (np.radians(60) if j == 1 else 0.0)
            out[i] = (ph % P2)/P2
        return out

    def _build_mask(self):
        b = self.blocks()
        m = np.zeros(self.dim, dtype=bool)
        for name, (i0, i1) in b.items():
            key = name if name in ("freq", "duty") else ("amp" if name.startswith("amp")
                  else ("phase" if name.startswith("phase") else "offset"))
            m[i0:i1] = self.opt_flags.get(key, True)
        self.mask = m
        self.dim_opt = int(m.sum())
        # 冻结相位时, 用预设步态填充 base_x
        if self.preset:
            for name, (i0, i1) in b.items():
                if name == "phase1":
                    self.base_x[i0:i1] = self.preset_phases(self.preset)
                elif name.startswith("phase"):
                    self.base_x[i0:i1] = 0.0

    def expand(self, x_red):
        """把'只优化的那几个参数'补全成完整参数向量."""
        x_red = np.asarray(x_red, float)
        if x_red.size == self.dim:
            return np.clip(x_red, 0, 1)
        full = self.base_x.copy()
        full[self.mask] = np.clip(x_red, 0, 1)
        return full

    def x0_opt(self):
        return self.base_x[self.mask].copy()

    def opt_summary(self):
        on = [k for k, v in self.opt_flags.items() if v]
        off = [k for k, v in self.opt_flags.items() if not v]
        s = f"[搜索空间] 优化 {self.dim_opt}/{self.dim} 维: {'+'.join(on)}"
        if off:
            s += f" | 冻结: {'+'.join(off)}"
            if self.preset and not self.opt_flags.get("phase", True):
                s += f" (相位固定为 '{self.preset}' 步态)"
        return s

    # ---- 优化器接口: 归一化 [0,1]^dim <-> 物理参数 ----
    def bounds(self):
        return np.zeros(self.dim), np.ones(self.dim)

    def x0(self):
        return np.full(self.dim, 0.5)

    def decode(self, x):
        x = np.clip(np.asarray(x, float), 0, 1)
        k = 0
        f = self.freq_range[0] + x[k]*(self.freq_range[1]-self.freq_range[0]); k += 1
        duty = self.duty_range[0] + x[k]*(self.duty_range[1]-self.duty_range[0]); k += 1
        A, PH = [], []
        for h in range(self.H):
            sc = 1.0 if h == 0 else self.amp2_scale
            A.append(self.amp_range[0] + x[k:k+self.n_amp]*(self.amp_range[1]-self.amp_range[0])*sc); k += self.n_amp
            PH.append(x[k:k+self.n_phase]*2*np.pi); k += self.n_phase
        offs = self.off_range[0] + x[k:k+self.n_off]*(self.off_range[1]-self.off_range[0])
        return dict(freq=float(f), duty=float(duty), A=A, PH=PH, offs=offs)

    def describe(self, x):
        p = self.decode(x)
        s = [f"freq={p['freq']:.3f} Hz  发力相占空比={p.get('duty',0.5)*100:.0f}%  (谐波 H={self.H})"]
        for i, nm in enumerate(self.names):
            harm = "  ".join(f"h{h+1}: amp={p['A'][h][self._amp_slot[i]]:+.3f} ph={np.degrees(p['PH'][h][i]):6.1f}°"
                             for h in range(self.H))
            s.append(f"  {nm:<16} off={p['offs'][self._off_slot[i]]:+.3f}  {harm}")
        return "\n".join(s)

    # ---- 仿真中调用 ----
    def ctrl(self, x_decoded, t):
        p = x_decoded
        ramp = min(1.0, t/self.ramp_t) if self.ramp_t > 0 else 1.0
        d = p.get("duty", 0.5)
        cyc = (p["freq"]*t) % 1.0                      # 周期内进度 0~1
        # 相位扭曲: 前 d 的时间走 0~π(发力半周), 其余走 π~2π(回收半周)
        th = (cyc/d)*np.pi if cyc < d else np.pi + ((cyc-d)/(1.0-d))*np.pi
        u = np.empty(self.n)
        for i in range(self.n):
            o = p["offs"][self._off_slot[i]]
            val = o
            for h in range(self.H):
                val += ramp*p["A"][h][self._amp_slot[i]]*np.sin((h+1)*th + p["PH"][h][i])
            u[i] = val
        return u

    def info(self):
        return (f"[gait] 执行器={self.n} 个, 谐波 H={self.H} → 参数维度={self.dim} "
                f"(freq 1 + H×(amp {self.n_amp} + phase {self.n_phase}) + offset {self.n_off})")
