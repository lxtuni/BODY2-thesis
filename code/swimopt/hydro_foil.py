# -*- coding: utf-8 -*-
"""
②b 物理层 — 升力型翼面（尾鳍）的准定常水动力模型。与 hydro.py 的 Morison 模型并用：
    Morison 管机身、连杆（阻力型），本模块管名字含 "foil" 的 geom（升力型）。
    在 config 里把 "foil" 加进 hydro.exclude，避免同一个 geom 被算两次。

每个翼面（一个 box geom：x = 弦向，y = 展向，z = 厚度/法向）：
    ¾ 弦点取相对来流 w = −v_¾           （Theodorsen 准定常：¾ 弦点速度已含俯仰角速度效应）
    攻角  α = atan2(w·ẑ, w·x̂)           （x̂ 从前缘指向后缘）
    升力  L = ½ρ|w|²·S·C_L(α)·f，方向 ⊥ w  ；阻力 D = ½ρ|w|²·S·C_D(α)·f，方向 ∥ w
    C_L(α) = [σ·C_Lα + (1−σ)·C_N]·sinα·cosα   σ = 1 附着流，|α| > 1.6·α_stall 时 σ = 0（平板法向力模型）
    C_D(α) = σ·(C_D0 + C_L²/(π·e·AR)) + (1−σ)·C_N·sin²α
    C_Lα = 2π·AR/(AR+2)（Helmbold，有限展弦比）
    合力作用在 ¼ 弦点（对称翼型 ¼ 弦点力矩 ≈ 0）→ 绕被动铰的俯仰力矩自动产生
    附加质量 m_a = C_a·ρ·π(c/2)²·b，用与 hydro.py 相同的"质量戏法"加进刚体，并补偿其重力
    浮力 ρgV·f，f = 浸没比例（与 hydro.py 同一算法）

参数由 config["foil"] 给：
    {"match": "foil", "CLa": null(自动), "alpha_stall_deg": 15, "CD0": 0.02, "e": 0.9, "CN": 1.8, "Ca": 1.0}
标定：同一 f / 摆幅 / 叶簧刚度下，与 UVLM（hydro/sculling_fsi.py）的平均推力对上（±20 %），
      主要调 CLa 与 alpha_stall_deg。
"""
import numpy as np
import mujoco

GT = mujoco.mjtGeom


class FoilModel:
    def __init__(self, model, cfg, rho=1000.0, water_z=0.0):
        self.m = model
        self.rho = float(cfg.get("rho", rho))
        self.water_z = float(cfg.get("water_z", water_z))
        self.g = abs(float(model.opt.gravity[2])) or 9.81
        match = cfg.get("match", "foil").lower()
        self.alpha_s = np.radians(float(cfg.get("alpha_stall_deg", 15.0)))
        self.CD0 = float(cfg.get("CD0", 0.02))
        self.e = float(cfg.get("e", 0.9))
        self.CN = float(cfg.get("CN", 1.8))
        self.Ca = float(cfg.get("Ca", 1.0))
        CLa_cfg = cfg.get("CLa", None)
        self.items = []
        for gid in range(model.ngeom):
            name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, gid) or ""
            if match not in name.lower():
                continue
            if model.geom_type[gid] != GT.mjGEOM_BOX:
                raise ValueError(f"foil geom {name} 必须是 box (x=弦, y=展, z=厚)")
            hc, hb, ht = [float(s) for s in model.geom_size[gid][:3]]
            c, b = 2 * hc, 2 * hb
            S, AR = c * b, b / c
            xc = float(model.geom_pos[gid][0])
            sgn = 1.0 if xc >= 0 else -1.0            # 弦向单位向量 = sgn·x̂（前缘靠近 body 原点=销）
            bid = int(model.geom_bodyid[gid])
            ma = self.Ca * self.rho * np.pi * (c / 2) ** 2 * b
            self.items.append(dict(
                gid=gid, bid=bid, name=name, c=c, b=b, S=S, AR=AR, sgn=sgn,
                CLa=float(CLa_cfg) if CLa_cfg else 2 * np.pi * AR / (AR + 2.0),
                r14=np.array([xc - sgn * hc / 2, 0.0, 0.0]),      # ¼ 弦点（body 系）
                r34=np.array([xc + sgn * hc / 2, 0.0, 0.0]),      # ¾ 弦点
                half=np.array([hc, hb, ht]), V=8 * hc * hb * ht, ma=ma,
            ))
        # 质量戏法：附加质量加进刚体（各向同性近似），惯量等比放大
        for it in self.items:
            bid = it["bid"]; m0 = float(model.body_mass[bid])
            if m0 > 0:
                model.body_mass[bid] = m0 + it["ma"]
                model.body_inertia[bid] *= (m0 + it["ma"]) / m0
        self.last = {}                                 # 调试：每个翼的 α、C_L、力

    def _coeffs(self, alpha, AR, CLa):
        a = abs(alpha)
        sig = np.clip((1.6 * self.alpha_s - a) / (0.6 * self.alpha_s), 0.0, 1.0)
        sc = np.sin(alpha) * np.cos(alpha)
        CL_att = CLa * sc
        CL = sig * CL_att + (1 - sig) * self.CN * sc
        CD = sig * (self.CD0 + CL_att ** 2 / (np.pi * self.e * AR)) + (1 - sig) * self.CN * np.sin(alpha) ** 2
        return CL, CD, sig

    def apply(self, data):
        """在 hydro.apply(data) 之后调用（那里会先把 xfrc_applied 清零）。"""
        m = self.m
        for it in self.items:
            gid, bid = it["gid"], it["bid"]
            R = data.xmat[bid].reshape(3, 3)            # geom 与 body 同姿态（geom 不带旋转）
            xhat, zhat = it["sgn"] * R[:, 0], R[:, 2]
            pc = data.geom_xpos[gid]
            # 浸没比例（box 竖直半跨度，含姿态）
            Rg = data.geom_xmat[gid].reshape(3, 3)
            hz = abs(Rg[2, 0]) * it["half"][0] + abs(Rg[2, 1]) * it["half"][1] + abs(Rg[2, 2]) * it["half"][2]
            f = float(np.clip((self.water_z - (pc[2] - hz)) / (2 * hz), 0.0, 1.0))
            F = np.zeros(3); T = np.zeros(3)
            F[2] += it["ma"] * self.g                   # 附加质量的重力补偿（始终）
            if f > 0.0:
                cv = data.cvel[bid]                      # [ω, v] 绕 subtree_com
                root = int(m.body_rootid[bid])
                p34 = data.xpos[bid] + R @ it["r34"]
                p14 = data.xpos[bid] + R @ it["r14"]
                v34 = cv[3:6] + np.cross(cv[0:3], p34 - data.subtree_com[root])
                w = -v34                                 # 相对来流
                wx, wz = float(w @ xhat), float(w @ zhat)
                U = np.hypot(wx, wz)
                if U > 1e-4:
                    alpha = np.arctan2(wz, wx)
                    CL, CD, sig = self._coeffs(alpha, it["AR"], it["CLa"])
                    q = 0.5 * self.rho * U * U * it["S"] * f
                    what = (wx * xhat + wz * zhat) / U
                    lhat = (-wz * xhat + wx * zhat) / U  # ⊥ w，朝 +z 一侧
                    Fh = q * (CL * lhat + CD * what)
                    F += Fh
                    T += np.cross(p14 - data.xipos[bid], Fh)
                    self.last[it["name"]] = dict(alpha_deg=float(np.degrees(alpha)), CL=float(CL), CD=float(CD), sigma=float(sig), U=U, F=Fh.copy())
                Fb = np.array([0.0, 0.0, self.rho * self.g * it["V"] * f])
                F += Fb
                T += np.cross(pc - data.xipos[bid], Fb)
            data.xfrc_applied[bid, :3] += F
            data.xfrc_applied[bid, 3:] += T

    def summary(self):
        return " | ".join(f"[foil] {i['name']}: c={i['c']*1e3:.1f}mm b={i['b']*1e3:.0f}mm S={i['S']*1e4:.1f}cm² AR={i['AR']:.2f} "
                          f"CLa={i['CLa']:.2f}/rad m_a={i['ma']*1e3:.1f}g" for i in self.items)
