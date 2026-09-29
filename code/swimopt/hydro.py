# -*- coding: utf-8 -*-
"""
② 物理层 — 通用水动力模型 (Morison 型), 与具体机器人无关.

对模型中每个 geom 自动计算:
    F = 浮力 ρgV·f  +  阻力 −c_d·ρ·A·v|v|·f  +  粘性 −c_v·v·f   (f = 浸没比例)
附加质量用"质量戏法": 初始化时把 m_a = C_a·ρ·V 加进刚体质量, 并施加
恒定上托力 m_a·g 抵消多出的重力 (数值绝对稳定, 各向同性近似).

系数按 geom/body 名字匹配规则指定 (见 params.yaml), 因此换机器人只需改配置.
"""
import numpy as np
import mujoco

GT = mujoco.mjtGeom


def geom_volume_areas(gtype, size):
    """返回 (体积, [Ax,Ay,Az]) —— A_i = 垂直于局部 i 轴的迎流面积."""
    gtype = int(gtype)
    if gtype == GT.mjGEOM_BOX:
        sx, sy, sz = size[:3]
        V = 8 * sx * sy * sz
        return V, np.array([4*sy*sz, 4*sx*sz, 4*sx*sy])
    if gtype == GT.mjGEOM_SPHERE:
        r = size[0]
        V = 4/3 * np.pi * r**3
        A = np.pi * r**2
        return V, np.array([A, A, A])
    if gtype == GT.mjGEOM_CAPSULE:
        r, hl = size[0], size[1]
        V = np.pi*r**2*(2*hl) + 4/3*np.pi*r**3
        Aside = 2*r*(2*hl) + np.pi*r**2
        return V, np.array([Aside, Aside, np.pi*r**2])
    if gtype == GT.mjGEOM_CYLINDER:
        r, hl = size[0], size[1]
        V = np.pi*r**2*(2*hl)
        return V, np.array([2*r*2*hl, 2*r*2*hl, np.pi*r**2])
    if gtype == GT.mjGEOM_ELLIPSOID:
        a, b, c = size[:3]
        V = 4/3*np.pi*a*b*c
        return V, np.array([np.pi*b*c, np.pi*a*c, np.pi*a*b])
    # mesh / plane / 其它: 用包围球近似
    r = max(float(size[0]), 1e-4)
    V = 4/3*np.pi*r**3
    A = np.pi*r**2
    return V, np.array([A, A, A])


def vertical_half_extent(gtype, size, R):
    """几何体在世界 z 方向的半跨度(考虑当前姿态) —— 用于算部分浸没比例."""
    gtype = int(gtype)   # numpy 整数与 pybind 枚举做 `in` 比较会失败(胶囊被误当成球), 先转 int
    if gtype == GT.mjGEOM_BOX:
        return abs(R[2,0])*size[0] + abs(R[2,1])*size[1] + abs(R[2,2])*size[2]
    if gtype == GT.mjGEOM_SPHERE:
        return size[0]
    if gtype in (GT.mjGEOM_CAPSULE, GT.mjGEOM_CYLINDER):
        return abs(R[2,2])*size[1] + size[0]
    if gtype == GT.mjGEOM_ELLIPSOID:
        return np.sqrt((R[2,0]*size[0])**2 + (R[2,1]*size[1])**2 + (R[2,2]*size[2])**2)
    return max(float(size[0]), 1e-4)


class HydroModel:
    """扫描模型中所有 geom, 按名字规则分配水动力系数, 每步施加外力."""

    def __init__(self, model, cfg):
        self.m = model
        self.rho = cfg.get("rho", 1000.0)
        self.g = abs(float(model.opt.gravity[2])) or 9.81
        self.water_z = cfg.get("water_z", 0.0)
        self.rules = cfg.get("rules", [])          # [{match, cd:[x,y,z], ca, cv}]
        self.default = cfg.get("default", {"cd": [0.3, 0.3, 0.3], "ca": 0.2, "cv": 0.02})
        self.exclude = cfg.get("exclude", [])      # 名字含这些词的 geom 不参与(如 visual/world)
        self.cfg_mesh_fill = cfg.get("mesh_fill", 0.55)   # 网格→等效长方体的填充率
        self.items = []                            # 每个浸水 geom 的预计算数据
        self._scan()
        self._apply_added_mass()
        self._vectorize()

    # ---------- 初始化 ----------
    def _gname(self, gid):
        n = mujoco.mj_id2name(self.m, mujoco.mjtObj.mjOBJ_GEOM, gid)
        return n if n else f"geom{gid}"

    def _match(self, name):
        for r in self.rules:
            if r["match"].lower() in name.lower():
                return r
        return self.default

    def _scan(self):
        for gid in range(self.m.ngeom):
            name = self._gname(gid)
            if any(x.lower() in name.lower() for x in self.exclude):
                continue
            bid = int(self.m.geom_bodyid[gid])
            if bid == 0:                        # world body 不参与
                continue
            gtype = self.m.geom_type[gid]
            if gtype in (GT.mjGEOM_PLANE, GT.mjGEOM_HFIELD):
                continue
            if gtype == GT.mjGEOM_MESH:
                V, A, half = self._mesh_box(gid)
                size_eff = half
            else:
                V, A = geom_volume_areas(gtype, self.m.geom_size[gid])
                size_eff = np.array(self.m.geom_size[gid], dtype=float)
            rule = self._match(name)
            self.items.append(dict(
                gid=gid, bid=bid, name=name, V=V, A=A,
                cd=np.array(rule.get("cd", self.default["cd"]), dtype=float),
                ca=float(rule.get("ca", self.default["ca"])),
                cv=float(rule.get("cv", self.default["cv"])),
                gtype=(int(GT.mjGEOM_BOX) if gtype == GT.mjGEOM_MESH else gtype),
                size=size_eff,
                ma=float(rule.get("ca", self.default["ca"])) * self.rho * V,
            ))

    def _mesh_box(self, gid):
        """网格几何体: 用其顶点包围盒当等效长方体(比包围球准得多).
        fill = 实体填充率(网格不会填满包围盒), 默认 0.55, 可在 config 里 hydro.mesh_fill 调."""
        mid = int(self.m.geom_dataid[gid])
        adr = int(self.m.mesh_vertadr[mid]); num = int(self.m.mesh_vertnum[mid])
        v = self.m.mesh_vert[adr:adr+num].reshape(-1, 3)
        half = np.maximum((v.max(0) - v.min(0)) / 2.0, 1e-4)
        fill = float(self.cfg_mesh_fill)
        V = 8*half[0]*half[1]*half[2]*fill
        A = np.array([4*half[1]*half[2], 4*half[0]*half[2], 4*half[0]*half[1]])*fill
        return V, A, half

    def _apply_added_mass(self):
        """质量戏法: 把附加质量加进刚体, 惯量按比例放大."""
        self.ma_body = np.zeros(self.m.nbody)
        for it in self.items:
            self.ma_body[it["bid"]] += it["ma"]
        for bid in range(1, self.m.nbody):
            ma = self.ma_body[bid]
            if ma <= 0:
                continue
            m0 = self.m.body_mass[bid]
            if m0 <= 0:
                continue
            scale = (m0 + ma) / m0
            self.m.body_mass[bid] = m0 + ma
            self.m.body_inertia[bid] *= scale

    def _vectorize(self):
        n=len(self.items)
        self.gids=np.array([i["gid"] for i in self.items],dtype=np.int32)
        self.bids=np.array([i["bid"] for i in self.items],dtype=np.int32)
        self.V=np.array([i["V"] for i in self.items])
        self.A=np.array([i["A"] for i in self.items])
        self.CD=np.array([i["cd"] for i in self.items])
        self.CV=np.array([i["cv"] for i in self.items])
        self.MA=np.array([i["ma"] for i in self.items])
        self.SZ=np.array([i["size"] for i in self.items])
        self.TY=np.array([int(i["gtype"]) for i in self.items])
        self.rootid=self.m.body_rootid[self.bids]
        self.is_box=self.TY==int(GT.mjGEOM_BOX)
        self.is_sph=self.TY==int(GT.mjGEOM_SPHERE)
        self.is_cap=(self.TY==int(GT.mjGEOM_CAPSULE))|(self.TY==int(GT.mjGEOM_CYLINDER))
        self.is_ell=self.TY==int(GT.mjGEOM_ELLIPSOID)
        self._buoy_ma=self.MA*self.g

    # ---------- 每步调用 ----------
    def apply(self, data):
        """向量化: 一次计算所有 geom 的水动力."""
        m=self.m
        data.xfrc_applied[:] = 0.0
        P = data.geom_xpos[self.gids]                       # (n,3)
        R = data.geom_xmat[self.gids].reshape(-1,3,3)       # (n,3,3)
        # 竖直半跨度(含姿态)
        hz = np.empty(len(P))
        a20,a21,a22 = np.abs(R[:,2,0]),np.abs(R[:,2,1]),np.abs(R[:,2,2])
        hz[self.is_box] = (a20*self.SZ[:,0]+a21*self.SZ[:,1]+a22*self.SZ[:,2])[self.is_box]
        hz[self.is_sph] = self.SZ[self.is_sph,0]
        hz[self.is_cap] = (a22*self.SZ[:,1]+self.SZ[:,0])[self.is_cap]
        hz[self.is_ell] = np.sqrt((R[:,2,0]*self.SZ[:,0])**2+(R[:,2,1]*self.SZ[:,1])**2+(R[:,2,2]*self.SZ[:,2])**2)[self.is_ell]
        other = ~(self.is_box|self.is_sph|self.is_cap|self.is_ell)
        hz[other] = np.maximum(self.SZ[other,0],1e-4)
        hz = np.maximum(np.nan_to_num(hz, nan=1e-6), 1e-6)
        f = np.clip((self.water_z-(P[:,2]-hz))/(2*hz),0.0,1.0)[:,None]   # (n,1)
        # geom 点速度: v = cvel_lin + omega × (p - subtree_com)
        cv6 = data.cvel[self.bids]
        off = P - data.subtree_com[self.rootid]
        Vw = cv6[:,3:6] + np.cross(cv6[:,0:3], off)
        # 各向异性二次阻力(局部系)
        Vl = np.einsum('nji,nj->ni', R, Vw)                  # R^T @ v
        Fl = -0.5*self.CD*self.rho*self.A*Vl*np.abs(Vl)*f
        F = np.einsum('nij,nj->ni', R, Fl)
        F -= self.CV[:,None]*Vw*f                            # 线性粘性
        F[:,2] += (self.rho*self.g*self.V)[None,:].ravel()*f.ravel()   # 浮力
        F[:,2] += self._buoy_ma                              # 附加质量重力补偿
        T = np.cross(P - data.xipos[self.bids], F)
        np.add.at(data.xfrc_applied, self.bids, np.hstack([F,T]))
        return

    def _apply_slow(self, data):
        data.xfrc_applied[:] = 0.0
        res = np.zeros(6)
        for it in self.items:
            gid, bid = it["gid"], it["bid"]
            p = data.geom_xpos[gid]
            R = data.geom_xmat[gid].reshape(3, 3)
            # 浸没比例 f: 用几何体在竖直方向的真实跨度(含姿态), 线性近似部分浸没
            hz = vertical_half_extent(it["gtype"], it["size"], R)
            f = np.clip((self.water_z - (p[2] - hz)) / (2 * hz), 0.0, 1.0) if hz > 1e-9 else 0.0
            F = np.zeros(3)
            # 附加质量的重力补偿(始终施加, 保证重量正确)
            F[2] += it["ma"] * self.g
            if f > 0.0:
                mujoco.mj_objectVelocity(self.m, data, mujoco.mjtObj.mjOBJ_GEOM, gid, res, 0)
                v_world = res[3:6]
                v_loc = R.T @ v_world
                # 各向异性二次阻力 (Webots 无½因子; 此处显式写½)
                f_loc = -0.5 * it["cd"] * self.rho * it["A"] * v_loc * np.abs(v_loc) * f
                F += R @ f_loc
                # 线性粘性
                F += -it["cv"] * v_world * f
                # 浮力
                F[2] += self.rho * self.g * it["V"] * f
            data.xfrc_applied[bid, :3] += F
            data.xfrc_applied[bid, 3:] += np.cross(p - data.xipos[bid], F)

    # ---------- 诊断 ----------
    def summary(self):
        tot_V = sum(i["V"] for i in self.items)
        tot_ma = sum(i["ma"] for i in self.items)
        tot_m = float(np.sum(self.m.body_mass[1:]))
        return (f"[hydro] geoms={len(self.items)} 总体积={tot_V*1e6:.1f}cm³ "
                f"满浸浮力={self.rho*tot_V*1000:.1f}g水 | 附加质量={tot_ma*1000:.1f}g "
                f"| 含附加质量后总质量={tot_m*1000:.1f}g")
