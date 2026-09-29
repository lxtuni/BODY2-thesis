# Vectorised drop-in for FoilModel.apply (same equations as hydro_foil.py, checked for equality).
import numpy as np
def install(foil):
    it = foil.items
    if not it: return
    foil._gid = np.array([i["gid"] for i in it]); foil._bid = np.array([i["bid"] for i in it])
    foil._sgn = np.array([i["sgn"] for i in it]); foil._half = np.array([i["half"] for i in it])
    foil._r14 = np.array([i["r14"] for i in it]); foil._r34 = np.array([i["r34"] for i in it])
    foil._S = np.array([i["S"] for i in it]); foil._AR = np.array([i["AR"] for i in it])
    foil._CLa = np.array([i["CLa"] for i in it]); foil._V = np.array([i["V"] for i in it]); foil._ma = np.array([i["ma"] for i in it])
    foil._root = np.array([int(foil.m.body_rootid[b]) for b in foil._bid])
    def apply(data):
        f_ = foil; b = f_._bid
        R = data.xmat[b].reshape(-1, 3, 3)
        xh = f_._sgn[:, None] * R[:, :, 0]; zh = R[:, :, 2]
        pc = data.geom_xpos[f_._gid]
        Rg = data.geom_xmat[f_._gid].reshape(-1, 3, 3)
        hz = np.abs(Rg[:, 2, 0])*f_._half[:, 0] + np.abs(Rg[:, 2, 1])*f_._half[:, 1] + np.abs(Rg[:, 2, 2])*f_._half[:, 2]
        fr = np.clip((f_.water_z - (pc[:, 2] - hz)) / (2*hz), 0.0, 1.0)
        F = np.zeros((len(b), 3)); T = np.zeros((len(b), 3))
        F[:, 2] += f_._ma * f_.g
        wet = fr > 0
        if wet.any():
            cv = data.cvel[b]
            p34 = data.xpos[b] + np.einsum('nij,nj->ni', R, f_._r34)
            p14 = data.xpos[b] + np.einsum('nij,nj->ni', R, f_._r14)
            v34 = cv[:, 3:6] + np.cross(cv[:, 0:3], p34 - data.subtree_com[f_._root])
            w = -v34
            wx = np.einsum('ni,ni->n', w, xh); wz = np.einsum('ni,ni->n', w, zh)
            U = np.hypot(wx, wz)
            act = wet & (U > 1e-4)
            if act.any():
                al = np.arctan2(wz, wx)
                a = np.abs(al); sig = np.clip((1.6*f_.alpha_s - a)/(0.6*f_.alpha_s), 0.0, 1.0)
                sc = np.sin(al)*np.cos(al); CLatt = f_._CLa*sc
                CL = sig*CLatt + (1-sig)*f_.CN*sc
                CD = sig*(f_.CD0 + CLatt**2/(np.pi*f_.e*f_._AR)) + (1-sig)*f_.CN*np.sin(al)**2
                Us = np.where(U > 1e-4, U, 1.0)
                q = 0.5*f_.rho*U*U*f_._S*fr
                what = (wx[:, None]*xh + wz[:, None]*zh)/Us[:, None]
                lhat = (-wz[:, None]*xh + wx[:, None]*zh)/Us[:, None]
                Fh = q[:, None]*(CL[:, None]*lhat + CD[:, None]*what)
                Fh[~act] = 0.0
                F += Fh; T += np.cross(p14 - data.xipos[b], Fh)
                f_.last_vec = dict(alpha=al, CL=CL, CD=CD, sig=sig, U=U, Fh=Fh.copy(), act=act)
            Fb = np.zeros((len(b), 3)); Fb[:, 2] = f_.rho*f_.g*f_._V*fr
            F += Fb; T += np.cross(pc - data.xipos[b], Fb)
        np.add.at(data.xfrc_applied, b, np.hstack([F, T]))
    foil.apply_orig = foil.apply
    foil.apply = apply
