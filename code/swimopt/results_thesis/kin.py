# Lift/drag decomposition of flipper thrust for Morison-model runs.
# For every immersed flipper sample: world force F (anisotropic drag of plate),
# relative-flow direction w = -v. Drag part = (F.ŵ)ŵ, lift part = F - drag part.
# Thrust shares = cycle-averaged x-components of each part.
from common import *
import mujoco
def decompose(sw, x, t_from=None, match="flip"):
    g=sw.gait; m,d=sw.model,sw.data
    mujoco.mj_resetData(m,d); mujoco.mj_forward(m,d)
    p=g.decode(g.expand(x)); dt=m.opt.timestep
    for _ in range(int(sw.settle/dt)): sw.hydro.apply(d); mujoco.mj_step(m,d)
    idx=[i for i,it in enumerate(sw.hydro.items) if it['name'].startswith(match)]
    t_from = sw.T/2 if t_from is None else t_from
    Dx=Lx=0.0; alph=[]; wts=[]; n=0; vz=vx=0.0
    for k in range(int(sw.T/dt)):
        t=k*dt; d.ctrl[:]=g.ctrl(p,t); sw.hydro.apply(d)
        if t>=t_from:
            for i in idx:
                it=sw.hydro.items[i]; gid=it['gid']
                R=d.geom_xmat[gid].reshape(3,3)
                hz=abs(R[2,0])*it['size'][0]+abs(R[2,1])*it['size'][1]+abs(R[2,2])*it['size'][2]
                f=np.clip((0-(d.geom_xpos[gid][2]-hz))/(2*hz),0,1)
                res=np.zeros(6); mujoco.mj_objectVelocity(m,d,mujoco.mjtObj.mjOBJ_GEOM,gid,res,0)
                v=res[3:6]; sp=np.linalg.norm(v)
                if sp<1e-4 or f<=0: continue
                vl=R.T@v; F=R@(-0.5*it['cd']*sw.hydro.rho*it['A']*vl*np.abs(vl)*f)
                wh=-v/sp; Fd=(F@wh)*wh; Fl=F-Fd
                Dx+=Fd[0]*dt; Lx+=Fl[0]*dt; n+=1
                a=np.degrees(np.arcsin(min(1,abs(vl[0])/sp)))
                alph.append(a); wts.append(max(F[0],0))
                vx+=abs(v[0])*dt; vz+=np.hypot(v[1],v[2])*dt
        mujoco.mj_step(m,d)
    T=sw.T-t_from; alph=np.array(alph); wts=np.array(wts)
    return dict(drag_thrust_N=Dx/T, lift_thrust_N=Lx/T,
                lift_share=Lx/(Lx+Dx) if (Lx+Dx)!=0 else float('nan'),
                alpha_thrust_weighted_deg=float((alph*wts).sum()/max(wts.sum(),1e-12)),
                path_transverse_ratio=vz/max(vx+vz,1e-12))
