import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import brentq
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':10,'axes.labelsize':10,
    'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
    'axes.edgecolor':'#555','axes.linewidth':0.8,'xtick.color':'#555','ytick.color':'#555','axes.labelcolor':'#222'})
C={'lift':'#2a78d6','cur':'#eb6834','v3':'#1baf7a','ink':'#222','muted':'#777','grid':'#e6e6e3'}
# ---- CFD data (single leg, cycle-averaged thrust, mN)
lift=np.array([[0,45.5],[0.08,36.5],[0.15,29.1],[0.2229,21.4]])
v3=np.array([[0,87.2],[0.08,77.9],[0.15,59.0],[0.2229,35.9]])
cur=np.array([[0,23.5],[0.15,-13.1]])
def line_fit(d): return np.polyfit(d[:,0],d[:,1],1)
pl=line_fit(lift); pc=line_fit(cur)
def v3f(U):  # piecewise linear, linear extrapolation of last segment
    return np.interp(U,v3[:,0],v3[:,1]) if U<=v3[-1,0] else v3[-1,1]+(U-v3[-1,0])*(v3[-1,1]-v3[-2,1])/(v3[-1,0]-v3[-2,0])
fl=lambda U: np.polyval(pl,U); fc=lambda U: np.polyval(pc,U)
D=lambda U: 77.6*(U/0.2)**2          # whole-robot drag, CFD static float at 0.2 m/s, scaled with U^2 (upper bound)

def style(ax):
    ax.grid(True,color=C['grid'],lw=0.6); ax.set_axisbelow(True)
    ax.axhline(0,color='#999',lw=0.8)

# ================= Fig. A: thrust vs speed
fig,ax=plt.subplots(figsize=(6.0,3.7)); style(ax)
Ug=np.linspace(0,0.24,50)
for d,k,m,lab in [(lift,'lift','o','Lift-based fluke (L0 kinematics)'),(v3,'v3','^','Drag-based paddle, optimized (V3trap)'),(cur,'cur','s','Drag-based paddle, robot gait (current)')]:
    ax.plot(d[:,0],d[:,1],color=C[k],lw=2,marker=m,ms=7,mec='white',mew=1.2,label=lab,zorder=3)
# direct labels
ax.text(0.226,35.9+4,'V3trap',color=C['ink'],fontsize=9,ha='left')
ax.text(0.226,21.4-1,'fluke',color=C['ink'],fontsize=9,ha='left')
ax.text(0.153,-13.1-1,'robot gait',color=C['ink'],fontsize=9,ha='left',va='top')
ax.set_xlim(-0.005,0.26); ax.set_ylim(-30,100)
ax.set_xlabel(r'Forward speed $U$ [m s$^{-1}$]'); ax.set_ylabel(r'Mean thrust per leg $\bar{T}$ [mN]')
ax.legend(loc='upper right',frameon=False,fontsize=8.5,handlelength=2.4)
fig.tight_layout(); fig.savefig('/home/claude/thesis_figs/fig_thrust_speed.pdf'); fig.savefig('/home/claude/thesis_figs/fig_thrust_speed.png',dpi=220); plt.close(fig)

# ================= Fig. B: self-propulsion balance 4T(U) vs D(U)
Uc={}
Uc['lift']=brentq(lambda U:4*fl(U)-D(U),0.05,0.4)
Uc['cur']=brentq(lambda U:4*fc(U)-D(U),0.01,0.2)
Uc['v3']=brentq(lambda U:4*v3f(U)-D(U),0.05,0.4)
print({k:round(v,3) for k,v in Uc.items()})
fig,ax=plt.subplots(figsize=(6.0,3.7)); style(ax)
Ug=np.linspace(0,0.30,200)
ax.plot(Ug,D(Ug),color='#555',lw=1.6,ls='--',label='Whole-robot drag $D(U)$ (upper bound)',zorder=2)
m_in=Ug<=0.2229
for k,f,m,lab in [('lift',fl,'o','Lift-based fluke'),('v3',v3f,'^','Drag-based, optimized (V3trap)'),('cur',fc,'s','Drag-based, robot gait')]:
    y=np.array([4*f(u) for u in Ug])
    lim=0.2229 if k!='cur' else 0.15
    ax.plot(Ug[Ug<=lim],y[Ug<=lim],color=C[k],lw=2,label=r'$4\bar{T}$, '+lab,zorder=3)
    ax.plot(Ug[Ug>=lim],y[Ug>=lim],color=C[k],lw=1.3,ls=(0,(1,1.6)),zorder=3)   # extrapolation
    u=Uc[k]; ax.plot([u],[D(u)],'o',color=C[k],ms=7,mec='white',mew=1.2,zorder=4)
    off={'cur':(0.012,28),'lift':(0.004,-48),'v3':(0.010,40)}[k]
    ax.annotate(f'{u:.2f} m/s',xy=(u,D(u)),xytext=(u+off[0],D(u)+off[1]),fontsize=9,color=C['ink'],
                arrowprops=dict(arrowstyle='-',color='#aaa',lw=0.7))
ax.set_xlim(0,0.30); ax.set_ylim(-60,380)
ax.set_xlabel(r'Forward speed $U$ [m s$^{-1}$]'); ax.set_ylabel('Force [mN]')
ax.legend(loc='upper left',frameon=False,fontsize=8.3,handlelength=2.4)
ax.text(0.298,-52,'solid: CFD data range   dotted: extrapolated',fontsize=8,color=C['muted'],ha='right')
fig.tight_layout(); fig.savefig('/home/claude/thesis_figs/fig_self_propulsion.pdf'); fig.savefig('/home/claude/thesis_figs/fig_self_propulsion.png',dpi=220); plt.close(fig)
