# -*- coding: utf-8 -*-
"""
force_rig.py  --  one paddle on a torque-feedback motor, for three hydro experiments.

  tau_motor = I_hinge*alpha + tau_hydro(Webots)  =>  tau_hydro = tau_motor - I_hinge*alpha

MODE = "drag"      Exp1 - DRAG CALIBRATION.
                   Set the Fluid.streamVelocity in force_rig.wbt to  (V_STREAM 0 0),
                   the paddle is HELD broadside; measured holding torque -> effective Cd.
MODE = "oscillate" Exp2 + Exp3 - DYNAMIC + ADDED MASS.
                   Sinusoid. At |omega| max (theta~0) the force is pure DRAG;
                   at omega~0 (theta=peak) it is pure ADDED MASS. The console reports
                   both bins, so you see (a) drag matches blade-element and
                   (b) the added-mass term is missing (or restored by the plugin).

Writes rig_log.csv in the controller folder. Geometry must match force_rig.wbt.
"""
from controller import Robot
import math, os

MODE     = "oscillate"   # "drag" or "oscillate"
# -- oscillate params --
F   = 0.6                # frequency [Hz]
AMP = 0.40               # amplitude [rad]
# -- drag params (must equal Fluid.streamVelocity you set in the world) --
V_STREAM = 0.20          # [m/s]

# geometry / fluid (match force_rig.wbt)
rho = 1000.0; Cd_set = 1.10
pw, pl, pt = 0.050, 0.070, 0.006
A = pw*pl; R = 0.10; m = 0.021
I_hinge = m/12.0*(pt*pt + pl*pl) + m*R*R
m_add   = rho*math.pi*(pw/2.0)**2*pl

robot = Robot(); dt = int(robot.getBasicTimeStep())
mot = robot.getDevice("rig_motor"); ps = robot.getDevice("rig_sensor")
mot.enableTorqueFeedback(dt); ps.enable(dt)
try: mot.setControlPID(120.0, 0.0, 3.0)
except Exception: pass

log = open(os.path.join(os.getcwd(), "rig_log.csv"), "w")
log.write("t,theta,omega,alpha,tau_motor,tau_inertia,tau_hydro_meas,tau_drag_analytic,tau_added_analytic\n")

W = 2.0*math.pi*F
omega_pk = AMP*W; alpha_pk = AMP*W*W
print("="*70)
print(" FORCE RIG  mode=%s   I_hinge=%.3e   m_add=%.1f g (%.1fx paddle)" % (MODE, I_hinge, m_add*1000, m_add/m))
if MODE == "drag":
    print(" DRAG CALIBRATION: ensure force_rig.wbt Fluid.streamVelocity = %.2f 0 0" % V_STREAM)
    mot.setPosition(0.0)
print("="*70)

# accumulators
db_meas=db_drag=db_n=0.0          # drag bin (|omega| near peak)
ab_meas=ab_add=ab_n=0.0           # added-mass bin (omega near 0)
dr_tau=dr_n=0.0                   # drag-mode holding torque

while robot.step(dt) != -1:
    t = robot.getTime()
    if MODE == "oscillate":
        amp = AMP*min(1.0, t/3.0)
        theta = amp*math.sin(W*t); mot.setPosition(theta)
        omega = amp*W*math.cos(W*t)
        alpha = -amp*W*W*math.sin(W*t)
    else:  # drag: held at 0, water flows past
        theta = 0.0; omega = 0.0; alpha = 0.0

    tau_motor = mot.getTorqueFeedback()
    tau_inertia = I_hinge*alpha
    tau_hydro = tau_motor - tau_inertia
    v = omega*R
    tau_drag  = 0.5*rho*Cd_set*A*v*abs(v)*R
    tau_added = m_add*alpha*R*R
    log.write("%.4f,%.5f,%.5f,%.5f,%.6f,%.6f,%.6f,%.6f,%.6f\n" %
              (t, theta, omega, alpha, tau_motor, tau_inertia, tau_hydro, tau_drag, tau_added)); log.flush()

    if t < 4.0: continue
    if MODE == "drag":
        # steady holding torque from the streamVelocity drag on the broadside paddle
        F_meas = abs(tau_motor)/R
        Cd_eff = F_meas/(0.5*rho*A*V_STREAM**2) if V_STREAM>0 else 0.0
        dr_tau += Cd_eff; dr_n += 1
        if int(t*1000) % 1000 < dt:
            print("[drag] t=%4.0fs  hold torque=%.4f Nm  F=%.4f N  ->  Cd_eff=%.2f  (set=%.2f)"
                  % (t, abs(tau_motor), F_meas, dr_tau/dr_n, Cd_set))
    else:
        if abs(omega) > 0.9*omega_pk:                 # near max speed -> DRAG
            db_meas += abs(tau_hydro); db_drag += abs(tau_drag); db_n += 1
        if abs(omega) < 0.15*omega_pk and abs(alpha) > 0.85*alpha_pk:  # near v=0 -> ADDED MASS
            ab_meas += abs(tau_hydro); ab_add += abs(tau_added); ab_n += 1
        if int(t*1000) % 2000 < dt and db_n>0 and ab_n>0:
            print("[osc] DRAG bin (v-peak): measured=%.4f  blade-element=%.4f  (ratio %.2f)"
                  % (db_meas/db_n, db_drag/db_n, (db_meas/db_n)/max(db_drag/db_n,1e-9)))
            print("      ADDED-MASS bin (v~0): measured=%.4f  analytic added-mass=%.4f"
                  % (ab_meas/ab_n, ab_add/ab_n))
            print("      -> 没插件时 measured(v~0)应≈0 (缺口); 加 physics \"added_mass\" 后应≈analytic")
            db_meas=db_drag=db_n=0.0; ab_meas=ab_add=ab_n=0.0
