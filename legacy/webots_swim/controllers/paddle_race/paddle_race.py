# -*- coding: utf-8 -*-
"""
paddle_race.py  --  one controller, gait chosen by controllerArgs.
Leg numbering (4 corners):  1=FL  2=FR  3=BL  4=BR
  "14_23" : legs 1&4 (FL,BR) together, 2&3 (FR,BL) together  -> diagonal / trot-like
  "12_34" : legs 1&2 (FL,FR) together, 3&4 (BL,BR) together  -> front-pair vs back-pair (bound-like)
Both gaits use the SAME paddle mechanism, so the comparison is fair: only the phasing differs.
Each robot prints its own name, gait, displacement from start, and speed.
"""
from controller import Robot
import sys, math

robot = Robot()
dt = int(robot.getBasicTimeStep())
gait = sys.argv[1] if len(sys.argv) > 1 else "14_23"
name = robot.getName()

if gait == "12_34":
    PHASE = {"FL": 0.0, "FR": 0.0, "BL": math.pi, "BR": math.pi}
else:  # "14_23"
    PHASE = {"FL": 0.0, "BR": 0.0, "FR": math.pi, "BL": math.pi}

LEGS    = ["FL", "FR", "BL", "BR"]
FREQ    = 0.8
HIP_AMP = 0.50
FEATHER = 1.20
RAMP_T  = 3.0

hip   = {L: robot.getDevice("%s_hip"   % L) for L in LEGS}
wrist = {L: robot.getDevice("%s_wrist" % L) for L in LEGS}
for L in LEGS:
    hip[L].setPosition(0.0); wrist[L].setPosition(0.0)
gps = robot.getDevice("gps"); gps.enable(dt)

x0 = None; t_rep = 0.0
print("[%s] gait=%s  start (freq=%.2fHz, 3s ease-in)" % (name, gait, FREQ))

W = 2.0 * math.pi * FREQ
while robot.step(dt) != -1:
    t = robot.getTime()
    ramp = min(1.0, t / RAMP_T)
    for L in LEGS:
        ph = W * t + PHASE[L]
        hip[L].setPosition(ramp * HIP_AMP * math.sin(ph))
        wrist[L].setPosition(ramp * FEATHER * max(0.0, math.cos(ph)))
    p = gps.getValues()
    if x0 is None:
        x0, y0 = p[0], p[1]; t0 = t
    elif (t - t_rep) >= 2.0:
        dx = p[0] - x0; dy = p[1] - y0
        dist = math.hypot(dx, dy)
        v = dist / (t - t0) if t > t0 else 0.0
        print("[%s] gait=%s  t=%5.1fs  travel=%.3f m  avg_speed=%.4f m/s" % (name, gait, t, dist, v))
        t_rep = t
