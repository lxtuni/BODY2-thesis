"""
paddle_swim.py  --  Webots controller for the robot0526 swim demo.

Drag-based "dog-paddle" gait (cf. Li 2019 / Qu 2025), diagonal pairs (1-4 / 2-3):
  FL+RR stroke together, then FR+RL  -> continuous, balanced thrust.

  - power stroke : leg sweeps BACKWARD, paddle held broadside -> high drag -> forward thrust
  - recovery     : paddle FEATHERS edge-on, leg swings forward -> low drag

Each leg has 2 motors:  "<LEG>_hip"  (swings the leg fore/aft, axis Y)
                        "<LEG>_wrist"(feathers the paddle,    axis Y)

Buoyancy + drag are produced by Webots (Fluid + per-link ImmersionProperties);
this controller only commands the joints. A 3 s ease-in ramp avoids any start-up jerk.
"""

from controller import Robot
import math

robot = Robot()
dt = int(robot.getBasicTimeStep())

# ---------------- gait parameters (tune these) ----------------
FREQ      = 0.7      # paddling frequency [Hz]
HIP_AMP   = 0.50     # hip swing amplitude [rad]  (~29 deg each way)
FEATHER   = 1.20     # paddle feather angle during recovery [rad] (~69 deg)
RAMP_T    = 3.0      # seconds to ease the amplitude from 0 -> full
# diagonal (trot-like) gait: legs 1&4 (FL,RR) together, 2&3 (FR,RL) together
PHASE = {"FL": 0.0, "RR": 0.0, "FR": math.pi, "RL": math.pi}
LEGS  = ["FL", "FR", "RL", "RR"]

# ---------------- grab devices ----------------
hip   = {L: robot.getDevice("%s_hip" % L)   for L in LEGS}
wrist = {L: robot.getDevice("%s_wrist" % L) for L in LEGS}
for L in LEGS:
    hip[L].setPosition(0.0)
    wrist[L].setPosition(0.0)

gps = robot.getDevice("gps")
if gps:
    gps.enable(dt)
imu = robot.getDevice("imu")
if imu:
    imu.enable(dt)

x_prev, t_prev = None, 0.0
print("[paddle_swim] start  freq=%.2fHz hip=%.2frad feather=%.2frad (3s ease-in)" % (FREQ, HIP_AMP, FEATHER))

# ---------------- main loop ----------------
while robot.step(dt) != -1:
    t = robot.getTime()
    ramp = min(1.0, t / RAMP_T)          # ease-in 0 -> 1 over RAMP_T seconds
    w = 2.0 * math.pi * FREQ
    for L in LEGS:
        ph = w * t + PHASE[L]
        hip[L].setPosition(ramp * HIP_AMP * math.sin(ph))
        # feather only during the forward (recovery) half-cycle: cos(ph) > 0
        wrist[L].setPosition(ramp * FEATHER * max(0.0, math.cos(ph)))

    # report forward speed once per second
    if gps:
        x = gps.getValues()[0]
        if x_prev is None:
            x_prev, t_prev = x, t
        elif (t - t_prev) >= 1.0:
            yaw = math.degrees(imu.getRollPitchYaw()[2]) if imu else 0.0
            v = (x - x_prev) / (t - t_prev)
            print("[paddle_swim] t=%5.1fs  x=%+.3f m  v=%+.3f m/s  yaw=%+.0f deg" % (t, x, v, yaw))
            x_prev, t_prev = x, t
