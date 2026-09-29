/*
 * added_mass.c  --  Webots physics plugin
 * Adds the ADDED-MASS (acceleration-reaction) force that the built-in Fluid model
 * omits, on each paddle, normal to its broad face:
 *
 *     F_added = -m_add * a_normal            (Han 2025 / Morison inertia term)
 *     m_add   = rho * pi * (w/2)^2 * L       (flat-plate added mass, normal to face)
 *
 * The paddle's broad-face normal is its local X axis (the paddle box is thin in X).
 * Works on any world that contains DEF FL_PADDLE / FR_PADDLE / RL_PADDLE / RR_PADDLE
 * (the swimmer) and/or DEF RIG_PADDLE (the force rig). Missing DEFs are skipped.
 *
 * Enable it by setting   WorldInfo { physics "added_mass" }   and compiling (see Makefile).
 */
#include <ode/ode.h>
#include <plugins/physics.h>
#include <math.h>

#define RHO       1000.0
#define PI        3.14159265358979
#define PADDLE_W  0.050      /* paddle width  (Y) */
#define PADDLE_L  0.070      /* paddle length (Z) */
#define GAIN      1.0        /* scale the added-mass force (lower if it gets stiff) */
#define ACC_LP    0.7        /* acceleration low-pass factor [0..1), higher = smoother */

#define NPAD 5
static const char *DEFS[NPAD] = {"FL_PADDLE", "FR_PADDLE", "RL_PADDLE", "RR_PADDLE", "RIG_PADDLE"};
static dBodyID body[NPAD];
static double  vn_prev[NPAD];
static double  an_filt[NPAD];
static int     have_prev[NPAD];
static double  m_add;
static double  t_prev_ms;

void webots_physics_init() {
  m_add = RHO * PI * (PADDLE_W * 0.5) * (PADDLE_W * 0.5) * PADDLE_L;
  for (int i = 0; i < NPAD; i++) {
    body[i]      = dWebotsGetBodyFromDEF(DEFS[i]);   /* NULL if not in this world */
    have_prev[i] = 0;
    vn_prev[i]   = 0.0;
    an_filt[i]   = 0.0;
  }
  t_prev_ms = dWebotsGetTime();
  dWebotsConsolePrintf("[added_mass] init: m_add = %.1f g per paddle (GAIN=%.2f)",
                       m_add * 1000.0, GAIN);
}

void webots_physics_step() {
  double t_ms = dWebotsGetTime();
  double dt = (t_ms - t_prev_ms) / 1000.0;      /* ms -> s */
  t_prev_ms = t_ms;
  if (dt <= 0.0) dt = 1e-3;

  for (int i = 0; i < NPAD; i++) {
    dBodyID b = body[i];
    if (!b)
      continue;
    const dReal *v = dBodyGetLinearVel(b);       /* world velocity of COM */
    const dReal *R = dBodyGetRotation(b);        /* body->world rotation (4x3) */
    double nx = R[0], ny = R[4], nz = R[8];       /* local X axis = broad-face normal */
    double vn = v[0] * nx + v[1] * ny + v[2] * nz;
    if (have_prev[i]) {
      double an_raw = (vn - vn_prev[i]) / dt;
      an_filt[i] = ACC_LP * an_filt[i] + (1.0 - ACC_LP) * an_raw;
      double F = -GAIN * m_add * an_filt[i];      /* opposes normal acceleration */
      dBodyAddForce(b, F * nx, F * ny, F * nz);
    } else {
      have_prev[i] = 1;
    }
    vn_prev[i] = vn;
  }
}

void webots_physics_cleanup() {
}
