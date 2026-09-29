# swimopt

**Model-agnostic gait optimization for paddling swimming robots, in MuJoCo.**

Give it any MuJoCo/URDF robot; it attaches a Morison-type hydrodynamic model to every
link, parameterizes the joint trajectories as a truncated Fourier series with a
duty-cycle warp, and runs CMA-ES to find the fastest / straightest / most efficient
swimming gait.

Built for the thesis *"Hydrodynamic Modeling and Trajectory Optimization of Paddling
Propulsion for a Quadruped Robot."*

中文使用说明见 [README.zh-CN.md](README.zh-CN.md)。

---

## Why

Most published paddling-robot studies hand-pick two or three gaits (diagonal,
front–back, wave) and compare them. The search space is much larger than that:
frequency, per-joint amplitude, phase offsets, mean angles, harmonic content and
power-stroke duty cycle all interact. `swimopt` searches that space automatically,
and is written so that **swapping the robot only touches layer ①** — useful when the
CAD model is still iterating.

## Architecture

```
robots/*.xml     ① model layer      ← the only thing that changes per robot
hydro.py         ② physics layer    robot-agnostic Morison hydrodynamics (vectorized)
gait.py          ③ control layer    robot-agnostic Fourier + duty-cycle gait
simulate.py      ──  glue: one rollout → one score
optimize.py      ④ optimization     CMA-ES
optimize_view.py ④ optimization     same, with live rendering
ui.py            Tkinter control panel (everything below, with checkboxes)
import_model.py  URDF → MJCF importer (portable meshes, auto actuators)
set_model.py     switch the active robot
view.py          replay a gait / the best result
config.json      all parameters, with inline comments
```

## Install

Requires Python 3.10–3.12.

```bash
git clone https://github.com/<you>/swimopt.git
cd swimopt
pip install -r requirements.txt
```

On Windows you can instead double-click `0_安装环境.bat`, which builds a local
virtualenv (`mjenv/`) and installs everything.

## Run

```bash
python view.py     config.json --demo    # watch a hand-made gait
python optimize.py config.json 250       # CMA-ES, 250 evaluations (~3 min)
python view.py     config.json --best    # watch the optimum
python ui.py                             # graphical control panel
```

Windows users: `1_看示例步态.bat` … `8_控制面板.bat` do the same thing, in order.

Results land in `results/` as `log.csv`, `best.json`, `convergence.json`.

## The physics

Per-link quasi-steady Morison-type force, scaled by immersion fraction `f`:

```
F = ρgV·f              buoyancy
  + ½ρ·C_d·A·v|v|·f    quadratic drag, anisotropic (per-axis C_d)
  + C_a·ρV·(dv/dt)     added mass
  + c_v·v·f            linear viscous
```

Notes on the implementation:

- **Partial immersion** is computed from each geom's true vertical extent under its
  current orientation, not a bounding sphere, so a flipper rolling at the surface
  gets the right `f`.
- **Anisotropic drag** is what produces thrust at all: the paddle face has a large
  `C_d`, the thin edge a small one, and the difference over a stroke cycle is the
  net force.
- **Added mass** uses the "mass trick" — add `m_a = C_a·ρV` to the body's inertial
  mass and cancel the resulting extra weight with a constant upward force. It is an
  isotropic approximation, but it is numerically unconditionally stable, which the
  explicit `dv/dt` form is not.
- **Not modelled:** turbulence, vortex shedding, wake interaction, free-surface
  waves. These matter for violent manoeuvres; for steady paddling their mean effect
  is absorbed into experimentally calibrated coefficients. Same modelling level as
  the beaver-robot literature this work builds on.

Links are matched to hydrodynamic rules by **name pattern** in `config.json`, so a
new robot needs no code changes — just rules like `"*flipper*": {cd: [1.2, 0.08, 0.08]}`.

## The gait parameterization

```
q_i(t) = offset_i + Σ_h  A_{i,h} · sin(2π·h·f·t + φ_{i,h})       h = 1, 2
```
plus a duty-cycle phase warp that lets the power stroke occupy a different fraction
of the cycle than the recovery stroke.

**The second harmonic is not optional.** With a first harmonic only, two legs driven
in antiphase produce exactly opposite thrust, all four legs cancel, and net thrust is
identically zero — verified numerically, not just argued. A second harmonic is
invariant under a π phase shift, which is the mathematical counterpart of the real
robot's passive flipper feathering (spread on the power stroke, furl on recovery).
It is the prerequisite for any net thrust at all.

Any parameter group can be frozen from the UI or config. Locking phase and offset
collapses a 36-dimensional search to 8 dimensions, which is how you compare a fixed
gait pattern (diagonal vs front–back vs wave) at its own best frequency and amplitude.

## The optimizer

CMA-ES (Hansen & Ostermeier 2001), via the `cma` package. It is a **gradient-free
evolution strategy, not reinforcement learning** — it optimizes a fixed-form
trajectory's parameters, not a state-feedback policy.

```
fitness = speed − w_yaw · yaw_rate − w_energy · mean_power
```

The weights are a ruler, not a result: fitness values are not comparable across
different weight settings, so report the physical quantities (m/s, BL/s, degrees, W),
never the fitness number. Presets in the UI cover *fastest*, *fast and straight*,
*strictly straight*, *power-thrifty*.

Convergence means the fitness curve flattens **and** the CMA-ES step size σ shrinks.
Re-running with a different seed should land nearby; if it doesn't, the budget was
too small.

## What is not in this repository

- The **BODY2 CAD model** (`robots/body2.xml` and its STL meshes) — lab property,
  gitignored. `robots/toy_quad.xml` is a small open four-legged swimmer that
  exercises the whole pipeline and is enough to reproduce every feature.
- Reference papers — copyrighted.
- `mjenv/` — machine-specific virtualenv.

To use your own robot: `python import_model.py my_robot.urdf robots/mine.xml --drive "2.1,1.1"`,
then `python set_model.py`. The importer copies meshes locally and rewrites paths so
the model stays portable across machines.

Note that URDF cannot express closed kinematic loops. A parallelogram or five-bar leg
will import as an open tree and produce no thrust; close it manually with MuJoCo
`<equality><connect>` after import.

## Status

With `<compiler angle="radian"/>` fixed in `toy_quad.xml`, the hand-made demo gait
simulates at 0.0706 m/s. A 250-evaluation run (seed 1) found 0.1073 m/s in 179 s,
with 42.7° yaw drift over 8 s. This run has not converged; the large drift makes it
unsuitable as a straight-swimming reference without further tuning.
Rollouts are bit-reproducible for identical parameters.

Coefficients are currently literature values. Experimental calibration against the
physical robot (coast-down, static draft, pendulum decay, biped paddling) is the
next step, and is what the accuracy ultimately rests on — not the choice of simulator.

## References

- Morison, J.R. et al. (1950). The force exerted by surface waves on piles. *J. Pet. Technol.* 2(5).
- Fossen, T.I. (1994). *Guidance and Control of Ocean Vehicles*. Wiley.
- Chen, G. et al. (2023). Modeling of swimming posture dynamics for a beaver-like robot. *Ocean Engineering*.
- Chen, G. et al. (2021). Hydrodynamic analysis of webbed feet. *Ocean Engineering* 234:109179.
- Hansen, N. & Ostermeier, A. (2001). Completely derandomized self-adaptation in evolution strategies. *Evolutionary Computation* 9(2):159–195.
- Lee, S. et al. (2025). Unified fluid–robot multiphysics simulation and optimization. arXiv:2506.05012. *(workflow reference; they use a gradient method, this project uses CMA-ES because MuJoCo is not differentiable)*

## License

MIT — see [LICENSE](LICENSE).
