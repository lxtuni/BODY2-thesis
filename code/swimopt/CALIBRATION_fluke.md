# Fluke model calibration status

The current `config_fluke.json` coefficients and `body2_fluke.xml` ballast are provisional. Optimization results are model predictions, not measured swimming speeds.

## Existing reference data

- Hull CFD (`robot_real/hull_study/results/hull_resistance.md`): M1 bare-hull drag is 3.3 mN at 0.1 m/s and 10.8 mN at 0.2 m/s at the report's fixed waterline. At the nominal fluke-model waterline the box hull is 78% immersed, so its effective frontal area is 0.004524 m². Fitting `F = 0.5 ρ Cd A U² + cv·0.78·U` to these two points gives `Cd_x ≈ 0.093`, `cv ≈ 0.0154 N·s/m`. These values are in `config_fluke.json`; they require rechecking once static trim and mass are measured. No 0.3 m/s hull CFD point was found.
- UVLM (`robot_real/leg_dynamics/results/papers/sculling_fsi.json`): this sweep includes frequencies, k, mean thrust and pitch amplitudes. It uses a different prescribed heave/pitch motion and foil geometry from `body2_fluke.xml`. A matched case with the same foil dimensions, frequency, shank angle and spring stiffness is needed before changing `foil.CLa` or `alpha_stall_deg` to meet a 20% tolerance.
- The XML currently assigns 0.45 kg to `trunk_hull` and 0.85 kg to `vis_ballast`. The task list does not provide a measured assembled mass or static O2 height. Set these values from a weighing and still-water measurement; then verify O2 relative to the waterline (target -0.011 m in the task list).

## Remaining measurements

1. Assembled robot mass including ballast, and the still-water O2 height.
2. Single-fluke mean thrust and passive pitch amplitude at a known frequency, shank amplitude and spring stiffness, or a UVLM case with matching geometry and motion.
3. Confirm which hull shape and waterline to use when fitting `trunk_hull` drag.

## Firmware integration

`export_fluke_gait.py` writes a checked seven-point trajectory and four phase offsets to `results_fluke_bounded/firmware_7point.json`. The firmware source present in the OpenFOAM workspace still exposes four `swim_cp*_yJ/zJ` points and the V2 IK. The V3 seven-point interface appears only in `firmware_v3/README_scull_v3.md` there. The export is therefore ready for a V3 firmware branch, but it cannot be pasted into the currently available V2 source as-is.

## Comparative experiments

The A/B/C/D paper comparison in `TODO_fluke.md` should use the same calibrated hull, foil and mass model and the same evaluation budget. Those calibration inputs are still missing, so no paper comparison is claimed here.
