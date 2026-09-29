"""Convert an optimized fluke gait into the firmware's seven control points."""
import json
import math
import sys

import numpy as np

from simulate import Swimmer, load_cfg
from leg_v3 import d, ik, O2E, EM

CAL = {
    "FL": (53.0, 1.0, 96.0, 1.0),
    "FR": (130.0, -1.0, 80.0, -1.0),
    "BL": (59.0, 1.0, 90.0, 1.0),
    "BR": (127.0, -1.0, 90.0, -1.0),
}
KNOTS = [0.0, 1/6, 1/3, 0.5, 0.53, 0.53 + 0.47/3, 0.53 + 2*0.47/3]


def main(config="config_fluke.json", best_path=None):
    cfg = load_cfg(config)
    best_path = best_path or cfg.get("outdir", "results_fluke") + "/best.json"
    sw = Swimmer(cfg["model"], cfg)
    with open(best_path, encoding="utf-8") as fh:
        p = sw.gait.decode(sw.gait.expand(json.load(fh)["x"]))
    if sw.gait.H != 1 or sw.gait.n != 4:
        raise ValueError("This exporter requires the four-shank, single-harmonic fluke model")
    freq = p["freq"]
    phases = {}
    offsets = {}
    amplitudes = {}
    for i, name in enumerate(sw.gait.names):
        leg = next((k for k in CAL if k in name), None)
        if leg is None:
            raise ValueError(f"Unknown leg in actuator {name}")
        phases[leg] = float(p["PH"][0][i] / (2*math.pi))
        offsets[leg] = float(p["offs"][sw.gait._off_slot[i]])
        amplitudes[leg] = float(p["A"][0][sw.gait._amp_slot[i]])
    if max(offsets.values()) - min(offsets.values()) > 1e-6 or max(amplitudes.values()) - min(amplitudes.values()) > 1e-6:
        raise ValueError("Firmware's shared seven points require shared amplitude and offset")
    mean = math.degrees(offsets["FL"])
    amp = math.degrees(amplitudes["FL"])
    E = O2E * d(-85.0)
    qrs = [mean - amp, mean - amp/2, mean + amp/2, mean + amp,
           mean + amp, mean + amp/2, mean - amp/2]
    if min(qrs) < float(cfg.get("fluke_qr_min_deg", 15.0)) or max(qrs) > math.degrees(1.25):
        raise ValueError("Optimized shank range would collide with hull or exceed joint limit")
    points = []
    for qr in qrs:
        M = E - EM * d(qr)
        s = ik(*M)
        if s is None or s["ACO2"] < 10 or s["O2EM"] < 45:
            raise ValueError(f"Four-bar IK invalid at q_r={qr:.1f}°")
        for leg, (o1, sign1, o2, sign2) in CAL.items():
            servo = (o1 + sign1*s["q1"], o2 + sign2*s["gamma"])
            if not all(0 <= angle <= 180 for angle in servo):
                raise ValueError(f"{leg} servo out of range at q_r={qr:.1f}°: {servo}")
        points.append([round(float(M[0]), 2), round(float(M[1]), 2)])
    # The firmware interpolates these points linearly at the stated phase knots.
    trace = []
    for j in range(1001):
        u = j/1000
        seg = max(i for i, knot in enumerate(KNOTS) if u >= knot)
        end = KNOTS[seg+1] if seg < 6 else 1.0
        a, b = np.array(points[seg]), np.array(points[(seg+1) % 7])
        M = a + (b-a)*(u-KNOTS[seg])/(end-KNOTS[seg])
        s = ik(*M)
        if s is None or s["ACO2"] < 10 or s["O2EM"] < 45:
            raise ValueError(f"Four-bar IK invalid at phase {u:.3f}")
        trace.append(s["q1"])
    rate = float(np.max(np.abs(np.diff(trace)))) * 1000 * freq
    if rate >= 820:
        raise ValueError(f"Peak servo rate {rate:.0f}°/s exceeds 820°/s")
    phase_offsets = {leg: (phase-phases["FL"]) % 1.0 for leg, phase in phases.items()}
    result = dict(freq_hz=freq, qr_mean_deg=mean, qr_amp_deg=amp,
                  points=points, knots=KNOTS, phase_offsets=phase_offsets,
                  max_servo_rate_deg_s=rate)
    out = cfg.get("outdir", "results_fluke") + "/firmware_7point.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
