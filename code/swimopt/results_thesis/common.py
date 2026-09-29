# Shared helpers for thesis experiments (all outputs stay in results_thesis/)
import os, sys, json, copy, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.chdir(ROOT)
os.environ.setdefault("MUJOCO_GL", "disable")
import numpy as np
from simulate import Swimmer, load_cfg

def toy_cfg(**kw):
    cfg = load_cfg("config.json")
    cfg["model"] = "robots/toy_quad.xml"
    cfg.pop("foil", None)
    cfg.update(kw)
    return cfg

def fluke_cfg(**kw):
    cfg = load_cfg("config_fluke.json")
    cfg.update(kw)
    return cfg

def save(name, obj):
    with open(os.path.join(HERE, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=float)
