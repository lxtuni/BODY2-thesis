#!/usr/bin/env python3
"""
按姿态拼装 CFD 几何 (posture assembly)：密封身体 + 盖板 + 放下的腿（只取会碰水的零件）
=====================================================================================
    python3 build_hull.py --legs 1,3                 # 1、3 号腿放下（默认 CAD 编号：1=左前 2=左后 3=右前 4=右后）
    python3 build_hull.py --legs FL,FR               # 同上，用位置代号写
    python3 build_hull.py --legs FL,RR               # 对角（左前 + 右后）
    python3 build_hull.py --legs 1,3 --map 1=FL,2=FR,3=RR,4=RL   # 你们的编号和 CAD 不一样时自己给映射
    python3 build_hull.py --legs 1,3 --rotate FL:-20 --rotate FR:-20   # 把放下的腿整体绕髋轴再转 -20°
    python3 build_hull.py --legs FL,FR,RL,RR         # 四条腿全放下
    python3 build_hull.py --legs FL,FR --all-parts   # 连水面以上的舵机、架子也带上（只用来看，别拿去划网格）

腿的代号 (leg IDs)：FL 左前  FR 右前  RL 左后  RR 右后 —— x 正方向是前（小头/艏），y 正方向是左。
收起来 (retracted) 的腿默认不放进几何：它们在空气里，对水阻力没有贡献。
--rotate LEG:deg 把该腿所有零件绕"髋轴"(hip axis, 平行 y 轴、过从动杆上端销孔) 刚体旋转，
    角度正 = 桨尖向前 (+x) 摆，负 = 向后摆。这是姿态近似，不是连杆机构的真实运动。
需要：parts/index.json、parts/*.stl、sets.json（解压 robot_real_parts.zip 即可）以及 numpy。
"""
import argparse, json, os, struct, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LEG_IDS = ("FL", "FR", "RL", "RR")
CN = {"FL": "左前", "FR": "右前", "RL": "左后", "RR": "右后"}
DEFAULT_MAP = "1=FL,2=RL,3=FR,4=RR"      # STEP 里各零件实例号 _1.._4 对应的位置


def read_stl(p):
    with open(p, "rb") as f:
        f.seek(80)
        n = struct.unpack("<I", f.read(4))[0]
        dt = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        return np.frombuffer(f.read(50 * n), dtype=dt, count=n)["v"].astype(np.float64)


def write_stl(tri, fn):
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
    rec = np.zeros(len(tri), dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
    rec["n"] = n
    rec["v"] = tri
    with open(fn, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(tri)))
        f.write(rec.tobytes())


def topology(tri):
    """(封闭?, 法向一致?, 边界边数, 非流形边数) —— 多个独立封闭零件拼在一起也算封闭"""
    v = tri.reshape(-1, 3)
    key = np.round(v * 1e7).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_inverse=True)
    f = idx.reshape(-1, 3)
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    _, cnt = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    _, dcnt = np.unique(e, axis=0, return_counts=True)
    b, nm = int((cnt == 1).sum()), int((cnt > 2).sum())
    return b == 0 and nm == 0, bool(np.all(dcnt == 1)), b, nm


def leg_of(tri):
    c = tri.reshape(-1, 3).mean(0)
    return ("F" if c[0] > 0 else "R") + ("L" if c[1] > 0 else "R")


def rotate_y(tri, deg, px, pz):
    """绕过 (px, *, pz) 的 y 轴旋转；deg>0 时 z 以下的点向 +x 摆（桨尖向前）"""
    t = np.deg2rad(deg)
    x, z = tri[..., 0] - px, tri[..., 2] - pz
    out = tri.copy()
    out[..., 0] = px + x * np.cos(t) - z * np.sin(t)
    out[..., 2] = pz + x * np.sin(t) + z * np.cos(t)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--legs", default="FL,FR,RL,RR", help="放下的腿：FL/FR/RL/RR 或 1-4（用 --map 解释数字）")
    ap.add_argument("--map", default=DEFAULT_MAP, help=f"数字编号→位置，默认 {DEFAULT_MAP}（STEP 里的实例号）")
    ap.add_argument("--rotate", action="append", default=[], metavar="LEG:deg",
                    help="把某条腿绕髋轴刚体旋转 deg 度（可重复给）")
    ap.add_argument("--pivot", action="append", default=[], metavar="LEG:x,z",
                    help="自己指定该腿的旋转轴位置 (x,z)，默认取从动杆上端销孔")
    ap.add_argument("--all-parts", action="store_true", help="放下的腿连舵机/架子等水面以上零件也包含（只用于查看）")
    ap.add_argument("--wl", type=float, default=None,
                    help="给定水线 z (m)，按 zmin < wl+5mm 重新判定哪些腿部零件碰水（不再用 sets.json 里按 -0.022 算的集合）")
    ap.add_argument("--out", default=None, help="输出文件名，默认 hull_<腿>.stl")
    a = ap.parse_args()

    num2leg = {}
    for kv in a.map.split(","):
        k, v = kv.split("=")
        num2leg[k.strip()] = v.strip().upper()
    legs = []
    for s in a.legs.split(","):
        s = s.strip().upper()
        if not s:
            continue
        s = num2leg.get(s, s)
        if s not in LEG_IDS:
            sys.exit(f"腿代号只能是 FL/FR/RL/RR 或 --map 里的数字，收到 {s!r}")
        if s not in legs:
            legs.append(s)
    rot = {}
    for r in a.rotate:
        k, v = r.split(":")
        k = num2leg.get(k.strip().upper(), k.strip().upper())
        rot[k] = float(v)
    piv = {}
    for r in a.pivot:
        k, v = r.split(":")
        k = num2leg.get(k.strip().upper(), k.strip().upper())
        piv[k] = tuple(float(x) for x in v.split(","))

    meta = json.load(open(os.path.join(HERE, "parts", "index.json"), encoding="utf-8"))
    sets = json.load(open(os.path.join(HERE, "sets.json"), encoding="utf-8"))
    WL = a.wl if a.wl is not None else sets.get("WL_est_190g")
    if a.wl is not None:
        wet = None            # 按水线动态判定
    else:
        wet = set(os.path.basename(p) for p in sets["wet"])

    parts = []            # (file, leg, tri, is_leg_part)
    for m in meta:
        base = os.path.basename(m["file"])
        tri = read_stl(os.path.join(HERE, m["file"]))
        nm = m["name"]
        if nm.startswith("密封身体"):
            parts.append((base, "身体", tri, False))
        else:
            is_leg = not (nm.startswith("0526servo") or nm.startswith("架子"))
            parts.append((base, leg_of(tri), tri, is_leg))

    # 髋轴：从动杆上端销孔（取上端 4 mm 内顶点的平均）
    hip = {}
    for base, leg, tri, _ in parts:
        if base.startswith("大腿大腿从动杆"):
            v = tri.reshape(-1, 3)
            top = v[v[:, 2] > v[:, 2].max() - 0.004]
            hip[leg] = (float(top[:, 0].mean()), float(v[:, 2].max() - 0.002))
    hip.update(piv)

    chosen, skipped = [], []
    for base, leg, tri, is_leg in parts:
        if leg == "身体":
            chosen.append((base, leg, tri)); continue
        if leg not in legs:
            skipped.append((base, leg, "腿收起，在空气里")); continue
        if not is_leg and not a.all_parts:
            skipped.append((base, leg, "舵机/架子，水面以上")); continue
        if is_leg and not a.all_parts:
            is_wet = (tri[..., 2].min() < WL + 0.005) if wet is None else (base in wet)
            if not is_wet:
                skipped.append((base, leg, "水面以上")); continue
        if leg in rot and rot[leg] != 0.0:
            tri = rotate_y(tri, rot[leg], *hip[leg])
        chosen.append((base, leg, tri))
    chosen.append(("lid_1.stl", "盖板", read_stl(os.path.join(HERE, sets["lid"]))))

    tag = "_".join(legs) + "".join(f"_{l}{rot[l]:+g}" for l in legs if l in rot and rot[l] != 0.0)
    out = a.out or os.path.join(HERE, f"hull_{tag}{'.all' if a.all_parts else ''}.stl")
    tri = np.concatenate([c[2] for c in chosen])
    write_stl(tri, out)

    print(f"放下的腿: {', '.join(f'{l}({CN[l]})' for l in legs)}   收起的腿: "
          f"{', '.join(f'{l}({CN[l]})' for l in LEG_IDS if l not in legs) or '无'}")
    for l in legs:
        if l in rot and rot[l] != 0.0:
            print(f"   {l} 绕髋轴 (x={hip[l][0]:.4f}, z={hip[l][1]:.4f}) 旋转 {rot[l]:+g}°")
    print(f"写入 {out}：{len(tri)} 个三角形，{len(chosen)} 个零件")
    for base, leg, _ in chosen:
        print(f"   + {leg:<4} {base}")
    print(f"未包含 {len(skipped)} 个：")
    for base, leg, why in skipped:
        print(f"   - {leg:<4} {base}  ({why})")
    v = tri.reshape(-1, 3)
    print("包围盒 min", np.round(v.min(0), 4), "max", np.round(v.max(0), 4))
    if WL is not None:
        print(f"水线 z = {WL}；最低点 z = {v[:, 2].min():.4f}，水下深度 {WL - v[:, 2].min():.3f} m")
    closed, oriented, nb, nnm = topology(tri)
    print("拓扑检查:", "封闭 OK" if closed else f"!! 不封闭（{nb} 条边界边，{nnm} 条非流形边）",
          "| 法向一致 OK" if oriented else "| !! 法向不一致")
    if not closed:
        print("   零件之间有贴合面时可能误报；先用 surfaceCheck 再决定。")


if __name__ == "__main__":
    main()
