#!/usr/bin/env python3
"""
OpenFOAM 船体阻力后处理
=======================
解析 forces 函数对象的输出，给出总阻力及其压差/摩擦分量，并画收敛曲线。

用法:
    python3 resistance.py ~/run/dtcHull
    python3 resistance.py ~/run/dtcHull --full-ship    # 半船模型 ×2 得整船

只用标准库就能算出数字；matplotlib 有的话额外存一张收敛曲线 PNG。
"""
import os
import re
import sys
import glob

TAIL_FRAC = 0.20   # 用最后这部分时间求平均，作为收敛值


def find_force_files(case_dir):
    """postProcessing 下可能有多个重启时间目录，全都收集。"""
    pats = [
        os.path.join(case_dir, "postProcessing", "**", "force.dat"),
        os.path.join(case_dir, "postProcessing", "**", "forces.dat"),
    ]
    files = []
    for p in pats:
        files.extend(glob.glob(p, recursive=True))
    return sorted(files)


def parse_header(lines):
    """从注释行里解析每组括号代表什么（total / pressure / viscous / porous）。"""
    for ln in lines:
        if not ln.startswith("#"):
            continue
        if "Time" not in ln:
            continue
        groups = re.findall(r"\(([^)]*)\)", ln)
        labels = []
        for g in groups:
            first = g.split()[0] if g.split() else ""
            labels.append(re.sub(r"_[xyz]$", "", first).lower())
        if not labels:                      # v2606 风格：total_x total_y total_z pressure_x ...
            labels = [m.lower() for m in re.findall(r"\b([A-Za-z]+)_x\b", ln)]
        if labels:
            return labels
    return []


def parse(files):
    times, data, labels = [], [], []
    for fn in files:
        with open(fn, "r", errors="replace") as fh:
            lines = fh.readlines()
        if not labels:
            labels = parse_header(lines)
        for ln in lines:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            nums = ln.replace("(", " ").replace(")", " ").split()
            try:
                vals = [float(x) for x in nums]
            except ValueError:
                continue
            if len(vals) < 4:
                continue
            times.append(vals[0])
            data.append(vals[1:])
    # 按时间排序去重（重启会产生重叠）
    seen, t_out, d_out = set(), [], []
    for t, d in sorted(zip(times, data), key=lambda x: x[0]):
        if t in seen:
            continue
        seen.add(t)
        t_out.append(t)
        d_out.append(d)
    return t_out, d_out, labels


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    full_ship = "--full-ship" in sys.argv
    case = os.path.expanduser(args[0]) if args else "."

    files = find_force_files(case)
    if not files:
        print(f"在 {case}/postProcessing 下找不到 force.dat / forces.dat。")
        print("可能原因：算例还没跑完，或 controlDict 里没有 forces 函数对象。")
        print("检查一下：  grep -n 'forces' %s/system/controlDict" % case)
        sys.exit(1)

    print("读取文件:")
    for f in files:
        print("   ", f.replace(os.path.expanduser("~"), "~"))

    t, d, labels = parse(files)
    if not t:
        print("\n文件里没有可解析的数据行。")
        sys.exit(1)

    ncomp = len(d[0]) // 3
    if not labels or len(labels) != ncomp:
        labels = ["total", "pressure", "viscous", "porous"][:ncomp]
        print(f"\n[注意] 表头没解析出来，按默认顺序假定为: {labels}")
    print(f"\n数据点 {len(t)} 个，时间 {t[0]:.4g} → {t[-1]:.4g} s，分量: {labels}")

    # 取尾段平均作为收敛值
    t_cut = t[-1] - (t[-1] - t[0]) * TAIL_FRAC
    idx = [i for i, tt in enumerate(t) if tt >= t_cut]
    if len(idx) < 2:
        idx = list(range(len(t)))
        t_cut = t[0]

    def mean_x(k):                       # 第 k 组的 x 分量（阻力方向）
        col = 3 * k
        return sum(d[i][col] for i in idx) / len(idx)

    def spread_x(k):
        col = 3 * k
        v = [d[i][col] for i in idx]
        return max(v) - min(v)

    vals = {lab: mean_x(k) for k, lab in enumerate(labels)}
    if "total" not in vals:
        vals["total"] = sum(vals.get(x, 0.0) for x in ("pressure", "viscous"))

    scale = 2.0 if full_ship else 1.0
    tag = "整船（半船模型 ×2）" if full_ship else "计算域内"

    print("\n" + "=" * 52)
    print(f"  阻力 Widerstand —— {tag}")
    print(f"  （对最后 {TAIL_FRAC:.0%} 时间段求平均: t >= {t_cut:.4g} s）")
    print("=" * 52)
    tot = vals["total"] * scale
    for lab in ("pressure", "viscous", "porous"):
        if lab in vals:
            v = vals[lab] * scale
            name = {"pressure": "压差阻力 (Druckwiderstand)",
                    "viscous":  "摩擦阻力 (Reibungswiderstand)",
                    "porous":   "多孔介质项"}[lab]
            frac = f"{v / tot * 100:5.1f}%" if tot else "   -- "
            print(f"  {name:<34} {v:12.3f} N   {frac}")
    print("-" * 52)
    print(f"  {'总阻力 (Gesamtwiderstand)':<34} {tot:12.3f} N")
    print("=" * 52)

    kt = labels.index("total") if "total" in labels else None
    if kt is not None and tot:
        col = 3 * kt
        vals_w = [d[i][col] * scale for i in idx]
        half = len(vals_w) // 2
        m1 = sum(vals_w[:half]) / max(half, 1)
        m2 = sum(vals_w[half:]) / max(len(vals_w) - half, 1)
        drift = (m2 - m1) / tot                    # 窗口前后半段均值之差 / 均值
        pp = (max(vals_w) - min(vals_w)) / abs(tot)  # 峰峰值
        print(f"\n收敛诊断（尾段窗口）：峰峰值 {pp*100:.1f}%，前后半段均值漂移 {drift*100:+.1f}%")
        if abs(drift) > 0.03:
            print("  → 均值还在漂移，没收敛。把 endTime 加大再跑，或检查库朗数。")
        elif pp > 0.05:
            print("  → 均值稳定但在振荡（钝头体尾流脱涡的典型表现）。取窗口平均值即可，")
            print("    想更准就把 endTime 延长、让窗口覆盖更多个振荡周期。")
        else:
            print("  → 收敛良好。")

    if not full_ship:
        print("\n提示：DTCHull 是沿中纵剖面切开的半船模型。要整船阻力加 --full-ship。")

    # 可选出图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n（想要收敛曲线图： sudo apt install -y python3-matplotlib）")
        return

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [2, 1]})
    styles = {"total": ("-", 2.0), "pressure": ("--", 1.3), "viscous": (":", 1.3)}
    series = {}
    for k, lab in enumerate(labels):
        if lab == "porous":
            continue
        series[lab] = [row[3 * k] * scale for row in d]
        ls, lw = styles.get(lab, ("-", 1.0))
        ax.plot(t, series[lab], ls, lw=lw, label=lab)
        ax2.plot(t, series[lab], ls, lw=lw, label=lab)
    # 左图：y 轴范围按起步瞬态之后（前 10% 时间除外）的数据定，不让冲击尖峰压扁曲线
    i0 = max(1, int(len(t) * 0.10))
    ys = [v for lab in series for v in series[lab][i0:]]
    if ys:
        lo, hi = min(ys), max(ys)
        pad = 0.15 * (hi - lo if hi > lo else abs(hi) or 1.0)
        ax.set_ylim(lo - pad, hi + pad)
    ax.axvline(t_cut, color="gray", lw=0.8, alpha=0.6)
    ax.text(t_cut, ax.get_ylim()[1], " averaging window", va="top", fontsize=8, color="gray")
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Force x [N]")
    ax.set_title(f"Resistance convergence — {os.path.basename(os.path.abspath(case))}")
    ax.legend(); ax.grid(alpha=0.3)
    # 右图：只看尾段窗口
    ax2.set_xlim(t_cut, t[-1])
    yw = [v for lab in series for v in series[lab][idx[0]:]]
    if yw:
        lo, hi = min(yw), max(yw)
        pad = 0.2 * (hi - lo if hi > lo else abs(hi) or 1.0)
        ax2.set_ylim(lo - pad, hi + pad)
    ax2.axhline(tot, color="k", lw=0.8, ls="--", alpha=0.6)
    ax2.set_title(f"tail window  (mean total = {tot:.4f} N)", fontsize=9)
    ax2.set_xlabel("Time [s]"); ax2.grid(alpha=0.3)
    out = os.path.join(case, "resistance.png")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"\n收敛曲线已存: {out}")
    print("看图：  explorer.exe .   然后双击 resistance.png")


if __name__ == "__main__":
    main()
