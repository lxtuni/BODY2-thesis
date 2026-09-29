# -*- coding: utf-8 -*-
"""
① 模型层导入工具 —— 把 URDF 转成本管线可用的 MJCF.

自动完成:
  1. 修复 package:// 网格路径 (MuJoCo 读不了这种路径)
  2. 注入 <mujoco><compiler> 设置 (网格目录/惯量修正)
  3. 为指定关节自动添加位置执行器 (URDF 里没有执行器概念)
  4. 输出报告: 有哪些刚体/关节/几何体, 以及水动力系数会怎么分配

用法:
  python import_model.py  你的模型.urdf  robots/输出名.xml  --drive "2.1,1.1"
      --drive  : 关节名里含这些关键字的 → 加执行器(逗号分隔). 不写则所有可动关节都加
      --meshdir: 网格文件夹(默认自动猜 urdf 同级的 meshes/)
"""
import os, re, sys, shutil, argparse
import numpy as np
import mujoco


def find_meshdir(urdf_path, text):
    """从 URDF 里挑一个网格文件名, 在磁盘上搜出它真正在哪, 反推 meshdir."""
    names = re.findall(r'filename="([^"]+)"', text)
    if not names:
        return "."
    rel = names[0].lstrip("./")                    # 如 meshes/base_link.STL
    fname = os.path.basename(rel)
    root = os.path.dirname(os.path.abspath(urdf_path))
    for base in [root, os.path.dirname(root), os.path.dirname(os.path.dirname(root))]:
        for dirpath, _, files in os.walk(base):
            if fname in files:
                # meshdir 应使 (meshdir + rel) 指向真实文件
                real = os.path.join(dirpath, fname)
                sub = os.path.dirname(rel)         # "meshes" 或 ""
                cand = os.path.abspath(os.path.join(dirpath, "..")) if sub else dirpath
                if os.path.isfile(os.path.join(cand, rel)):
                    return cand
                if os.path.isfile(os.path.join(dirpath, fname)) and not sub:
                    return dirpath
        # 只搜前两层, 避免太慢
    return root


def fix_urdf(src, meshdir):
    """返回修好的 URDF 文本: package:// → 相对路径, 并注入 mujoco 编译选项."""
    s = open(src, "r", encoding="utf-8", errors="ignore").read()
    # package://任意包名/xxx  →  xxx (只保留 meshes/ 之后的部分)
    s = re.sub(r'package://[^/]+/', '', s)
    if "<mujoco>" not in s:
        inject = (f'  <mujoco>\n'
                  f'    <compiler meshdir="{meshdir}" balanceinertia="true" '
                  f'discardvisual="false" fusestatic="false" strippath="false"/>\n'
                  f'  </mujoco>\n')
        s = re.sub(r'(<robot[^>]*>)', r'\1\n' + inject, s, count=1)
    return s


def add_actuators(mjcf_path, drive_keys, kp=60.0, forcerange=12.0, ctrlrange=1.4):
    """在 MJCF 里为匹配的关节添加 position 执行器."""
    s = open(mjcf_path, encoding="utf-8").read()
    joints = re.findall(r'<joint[^>]*name="([^"]+)"[^>]*>', s)
    # 排除 free/ball 关节
    sel = []
    for j in joints:
        if drive_keys and not any(k.strip().lower() in j.lower() for k in drive_keys):
            continue
        sel.append(j)
    if not sel:
        print("[!] 没有匹配到任何关节, 不添加执行器"); return s, []
    acts = "\n".join(
        f'    <position name="{j}" joint="{j}" kp="{kp}" dampratio="1" '
        f'ctrlrange="-{ctrlrange} {ctrlrange}" forcerange="-{forcerange} {forcerange}"/>'
        for j in sel)
    block = f"  <actuator>\n{acts}\n  </actuator>\n"
    if "<actuator>" in s:
        s = re.sub(r'</actuator>', acts + "\n  </actuator>", s, count=1)
    else:
        s = s.replace("</mujoco>", block + "</mujoco>")
    open(mjcf_path, "w", encoding="utf-8").write(s)
    return s, sel


def make_portable(mjcf_path):
    """把 MJCF 里引用的网格文件复制到模型旁边的 meshes/, 并改成相对路径.
    否则模型换台电脑就打不开(MuJoCo 保存的是绝对路径)."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(mjcf_path); root = tree.getroot()
    outdir = os.path.dirname(os.path.abspath(mjcf_path))
    mdir = os.path.join(outdir, "meshes")
    comp0 = root.find("compiler")
    old_meshdir = (comp0.get("meshdir") if comp0 is not None else "") or ""
    n_copy = 0
    for mesh in root.iter("mesh"):
        f = mesh.get("file")
        if not f:
            continue
        if os.path.isabs(f):
            src = f
        else:
            # 依次尝试: 旧 meshdir / 模型目录
            cands = [os.path.join(old_meshdir, f), os.path.join(outdir, f)]
            src = next((c for c in cands if os.path.isfile(c)), cands[0])
        base = os.path.basename(f)
        dst = os.path.join(mdir, base)
        if os.path.isfile(src):
            os.makedirs(mdir, exist_ok=True)
            if not os.path.isfile(dst) or os.path.getsize(dst) != os.path.getsize(src):
                shutil.copy2(src, dst); n_copy += 1
        mesh.set("file", base)                      # 改成裸文件名
    comp = root.find("compiler")
    if comp is None:
        comp = ET.Element("compiler"); root.insert(0, comp)
    comp.set("meshdir", "meshes")                   # 配合裸文件名
    comp.set("strippath", "false")
    tree.write(mjcf_path, encoding="utf-8", xml_declaration=True)
    return n_copy


def add_defaults(mjcf_path, armature, damping):
    """给所有关节加 armature/damping —— 小惯量连杆链的数值稳定性关键."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(mjcf_path); root = tree.getroot()
    d = root.find("default")
    if d is None:
        d = ET.Element("default"); root.insert(0, d)
    j = d.find("joint")
    if j is None:
        j = ET.SubElement(d, "joint")
    j.set("armature", str(armature)); j.set("damping", str(damping))
    tree.write(mjcf_path, encoding="utf-8", xml_declaration=True)


def name_geoms(mjcf_path):
    """给没有名字的 geom 按其所属刚体命名 (URDF 转换后 geom 通常无名, 导致水动力规则匹配不上)."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(mjcf_path); root = tree.getroot()
    n = 0
    def walk(el, body_name):
        nonlocal n
        for ch in el:
            if ch.tag == "body":
                walk(ch, ch.get("name", body_name))
            else:
                if ch.tag == "geom" and not ch.get("name"):
                    ch.set("name", f"{body_name}_g{n}"); n += 1
                walk(ch, body_name)
    wb = root.find("worldbody")
    if wb is not None:
        walk(wb, "world")
    tree.write(mjcf_path, encoding="utf-8", xml_declaration=True)
    return n


def report(model, rules_hint=("2.4", "2.3", "2.2", "2.1", "1.", "base")):
    print("\n" + "="*64)
    print("模型报告")
    print("="*64)
    print(f"刚体 {model.nbody-1} | 几何体 {model.ngeom} | 自由度 {model.nq} | 执行器 {model.nu}")
    print("\n[执行器] (步态参数会按这些自动生成)")
    for i in range(model.nu):
        print("   ", mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i))
    print("\n[几何体 → 水动力系数分配] (按名字关键字匹配 config.json 的 hydro.rules)")
    hit = {k: [] for k in rules_hint}; miss = []
    for g in range(model.ngeom):
        n = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, g) or f"geom{g}"
        for k in rules_hint:
            if k in n.lower(): hit[k].append(n); break
        else: miss.append(n)
    for k, v in hit.items():
        if v: print(f"    含 '{k}' → {len(v)} 个: {', '.join(v[:4])}{' ...' if len(v)>4 else ''}")
    if miss:
        print(f"    ⚠ 未匹配(将用 default 系数) → {len(miss)} 个: {', '.join(miss[:6])}{' ...' if len(miss)>6 else ''}")
        print("      建议: 在 config.json 的 hydro.rules 里加规则, 或给 geom 改名")
    tot_m = float(np.sum(model.body_mass[1:]))
    print(f"\n总质量 {tot_m*1000:.1f} g   (若与真机不符, 需在模型里修正 mass)")
    print("="*64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("urdf"); ap.add_argument("out")
    ap.add_argument("--drive", default="")
    ap.add_argument("--armature", type=float, default=1e-4, help="关节等效转子惯量(稳定小连杆链的关键)")
    ap.add_argument("--damping", type=float, default=0.01, help="关节阻尼")
    ap.add_argument("--kp", type=float, default=20.0, help="位置执行器刚度")
    ap.add_argument("--force", type=float, default=5.0, help="执行器最大力矩 N·m")
    ap.add_argument("--meshdir", default="")
    a = ap.parse_args()

    urdf_dir = os.path.dirname(os.path.abspath(a.urdf))
    raw = open(a.urdf, encoding="utf-8", errors="ignore").read()
    raw_nopkg = re.sub(r'package://[^/]+/', '', raw)
    meshdir = a.meshdir or find_meshdir(a.urdf, raw_nopkg)
    print(f"[1/4] 读取 {a.urdf}\n      网格目录(自动定位): {meshdir}")
    fixed = fix_urdf(a.urdf, meshdir)
    outdir = os.path.dirname(os.path.abspath(a.out)) or "."
    os.makedirs(outdir, exist_ok=True)
    tmp = os.path.join(outdir, "_mj_tmp.urdf")
    open(tmp, "w", encoding="utf-8").write(fixed)

    print("[2/4] MuJoCo 编译 URDF ...")
    try:
        m = mujoco.MjModel.from_xml_path(tmp)
    except Exception as e:
        print("  [X] 编译失败:", str(e)[:400])
        print("  常见原因: 网格路径不对(用 --meshdir 指定) / STL 文件缺失 / 惯量非法")
        try: os.remove(tmp)
        except Exception: pass
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    mujoco.mj_saveLastXML(a.out, m)
    try: os.remove(tmp)
    except Exception: pass
    print(f"[3/4] 已存 MJCF: {a.out}")

    nc = make_portable(a.out)
    print(f"      网格已复制 {nc} 个到 robots/meshes/ (模型改用相对路径, 可换电脑)")
    add_defaults(a.out, a.armature, a.damping)
    nn = name_geoms(a.out)
    print(f"      已为 {nn} 个几何体按刚体名命名 (水动力规则才能匹配)")

    print("[4/4] 添加执行器 ...")
    keys = [k for k in a.drive.split(",") if k.strip()]
    _, sel = add_actuators(a.out, keys, kp=a.kp, forcerange=a.force)
    print(f"      添加了 {len(sel)} 个位置执行器" + (f" (匹配 {keys})" if keys else " (所有可动关节)"))
    m2 = mujoco.MjModel.from_xml_path(a.out)
    report(m2)
    print(f"\n下一步: 把 config.json 里的 \"model\" 改成 \"{a.out}\", 然后双击 1_看示例步态.bat")


if __name__ == "__main__":
    main()
