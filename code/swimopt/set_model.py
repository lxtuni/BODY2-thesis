# -*- coding: utf-8 -*-
"""
切换当前使用的模型 —— 自动改 config.json 里的 "model" 和 "trunk_body".
用法: python set_model.py            (列出 robots/ 里的模型让你选)
      python set_model.py robots/body2.xml   (直接指定)
"""
import os, re, sys, glob
import mujoco


def detect_trunk(model):
    """找机身刚体: 优先带自由关节的那个, 否则取第 1 个刚体."""
    for b in range(1, model.nbody):
        jadr, jnum = model.body_jntadr[b], model.body_jntnum[b]
        for j in range(jadr, jadr+jnum):
            if j >= 0 and model.jnt_type[j] == mujoco.mjtJoint.mjJNT_FREE:
                return mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, b)
    return mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, 1)


def set_in_config(cfg_path, model_path, trunk):
    s = open(cfg_path, encoding="utf-8").read()
    # 只替换第一个未被注释掉的 "model": "..."
    name = os.path.basename(model_path)
    s = re.sub(r'(?m)^(\s*)"model"\s*:\s*"[^"]*"\s*,?.*$',
               rf'\1"model": "{model_path}",     // 当前模型: {name}', s, count=1)
    s = re.sub(r'(?m)^(\s*)"trunk_body"\s*:\s*"[^"]*"\s*,?.*$',
               rf'\1"trunk_body":   "{trunk}",  // 机身刚体(自动识别)', s, count=1)
    open(cfg_path, "w", encoding="utf-8").write(s)


def main():
    cfg_path = "config.json"
    choice = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1].endswith(".xml") else None
    if choice is None:
        files = sorted(glob.glob("robots/*.xml"))
        if not files:
            print("robots/ 里没有 .xml 模型。先用 4_导入我的模型.bat 转换 URDF。"); return
        cur = ""
        try:
            m = re.search(r'(?m)^\s*"model"\s*:\s*"([^"]*)"', open(cfg_path, encoding="utf-8").read())
            cur = m.group(1) if m else ""
        except Exception: pass
        print("robots/ 里可用的模型:\n")
        for i, f in enumerate(files, 1):
            mark = "  ← 当前使用" if f.replace("\\", "/") == cur.replace("\\", "/") else ""
            print(f"  [{i}] {f}{mark}")
        print()
        k = input("选一个 (输入编号后回车): ").strip()
        if not k.isdigit() or not (1 <= int(k) <= len(files)):
            print("没选有效编号, 退出。"); return
        choice = files[int(k)-1]

    choice = choice.replace("\\", "/")
    print(f"\n加载 {choice} ...")
    model = mujoco.MjModel.from_xml_path(choice)
    trunk = detect_trunk(model)
    set_in_config(cfg_path, choice, trunk)
    print(f"已写入 config.json:")
    print(f'   "model"      : "{choice}"')
    print(f'   "trunk_body" : "{trunk}"   (自动识别的机身刚体)')
    print(f"\n模型信息: 刚体 {model.nbody-1} | 几何体 {model.ngeom} | 执行器 {model.nu}")
    if model.nu == 0:
        print("⚠ 这个模型没有执行器! 用 4_导入我的模型.bat 重新导入(会自动添加)。")
    else:
        print("执行器:", ", ".join(mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
                                   for i in range(min(model.nu, 12))) + (" ..." if model.nu > 12 else ""))
    # 提醒水动力规则是否匹配
    import json
    from simulate import load_cfg
    rules = [r["match"] for r in load_cfg(cfg_path).get("hydro", {}).get("rules", [])]
    names = [(mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, g) or "") for g in range(model.ngeom)]
    unmatched = [n for n in names if n and not any(r.lower() in n.lower() for r in rules)]
    if unmatched:
        print(f"\n⚠ 有 {len(unmatched)}/{len(names)} 个几何体没匹配上水动力规则(会用 default 系数):")
        print("   ", ", ".join(unmatched[:6]) + (" ..." if len(unmatched) > 6 else ""))
        print("    → 打开 config.json 的 hydro.rules, 把 match 改成你模型里的名字关键字")
    else:
        print("\n✓ 所有几何体都能匹配到水动力规则")
    print("\n现在双击 1_看示例步态.bat 或 5_边训练边看.bat 即可。")


if __name__ == "__main__":
    main()
