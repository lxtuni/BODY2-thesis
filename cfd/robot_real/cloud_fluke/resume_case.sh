#!/usr/bin/env bash
# 停掉正在跑的算例 → 改数值设置 → 从最近写出的时刻续跑（startFrom latestTime，已算的不浪费）→ 跑完自动重构
#   bash resume_case.sh <算例目录> [--tracking] [--limitU 1.2] [--maxCo 1.0] [--maxAlphaCo 1.0] [--np 8] [--no-restart]
#   --tracking   overset 洞切改 trackingInverseDistance（v2006+；增量更新供体，mesh.update 快数倍；若第一步就报错或洞数不对，再跑一次不带它即可切回）
#   --no-restart 只改设置不启动（例如想手动看一眼）
#   --mode M     洞切设置整体换成 run_fluke.sh 的同名模式：id / id1mm / track2（tracking，每区一组体素，推荐）/ track / cvw
#                （比 --tracking 更稳：会按 tracking 要求的格式重写 searchBoxDivisions）
# 改动：system/fvOptions limitVelocity（压住伪速度）；controlDict maxCo / maxAlphaCo；fvSolution p_rgh 容差 1e-7/relTol 0.05；controlDict 加 maxU 监控（若没有）
set -e
CASE=""; TRACK=""; MODE=""; LIMU=1.2; MAXCO=1.0; MAXACO=1.0; NP=""; START=1
while [ $# -gt 0 ]; do case "$1" in
    --tracking) TRACK=trackingInverseDistance ;; --inverse) TRACK=inverseDistance ;; --mode) MODE="$2"; shift ;;
    --limitU) LIMU="$2"; shift ;; --maxCo) MAXCO="$2"; shift ;; --maxAlphaCo) MAXACO="$2"; shift ;; --np) NP="$2"; shift ;;
    --no-restart) START=0 ;; *) CASE="$1" ;; esac; shift; done
CASE="${CASE:?用法: bash resume_case.sh <算例目录> [--tracking] ...}"; CASE="${CASE/#\~/$HOME}"; CASE="$(cd "$CASE" && pwd)"
[ -f "$CASE/system/controlDict" ] || { echo "不是算例目录：$CASE"; exit 1; }
[ -n "$NP" ] || NP=$(grep -a -m1 numberOfSubdomains "$CASE/system/decomposeParDict" | tr -dc '0-9')
APP=$(grep -a -m1 "^application" "$CASE/system/controlDict" | awk '{print $2}' | tr -d ';'); LOG="$CASE/log.$APP"

# 1. 停掉这个目录里在跑的 Allrun / mpirun / 求解器（按进程的工作目录认，不误伤别的算例）
stop_case() {
    local n=0
    for p in $(pgrep -f "Allrun|mpirun|prterun|prted|$APP|redistributePar" 2>/dev/null); do
        [ "$p" = "$$" ] || [ "$p" = "$PPID" ] && continue
        [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$CASE" ] && { kill "$p" 2>/dev/null && n=$((n+1)); }
    done
    [ $n -gt 0 ] && { echo "已停止 $n 个进程，等 5 s ..."; sleep 5; }
    for p in $(pgrep -f "$APP|redistributePar" 2>/dev/null); do
        [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$CASE" ] && kill -9 "$p" 2>/dev/null || true
    done
    return 0
}
stop_case

# 2. 改设置（python 就地改，改前各留一份 .bak）
python3 - "$CASE" "$LIMU" "$MAXCO" "$MAXACO" "$TRACK" "$MODE" <<'EOF'
import re, sys, os, shutil
case, limU, maxCo, maxACo, track, mode = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
def edit(rel, fn):
    p = os.path.join(case, rel); s = open(p).read(); s2 = fn(s)
    if s2 != s:
        if not os.path.exists(p + ".bak"): shutil.copyfile(p, p + ".bak")
        open(p, "w").write(s2); print(f"  改了 {rel}")
def ctrl(s):
    s = re.sub(r"^maxCo\s+[\d.e-]+;", f"maxCo           {maxCo};", s, flags=re.M)
    if re.search(r"^maxAlphaCo", s, flags=re.M): s = re.sub(r"^maxAlphaCo\s+[\d.e-]+;", f"maxAlphaCo      {maxACo};", s, flags=re.M)
    if "fieldMinMax" not in s:
        blk = '    maxU { type fieldMinMax; libs ("libfieldFunctionObjects.so"); mode magnitude; fields (U); location yes; writeControl timeStep; writeInterval 1; log true; }\n'
        i = s.rstrip().rfind("}"); s = s[:i] + blk + s[i:]
    return s
edit("system/controlDict", ctrl)
fvo = os.path.join(case, "system", "fvOptions")
if limU > 0:
    txt = f"""FoamFile {{ version 2.0; format ascii; class dictionary; location "system"; object fvOptions; }}
// 速度上限：桨尖 ~0.4 m/s，水里真实速度 < 1 m/s；空气区 / 桨面薄片单元里的伪速度可到几 m/s，会把 Δt 卡在 0.1 ms
limitU {{ type limitVelocity; active yes; selectionMode all; max {limU}; }}
"""
    if not os.path.exists(fvo) or open(fvo).read() != txt: open(fvo, "w").write(txt); print(f"  写了 system/fvOptions（limitVelocity {limU} m/s）")
elif os.path.exists(fvo): os.rename(fvo, fvo + ".off"); print("  关掉 fvOptions")
edit("system/fvSolution", lambda s: re.sub(r"p_rgh\s*\{\s*solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; \}",
                                          "p_rgh       { solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.05; }", s))
if track: edit("system/fvSchemes", lambda s: re.sub(r"method\s+\w+;", f"method              {track};", s, count=1))
def set_mode(s):                       # 与 run_fluke.sh 的 set_overset 相同
    m = re.search(r"oversetInterpolation\s*\{(.*?)\n\}", s, re.S); blk = m.group(1)
    box = re.search(r"searchBox\s+(\([^)]*\)\s*\([^)]*\));", blk).group(1)
    o = re.search(r"//\s*ORIG_DIV\s*\(([\d ]+)\)", s) or re.search(r"searchBoxDivisions\s*\(\s*([\d ]+)\)", blk)
    n0 = [int(x) for x in o.group(1).split()]
    if mode == "id1mm":
        k = 120.0 / max(n0); n = [max(40, int(round(x * k))) for x in n0] if max(n0) > 120 else n0
        body = f"    method              inverseDistance;\n    searchBox           {box};\n    searchBoxDivisions  ({n[0]} {n[1]} {n[2]});"
    elif mode == "id":
        body = f"    method              inverseDistance;\n    searchBox           {box};\n    searchBoxDivisions  ({n0[0]} {n0[1]} {n0[2]});"
    elif mode == "track2":
        v = f"({n0[0]} {n0[1]} {n0[2]})"
        body = f"    method              trackingInverseDistance;\n    searchBox           {box};\n    searchBoxDivisions  ({v} {v});"
    elif mode == "track":
        body = f"    method              trackingInverseDistance;\n    searchBox           {box};\n    searchBoxDivisions  ({n0[0]} {n0[1]} {n0[2]});"
    elif mode == "cvw":
        body = "    method              cellVolumeWeight;"
    else:
        sys.exit("未知 --mode " + mode)
    s = s[:m.start(1)] + "\n" + body + s[m.end(1):]
    if "ORIG_DIV" not in s: s = s.replace("oversetInterpolation", f"// ORIG_DIV ({n0[0]} {n0[1]} {n0[2]})\noversetInterpolation", 1)
    return s
if mode:
    edit("system/fvSchemes", set_mode); open(os.path.join(case, "overset_mode.txt"), "w").write(mode + "\n")
print("  当前：maxCo", maxCo, " maxAlphaCo", maxACo, " limitU", limU, " overset", (track or re.search(r"method\s+(\w+);", open(os.path.join(case,'system/fvSchemes')).read()).group(1)))
EOF

# 3. 旧 log 改名（保留），新 log 从头写；续跑 + 跑完自动重构
[ "$START" = 1 ] || { echo "设置已改，没有启动（--no-restart）"; exit 0; }
if [ -f "$LOG" ]; then n=1; while [ -f "$LOG.part$n" ]; do n=$((n+1)); done; mv "$LOG" "$LOG.part$n"; echo "  旧 log → $(basename "$LOG").part$n"; fi
rm -f "$CASE/log.redistributePar.reconstruct"
LATEST=$(ls -d "$CASE"/processor0/[0-9]* 2>/dev/null | sed 's#.*/##' | sort -g | tail -1)
echo "  从 t = ${LATEST:-0} 续跑，$NP 核（进度：bash progress.sh $CASE）"
cd "$CASE"
nohup bash -c "mpirun -np $NP $APP -parallel > '$LOG' 2>&1; mpirun -np $NP redistributePar -reconstruct -parallel > log.redistributePar.reconstruct 2>&1; touch case.foam" > resume.out 2>&1 &
sleep 20; grep -a -m1 -A3 "calculated" "$LOG" 2>/dev/null || tail -5 "$LOG"
