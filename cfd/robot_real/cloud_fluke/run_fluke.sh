#!/usr/bin/env bash
# 升力型尾鳍单腿 CFD + 阻力型航速对照 —— 一条命令按顺序跑完（april：WSL + OpenFOAM v2606）
#   cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics
#   nohup bash run_fluke.sh > ~/run/run_fluke.out 2>&1 &        # 挂后台；tail -f ~/run/run_fluke.out
#   bash progress.sh ~/run/fluke_L0                             # 看单个算例
# 可重复执行：跑完的（log 有 End）自动跳过并收集；中途被杀的从 latestTime 续算；只到网格阶段的重新生成。
# 环境变量：WRITES=<每周期写几帧流场，默认 20>  KEEP_PROC=0（跑完删 processor*，不在 KEEP3D 里的再删三维流场）  KEEP3D="L0 L1a28"  ORDER="..."（默认下面的顺序）  NP=<进程数，默认 = 物理核数（不算超线程），最多 16>  RUN=~/run
#
# 算例（说明见 fluke/README_fluke.md）：
#   L0       尾鳍 A2，U 0.223，α 20°，运动 = Flare 单翼基准实测关节角        → 与 Flare 交叉验证（9.7 mN）
#   L2U000   同 L0 的运动，U = 0     （系泊，与阻力型 U=0 数据直接比）
#   L2U008   同 L0 的运动，U = 0.08  （Flare 自航速度）
#   L2U015   同 L0 的运动，U = 0.15
#   L1a10    攻角律 α 10°，U 0.223
#   L1a28    攻角律 α 28°，U 0.223
#   D1U008   阻力型 V3trap（open + closed），U = 0.08
#   D1U015   阻力型 V3trap，U = 0.15
#   D1U022   阻力型 V3trap，U = 0.223
#   V1       翼型基准：NACA0012 c 0.1 m、展 0.6 m，h0/c 0.75，轴 c/3，ψ 90°，St 0.25，α_max 15°，U 0.4（Schouveiler 2005）
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; RUN="${RUN:-$HOME/run}"; mkdir -p "$RUN"; FL="$HERE/fluke"
ORDER="${ORDER:-L0 L2U000 D1U008 D1U015 D1U022 L1a10 L1a28 L2U008 L2U015 V1}"; APP=overPimpleDyMFoam
NCPU=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l); [ "${NCPU:-0}" -gt 0 ] || NCPU=$(nproc 2>/dev/null || echo 8)
NP="${NP:-$(( NCPU < 16 ? NCPU : 16 ))}"          # 物理核数：超线程开 16 进程 = 超订，overset 通信反而变慢
say(){ printf "\n===== [%s] %s =====\n" "$(date '+%m-%d %H:%M:%S')" "$*"; }
if [ -z "${WM_PROJECT_DIR:-}" ]; then set +u
  for rc in /usr/lib/openfoam/openfoam2606/etc/bashrc /opt/openfoam2606/etc/bashrc "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc" /usr/lib/openfoam/openfoam2512/etc/bashrc; do
    [ -f "$rc" ] && { source "$rc" 2>/dev/null; break; }; done; set -u; fi
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "!! OpenFOAM 环境没加载：先 source /usr/lib/openfoam/openfoam2606/etc/bashrc"; exit 1; }
command -v $APP >/dev/null || { echo "!! 找不到 $APP"; exit 1; }
python3 -c "import numpy" 2>/dev/null || { echo "!! 需要 python3-numpy：sudo apt install python3-numpy"; exit 1; }   # 运动表已随 fluke/ 附上；重新生成 V1 的才需要 scipy

UB=0.2229                                          # Flare 单翼基准的拖曳速度 = 2 Hz × 39 mm / St 0.35
V3TRAP="T=0.72 DUTY=0.382 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5 VEL_LAW=trap RAMP=0.15 CLOSED_W=0.012"
FLUKE_MESH="--stl $FL/A2_blade_cfd.stl --ref-deg 0 --comp 0.0012 --bg 0.0025 --cycles 3 --ramp 0.5 --writes ${WRITES:-20}"   # 干净叶片（无销孔/T 臂、尾缘 95 % 弦截平）；起步 0.5 周期平滑
# 方案表：类型 | 运动表文件 | 运动表生成参数 | U | make_*_cfd 额外参数
case_of(){ case "$1" in
  L0)      echo "fluke|motion_L0.csv|--bench $FL/bench_LF_active.csv --U $UB --freq 2|$UB|$FLUKE_MESH" ;;
  L2U000)  echo "fluke|motion_L0.csv|--bench $FL/bench_LF_active.csv --U $UB --freq 2|0|$FLUKE_MESH" ;;
  L2U008)  echo "fluke|motion_L0.csv|--bench $FL/bench_LF_active.csv --U $UB --freq 2|0.08|$FLUKE_MESH" ;;
  L2U015)  echo "fluke|motion_L0.csv|--bench $FL/bench_LF_active.csv --U $UB --freq 2|0.15|$FLUKE_MESH" ;;
  L2U030)  echo "fluke|motion_L0.csv|--bench $FL/bench_LF_active.csv --U $UB --freq 2|0.30|$FLUKE_MESH" ;;   # 交叉点验证：升力型 U 0.30
  L0c)     echo "fluke|motion_L0.csv|--bench $FL/bench_LF_active.csv --U $UB --freq 2|$UB|${FLUKE_MESH/--comp 0.0012 --bg 0.0025/--comp 0.0016 --bg 0.0033}" ;;   # 网格无关性：L0 粗网格（×1.33）
  L1a10)   echo "fluke|motion_a10.csv|--alpha 10 --U $UB --freq 2 --nh 6|$UB|$FLUKE_MESH" ;;
  L1a28)   echo "fluke|motion_a28.csv|--alpha 28 --U $UB --freq 2 --nh 6 --nh-hinge 12|$UB|$FLUKE_MESH" ;;   # 铰链 12 阶平滑（2026-09-26）
  V1)      echo "fluke|motion_V1.csv|--heave-pitch 0.075,auto,90 --U 0.4 --St 0.25 --alpha-max 15|0.4|--naca 0012,0.1,0.6,0.3333 --origin 0,0 --comp 0.0025 --bg 0.005 --vox 0.0025 --tmin 0.009 --Lref 0.3 --cycles 3" ;;
  D1U008)  echo "drag|||0.08|" ;;
  D1U015)  echo "drag|||0.15|" ;;
  D1U022)  echo "drag|||$UB|" ;;
  D1U030)  echo "drag|||0.30|" ;;   # 交叉点验证：阻力型 U 0.30
  D0U015)  echo "cur|||0.15|" ;;    # 阻力型 真机现步态（T 1.25, DUTY 0.5, 收拢宽 30 mm）U 0.15 —— 对应真机对比试验
  D0U008)  echo "cur|||0.08|" ;;    # 阻力型 真机现步态 U 0.08 —— 确认零推力航速
  L1a20U040) echo "fluke|motion_a20U040.csv|--alpha 20 --U 0.40 --freq 2 --nh 6 --nh-hinge 12|0.40|$FLUKE_MESH" ;;   # 升力型 α 律 20°，U 0.40 > 铰链最大速度 0.29 m/s —— 升力型不受附肢速度限制
  L3U022)  echo "fluke|motion_L3U022.csv|--upload-motion_L3U022.csv-first|0.2229|$FLUKE_MESH" ;;   # 升力型 变攻角律 θ0 = γ_max/2（铰链同 L0），U 0.223 —— 对照 L0 21.4 / V3trap 35.9
  L3U030)  echo "fluke|motion_L3U030.csv|--upload-motion_L3U030.csv-first|0.30|$FLUKE_MESH" ;;     # 同上，U 0.30 —— 对照 L2U030 12.7 / V3trap 11.8（决定自航排序）
  L3U040)  echo "fluke|motion_L3U040.csv|--upload-motion_L3U040.csv-first|0.40|$FLUKE_MESH" ;;     # 同上，U 0.40 —— 变攻角律的自航交点 + 与 B1（α 律 20°）同航速对照
  *) echo "" ;; esac; }
for V in $ORDER; do [ -n "$(case_of "$V")" ] || { echo "!! 未定义算例 $V"; exit 1; }; done
say "开始：$ORDER   $NP 核（机器 $NCPU 物理核，OpenFOAM $WM_PROJECT_VERSION）"
busy(){ [ -n "$(find "$1" -maxdepth 1 -name 'log.*' -mmin -3 2>/dev/null | head -1)" ]; }

# ---- 洞切检查：第一步的 Overset analysis 里洞单元 > 3 % 总单元或 < 50 = 洞切失败（v2606 + 1.3 mm 体素实测洞 0）（v2606 上 L0 出现过 19.5 万洞）
#      → 停掉，换下一种洞切设置重启（网格不用重做）。顺序：体素 0.68 mm → 体素 ≈1.3 mm → tracking（每区一组）→ tracking → cellVolumeWeight
set_overset(){ python3 - "$1/system/fvSchemes" "$2" <<'PYEOF'
import re, sys
fn, mode = sys.argv[1], sys.argv[2]; s = open(fn).read()
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
else:
    body = "    method              cellVolumeWeight;"
s = s[:m.start(1)] + "\n" + body + s[m.end(1):]
if "ORIG_DIV" not in s: s = s.replace("oversetInterpolation", f"// ORIG_DIV ({n0[0]} {n0[1]} {n0[2]})\noversetInterpolation", 1)
open(fn, "w").write(s); print("    洞切设置 →", mode, "|", body.split(chr(10))[0].split()[-1], (body.split(chr(10))[-1].strip() if mode != "cvw" else ""))
PYEOF
}
first_holes(){  # 输出 "总单元 洞"；还没算到第一步则输出空
  awk '/Overset analysis : nCells :/{n=$NF; f=1; next} f&&/hole/{gsub(/[^0-9]/," ",$0); split($0,a," "); print n, a[1]; exit}' "$1" 2>/dev/null; }
solve_guarded(){
  local CASE="$1" LOG="$1/log.$APP" MODES="${OVERSET_MODES:-id id1mm track2 track cvw}" M PID NH R W
  for M in $MODES; do
    set_overset "$CASE" "$M" || { echo "!! 改 fvSchemes 失败"; return 1; }
    rm -rf "$CASE"/processor*/[1-9]* "$CASE"/processor*/0.[0-9]* "$CASE/postProcessing"
    ( cd "$CASE" && exec mpirun -np "$NP" $APP -parallel > "$LOG" 2>&1 ) & PID=$!
    W=0; R=""
    while kill -0 $PID 2>/dev/null && [ $W -lt 1800 ]; do
      NH=$(first_holes "$LOG"); [ -n "$NH" ] && break
      grep -a -q "FOAM FATAL" "$LOG" && break; sleep 10; W=$((W+10)); done
    NH=$(first_holes "$LOG")
    if [ -z "$NH" ]; then
      echo "    [$M] 第一步前就停了 / 报错："; { grep -a -m1 -A4 "FOAM FATAL" "$LOG" || tail -4 "$LOG"; } | sed 's/^/      /'
      kill $PID 2>/dev/null; wait $PID 2>/dev/null; mv "$LOG" "$LOG.fail_$M"; continue; fi
    set -- $NH; R=$(awk -v n=$1 -v h=$2 'BEGIN{printf "%.2f", 100*h/n}')
    echo "    [$M] 第一步：总单元 $1，洞 $2（$R %）"
    if awk -v r=$R -v h=$2 'BEGIN{exit !(r>3 || h<50)}'; then
      echo "    [$M] 洞太多或 = 0（物体内部没挖掉，背景看不到尾鳍）= 洞切失败，换下一种"; kill $PID 2>/dev/null; sleep 3; pkill -P $PID 2>/dev/null; wait $PID 2>/dev/null
      mv "$LOG" "$LOG.fail_$M"; continue; fi
    echo "$M" > "$CASE/overset_mode.txt"; wait $PID
    ( cd "$CASE" && mpirun -np "$NP" redistributePar -reconstruct -parallel > log.redistributePar.reconstruct 2>&1; touch case.foam )
    return 0
  done
  echo "!! 所有洞切设置都失败（$MODES）。把下面几行发给 Claude："; for f in "$LOG".fail_*; do echo "== $f"; grep -a -m1 -A4 "Overset analysis\|FOAM FATAL" "$f"; done
  return 1
}

# 省硬盘（KEEP_PROC=0 时）：删 processor*；不在 KEEP3D 里的算例再删三维流场（只留 0/、constant/、postProcessing/、log）。
# postProcessing 里有每步的力、尾鳍表面压力、展向中截面（每周期 40 帧），画机理图够用；三维流场只给 KEEP3D（默认 L0 L1a28）留着
KEEP3D="${KEEP3D:-L0 L1a28}"
slim(){ [ "${KEEP_PROC:-1}" = 0 ] || return 0
  rm -rf "$1"/processor*
  case " $KEEP3D " in *" $2 "*) echo "    省硬盘：删了 processor*，三维流场留着（$2 在 KEEP3D 里）" ;;
    *) find "$1" -maxdepth 1 -type d -regex '.*/[0-9.e+-]+' ! -name 0 -exec rm -rf {} + ; echo "    省硬盘：删了 processor* 和三维流场（postProcessing 的力 / 表面 / 截面都在）" ;; esac
  echo "    $(du -sh "$1" | cut -f1)   硬盘剩 $(df -h ~ | awk 'NR==2{print $4}')"; }

# 生成 / 续算 / 跳过一个算例目录。$1 目录  $2 生成命令（字符串）
run_case(){   # $3 = guard → 网格与求解分开，求解带洞切检查（尾鳍算例）；不给 = 原来的 Allrun 一条龙（阻力型）
  local CASE="$1" GEN="$2" GUARD="${3:-}" LOG="$1/log.$APP"
  if [ -f "$LOG" ] && grep -q "^End" "$LOG"; then say "$(basename "$CASE") 已跑完"; return 0; fi
  if [ -d "$CASE" ] && busy "$CASE"; then say "$(basename "$CASE") 似乎正在运行（log 3 分钟内有更新），跳过"; return 1; fi
  if [ -f "$LOG" ] && [ -d "$CASE/processor0" ]; then
    say "$(basename "$CASE") 上次中断，从 latestTime 续算"
    ( cd "$CASE" && mv "$LOG" "$LOG.part$(ls "$LOG".part* 2>/dev/null | wc -l)"
      mpirun -np "$NP" $APP -parallel > "$LOG" 2>&1
      mpirun -np "$NP" redistributePar -reconstruct -parallel > log.redistributePar.reconstruct 2>&1; touch case.foam )
  else
    say "生成 $(basename "$CASE")"; rm -rf "$CASE"
    eval "$GEN" > "$CASE.gen.out" 2>&1 || { echo "!! 生成失败，看 $CASE.gen.out"; tail -8 "$CASE.gen.out"; return 1; }
    grep -E "算例：|部件网格盒|背景域|最大速度|洞切" "$CASE.gen.out" | sed 's/^/    /'
    if [ -z "$GUARD" ]; then say "跑 $CASE（网格 + 求解；看进度 bash progress.sh $CASE）"; ( cd "$CASE" && ./Allrun > allrun.out 2>&1 ); else
    say "网格 $CASE"
    ( cd "$CASE" && MESH_ONLY=1 ./Allrun > allrun.out 2>&1 )
    [ -d "$CASE/processor0" ] || { echo "!! 网格阶段失败：看 $CASE/allrun.out"; tail -5 "$CASE/allrun.out"; return 1; }
    # 背景加密检查：v2606 上 L0 的背景网格 = 194824 = blockMesh 原样（snappy 没加密）→ 10 mm 粗单元整片被判成洞
    local NB NR; NB=$(grep -a -m1 "nCells:" "$CASE/log.blockMesh" | tr -dc '0-9'); NR=$(grep -a "Refined mesh : cells:" "$CASE/log.snappyHexMesh" | tail -1 | sed 's/.*cells:\([0-9]*\).*/\1/')
    echo "    背景网格：blockMesh ${NB:-?} → snappy 加密后 ${NR:-没有}"
    if [ -z "$NR" ] || [ "${NR:-0}" -lt $(( ${NB:-1} * 3 / 2 )) ]; then
      echo "!! 背景网格没有加密（snappyHexMesh 失败或没写回 constant/polyMesh）。把这几行发给 Claude："
      grep -a -m3 -A6 "FOAM FATAL\|FOAM Warning" "$CASE/log.snappyHexMesh"; tail -15 "$CASE/log.snappyHexMesh"; ls "$CASE"; return 1; fi
    say "求解 $CASE（看进度 bash progress.sh $CASE）"
    solve_guarded "$CASE" || return 1; fi
  fi
  if ! grep -q "^End" "$LOG" 2>/dev/null; then echo "!! $(basename "$CASE") 没有正常结束：看 $LOG 末尾与 $CASE/allrun.out"; tail -5 "$CASE/allrun.out" 2>/dev/null; grep -a -m3 "FOAM FATAL" -A3 "$LOG" 2>/dev/null; return 1; fi
  echo "    用时：$(grep -a ClockTime "$LOG" | tail -1)   步数：$(grep -ac '^Time = ' "$LOG")"
}

for V in $ORDER; do
  SPEC="$(case_of "$V")"; IFS='|' read -r KIND MOT MGEN UU EXTRA <<< "$SPEC"
  if [ "$KIND" = fluke ]; then
    [ -f "$FL/$MOT" ] || { python3 "$HERE/fluke_traj.py" $MGEN --out "$FL/$MOT" || { echo "!! 运动表生成失败 $MOT"; continue; }; }
    CASE="$RUN/fluke_$V"
    run_case "$CASE" "python3 '$HERE/make_fluke_cfd.py' --case '$CASE' --motion '$FL/$MOT' --U $UU --np $NP $EXTRA" guard \
      && python3 "$HERE/fluke_forces.py" collect "$CASE" "$V" \
      && slim "$CASE" "$V"
  else
    if [ "$KIND" = cur ]; then GENV=""; TAG=cur; else GENV="$V3TRAP"; TAG=V3trap; fi   # cur = 真机现步态（leg_dynamics 默认参数）
    for FAN in open closed; do
      CASE="$RUN/leg_${TAG}_U${UU/./}_$FAN"
      run_case "$CASE" "env $GENV python3 '$HERE/make_leg_cfd.py' --case '$CASE' --fan $FAN --U $UU --cycles 3 --writes ${WRITES:-20} --np $NP" \
        && python3 "$HERE/compare_verify.py" collect "$CASE" "${TAG}_U${UU/./}" "$FAN" \
        && slim "$CASE" "D1"
    done
  fi
  say "$V 完成，汇总"; python3 "$HERE/fluke_forces.py" summary 2>/dev/null | sed 's/^/    /'
done
say "全部完成：results_fluke/summary.json（尾鳍）与 results_verify/V3trap_U*/（阻力型）"
