#!/usr/bin/env bash
# 云端机器（ric-ms-7d25 等 Linux + OpenFOAM v2606）批处理：脚蹼/速度律/目标轨迹 CFD 验证清单（results/papers/04_目标轨迹与CFD验证清单.md）
#   cd ~/leg_dynamics                      # 整个 leg_dynamics 目录 scp 到云端后
#   nohup bash run_cloud.sh > ~/run/run_cloud.out 2>&1 &     # 挂后台；tail -f ~/run/run_cloud.out 看进度；bash progress.sh ~/run/leg_V3trap_open 看单个
# 可重复执行：已跑完（log 有 End）的算例自动跳过并收集；中途被杀的从 latestTime 续算；只到网格阶段的重新生成。
# 环境变量：ORDER="V3c12 V3c06 V3trap opt optL63 optL126"  U=0  CYCLES=3  NP=<核数，默认 min(nproc,16)>  RUN=~/run
# 每个算例结束后：force/moment/gait 拷到 leg_dynamics/results_verify/<方案>/（closed-only 方案自动配上 results_verify/V3 的 open），并跑一次汇总
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; RUN="${RUN:-$HOME/run}"; mkdir -p "$RUN"
ORDER="${ORDER:-V3c12 V3c06 V3trap opt optL63 optL126}"; U="${U:-0}"; CYCLES="${CYCLES:-3}"; APP=overPimpleDyMFoam
NCPU=$(nproc 2>/dev/null || echo 8); NP="${NP:-$(( NCPU < 16 ? NCPU : 16 ))}"
say(){ printf "\n===== [%s] %s =====\n" "$(date '+%m-%d %H:%M:%S')" "$*"; }
if [ -z "${WM_PROJECT_DIR:-}" ]; then set +u
  for rc in /usr/lib/openfoam/openfoam2606/etc/bashrc /opt/openfoam2606/etc/bashrc "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc" /usr/lib/openfoam/openfoam2512/etc/bashrc /usr/lib/openfoam/openfoam2506/etc/bashrc; do
    [ -f "$rc" ] && { source "$rc" 2>/dev/null; break; }; done; set -u; fi
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "!! OpenFOAM 环境没加载：先 source /usr/lib/openfoam/openfoam2606/etc/bashrc（或装：curl https://dl.openfoam.com/add-debian-repo.sh | sudo bash && sudo apt install openfoam2606-default）"; exit 1; }
command -v $APP >/dev/null || { echo "!! 找不到 $APP（OpenFOAM 版本 $WM_PROJECT_VERSION 没有 overset 求解器？）"; exit 1; }
python3 -c "import numpy" 2>/dev/null || { echo "!! 需要 python3-numpy：pip3 install numpy 或 sudo apt install python3-numpy"; exit 1; }
V3="T=0.81 DUTY=0.45 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5"
# 方案表：fans | 环境变量 | make_leg_cfd 额外参数
V3TRAP="T=0.72 DUTY=0.382 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5 VEL_LAW=trap RAMP=0.15 CLOSED_W=0.012"
case_of(){ case "$1" in
  V3c12)   echo "closed|$V3 CLOSED_W=0.012|" ;;                                              # C6 真实收拢 12 mm（配 V3 的 open）
  V3c06)   echo "closed|$V3 CLOSED_W=0.006|" ;;                                              # C6 羽化近似：迎流 5.6 mm
  V3trap)  echo "open closed|T=0.72 DUTY=0.382 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5 VEL_LAW=trap RAMP=0.15 CLOSED_W=0.012|" ;;   # C1 梯形律：划水 0.275 s 峰值 214 °/s、回收同 V3 0.445 s
  opt)     echo "open closed|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv" ;;              # 结构无关目标轨迹（L 84 mm）
  optL63)  echo "open closed|CLOSED_W=0.012|--traj $HERE/results/papers/traj_L63/optimal_trajectory.csv" ;;     # C2 行程 1.5 弦
  optL126) echo "open closed|CLOSED_W=0.012|--traj $HERE/results/papers/traj_L126/optimal_trajectory.csv" ;;    # C2 行程 3 弦
  # ---- 第二套：刚性矩形板 + 脚踝羽化（绕桨轴转 90°，回程边缘迎流），一次算例含整周期 ----
  D2a)     echo "plate|$V3TRAP|--feather 0.33,-0.06,0.04" ;;                                                       # D2-0 V3trap 圆弧运动学 + 等面积板 44 mm，40 ms 羽化 → 与 V3trap 两态合成比
  D2b)     echo "plate|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv --feather 0.351,-0.066,0.04" ;;   # D2-1 结构无关直线轨迹（84 mm 梯形 + 抬升回程）
  D2s)     echo "plate|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv --feather 0.351,-0.066,0.11,lin" ;;   # D2-2 实际舵机 7465W @8.4V：0.07 s/60° 带载≈0.075 → 90° 0.11 s 匀速
  D2c100)  echo "plate|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv --feather 0.351,-0.066,0.10" ;;   # D2-2 羽化 100 ms（900 °/s）
  D2c200)  echo "plate|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv --feather 0.351,-0.066,0.20" ;;   # D2-2 羽化 200 ms（450 °/s）
  D2w73)   echo "plate|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv --feather 0.351,-0.066,0.04 --plate-w 0.073" ;;  # D2-3 加宽到 73 mm（30.7 cm²）
  D2nof)   echo "plate|CLOSED_W=0.012|--traj $HERE/results/papers/optimal_trajectory.csv" ;;                          # D2-4 无脚踝（板全程正对）：羽化关节的全部价值
  *) echo "" ;; esac; }
pair_open(){ case "$1" in V3c12|V3c06) echo V3 ;; *) echo "" ;; esac; }   # closed-only 方案配哪个方案的 open
for V in $ORDER; do [ -n "$(case_of "$V")" ] || { echo "!! 未定义方案 $V"; exit 1; }; done
# 目标轨迹 CSV（没有就生成）
gen_traj(){ local d="$1"; shift; [ -f "$d/optimal_trajectory.csv" ] || { mkdir -p "$d"; python3 "$HERE/results/papers/optimal_trajectory.py" --out "$d" --no-anim "$@" > "$d/gen.out" 2>&1 || { echo "!! 轨迹生成失败：$d/gen.out"; tail -5 "$d/gen.out"; }; }; }
gen_traj "$HERE/results/papers"; gen_traj "$HERE/results/papers/traj_L63" --lp 0.063; gen_traj "$HERE/results/papers/traj_L126" --lp 0.126
say "开始：方案 $ORDER   U=$U  周期数 $CYCLES  $NP 核（机器 $NCPU 核，OpenFOAM $WM_PROJECT_VERSION）"
busy(){ [ -n "$(find "$1" -maxdepth 1 -name 'log.*' -mmin -3 2>/dev/null | head -1)" ]; }
for V in $ORDER; do
  SPEC="$(case_of "$V")"; FANS="${SPEC%%|*}"; REST="${SPEC#*|}"; ENVS="${REST%%|*}"; EXTRA="${REST#*|}"
  for FAN in $FANS; do
    CASE="$RUN/leg_${V}_$FAN"; LOG="$CASE/log.$APP"
    if [ -f "$LOG" ] && grep -q "^End" "$LOG"; then say "$V/$FAN 已跑完，收集"; python3 "$HERE/compare_verify.py" collect "$CASE" "$V" "$FAN"; continue; fi
    if [ -d "$CASE" ] && busy "$CASE"; then say "$V/$FAN 似乎正在运行（log 3 分钟内有更新），跳过"; continue; fi
    if [ -f "$LOG" ] && [ -d "$CASE/processor0" ]; then
        say "$V/$FAN 上次中断，从 latestTime 续算"
        ( cd "$CASE" && foamDictionary -entry startFrom -set latestTime system/controlDict >/dev/null 2>&1
          mv "$LOG" "$LOG.part$(ls "$LOG".part* 2>/dev/null | wc -l)"
          mpirun -np "$NP" $APP -parallel > "$LOG" 2>&1
          mpirun -np "$NP" redistributePar -reconstruct -parallel > log.redistributePar.reconstruct 2>&1; touch case.foam )
    else
        say "生成 $V/$FAN → $CASE   ($ENVS $EXTRA)"
        rm -rf "$CASE"
        env $ENVS python3 "$HERE/make_leg_cfd.py" --case "$CASE" --fan "$FAN" --U "$U" --cycles "$CYCLES" --np "$NP" $EXTRA > "$CASE.gen.out" 2>&1 \
            || { echo "!! 生成失败，看 $CASE.gen.out"; tail -8 "$CASE.gen.out"; continue; }
        grep -E "算例：|桨尖最大速度|部件网格盒|背景域|--traj" "$CASE.gen.out" | sed 's/^/    /'
        say "跑 $CASE（网格 + 求解，看 bash progress.sh $CASE）"
        ( cd "$CASE" && ./Allrun > allrun.out 2>&1 )
    fi
    if ! grep -q "^End" "$LOG" 2>/dev/null; then echo "!! $V/$FAN 没有正常结束：看 $LOG 末尾与 $CASE/allrun.out"; tail -5 "$CASE/allrun.out" 2>/dev/null; grep -a -m3 "FOAM FATAL" -A3 "$LOG" 2>/dev/null; continue; fi
    echo "    用时：$(grep -a ClockTime "$LOG" | tail -1)   步数：$(grep -ac '^Time = ' "$LOG")"
    python3 "$HERE/compare_verify.py" collect "$CASE" "$V" "$FAN"
  done
  P="$(pair_open "$V")"
  if [ -n "$P" ] && [ -d "$HERE/results_verify/$V" ] && [ -f "$HERE/results_verify/$P/force_open.dat" ]; then
    for f in force_open.dat moment_open.dat gait_open.json run_open.json; do [ -f "$HERE/results_verify/$P/$f" ] && cp "$HERE/results_verify/$P/$f" "$HERE/results_verify/$V/$f"; done
    echo "    $V 的 open 用 $P 的（closed-only 方案）"
  fi
  say "$V 收集完毕，汇总"; python3 "$HERE/results/papers/fold_timing_sweep.py" --extra "$V=$HERE/results_verify/$V" 2>/dev/null | grep -E "^$V|^V3 " | sed 's/^/    /'
done
say "全部完成：results_verify/<方案>/ 与 results/papers/fold_timing_sweep.json"
