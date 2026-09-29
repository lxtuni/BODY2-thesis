#!/usr/bin/env bash
# 单腿步态验证批处理：V1（前沿点）/ V2（狗式占空比）/ V3（低中心角）× open / closed 共 6 个 overset 算例，顺序跑完并汇总
#   cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics
#   nohup bash run_verify.sh > ~/run/run_verify.out 2>&1 &         # 挂后台；tail -f ~/run/run_verify.out 看进度
#   bash progress.sh ~/run/leg_V1_open                              # 看单个算例进度
# 可重复执行：已跑完（log 有 End）的算例自动跳过；中途被杀的算例从 latestTime 续算；只到网格阶段的算例重新生成。
# 环境变量：ORDER="V1 V2 V3"  FANS="open closed"  U=0  CYCLES=3  NP=8  RUN=~/run
#           每个算例结束后把 force/moment/gait 拷到 leg_dynamics/results_verify/<方案>/，并更新 verify_summary.md / fig_verify.png（含现步态基准）
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; RUN="${RUN:-$HOME/run}"; mkdir -p "$RUN"
ORDER="${ORDER:-V1 V2 V3}"; FANS="${FANS:-open closed}"; U="${U:-0}"; CYCLES="${CYCLES:-3}"; NP="${NP:-8}"; APP=overPimpleDyMFoam
say(){ printf "\n===== [%s] %s =====\n" "$(date '+%m-%d %H:%M:%S')" "$*"; }
if [ -z "${WM_PROJECT_DIR:-}" ]; then set +u; source /usr/lib/openfoam/openfoam2606/etc/bashrc 2>/dev/null; set -u; fi
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "!! OpenFOAM 环境没加载：先 source /usr/lib/openfoam/openfoam2606/etc/bashrc"; exit 1; }
# 三个方案的步态参数（与 scan_gait_pareto.py / 2 页要点报告一致）
gait_of(){ case "$1" in
    V1) echo "T=0.89 DUTY=0.49 SWEEP=60 PSI_START=135 SWITCH_FRAC=0.5" ;;   # 前沿点：划水 0.44 s（速率 ×1.43）+ 回收 0.45 s 慢蜷缩，模型 47.7 mN
    V2) echo "T=1.25 DUTY=0.35 SWEEP=60 PSI_START=136 SWITCH_FRAC=0.5" ;;   # 狗式占空比 0.35：划水 0.44 s + 回收 0.81 s 慢蜷缩，模型 36.5 mN
    V3) echo "T=0.81 DUTY=0.45 SWEEP=50 PSI_START=120 SWITCH_FRAC=0.5" ;;   # 低中心角 120→70°：降桨板力的垂向分量，模型 42.6 mN、F_z/F_x 0.30
    *) echo "" ;; esac; }
for V in $ORDER; do [ -n "$(gait_of "$V")" ] || { echo "!! 未定义方案 $V"; exit 1; }; done
say "开始：方案 $ORDER × $FANS   U=$U  周期数 $CYCLES  $NP 核"
for V in $ORDER; do echo "--- $V: $(gait_of "$V")"; env $(gait_of "$V") python3 "$HERE/check_clearance.py" 2>/dev/null | grep -E "划水相|回收相|结论|髋角" | sed 's/^/    /'; done

busy(){ [ -n "$(find "$1" -maxdepth 1 -name 'log.*' -mmin -3 2>/dev/null | head -1)" ]; }
for V in $ORDER; do
  for FAN in $FANS; do
    CASE="$RUN/leg_${V}_$FAN"; LOG="$CASE/log.$APP"
    if [ -f "$LOG" ] && grep -q "^End" "$LOG"; then
        say "$V/$FAN 已跑完，收集结果"; python3 "$HERE/compare_verify.py" collect "$CASE" "$V" "$FAN"; continue
    fi
    if [ -d "$CASE" ] && busy "$CASE"; then say "$V/$FAN 似乎正在运行（log 3 分钟内有更新），跳过"; continue; fi
    if [ -f "$LOG" ] && [ -d "$CASE/processor0" ]; then
        say "$V/$FAN 上次中断，从 latestTime 续算"
        ( cd "$CASE" && foamDictionary -entry startFrom -set latestTime system/controlDict >/dev/null 2>&1
          mv "$LOG" "$LOG.part$(ls "$LOG".part* 2>/dev/null | wc -l)"
          mpirun -np "$NP" $APP -parallel > "$LOG" 2>&1
          mpirun -np "$NP" redistributePar -reconstruct -parallel > log.redistributePar.reconstruct 2>&1; touch case.foam )
    else
        say "生成 $V/$FAN → $CASE   ($(gait_of "$V"))"
        rm -rf "$CASE"
        env $(gait_of "$V") python3 "$HERE/make_leg_cfd.py" --case "$CASE" --fan "$FAN" --U "$U" --cycles "$CYCLES" --np "$NP" > "$CASE.gen.out" 2>&1 \
            || { echo "!! 生成失败，看 $CASE.gen.out"; tail -8 "$CASE.gen.out"; continue; }
        grep -E "算例：|桨尖最大速度|部件网格盒|背景域" "$CASE.gen.out" | sed 's/^/    /'
        say "跑 $CASE（网格 + 求解，看 bash progress.sh $CASE）"
        ( cd "$CASE" && ./Allrun > allrun.out 2>&1 )
    fi
    if ! grep -q "^End" "$LOG" 2>/dev/null; then echo "!! $V/$FAN 没有正常结束：看 $LOG 末尾与 $CASE/allrun.out"; tail -5 "$CASE/allrun.out" 2>/dev/null; grep -a -m3 "FOAM FATAL" -A3 "$LOG" 2>/dev/null; continue; fi
    echo "    用时：$(grep -a ClockTime "$LOG" | tail -1)   步数：$(grep -ac '^Time = ' "$LOG")"
    python3 "$HERE/compare_verify.py" collect "$CASE" "$V" "$FAN"
  done
  say "$V 两态收集完毕，更新汇总"; python3 "$HERE/compare_verify.py" analyze --baseline "$HERE/results" 2>/dev/null | grep -E "^\[|^\|" | head -12
done
say "全部完成：results_verify/verify_summary.md、fig_verify.png（含现步态基准）"
