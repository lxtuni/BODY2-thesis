#!/usr/bin/env bash
# 船体形状比选：三个裸船体（M1 现船体 / M2 圆滑底 / M3 流线型圆舭）× 两个航速的阻力算例，一键生成 + 运行 + 汇总
#   cd /mnt/c/Users/L/Desktop/openfoam/robot_real/hull_study
#   nohup bash run_hulls.sh > ~/run/run_hulls.out 2>&1 &        # 挂后台；tail -f ~/run/run_hulls.out 看进度
# 每个算例：setup_robot.sh（DTCHull 模板；与静止浮态报告相同设置：FULL=1 NLEVELS=5 SNAPLEVEL=1 LAYERS=0 COARSE=1 LTS 1500 步，8 核）
#           → ./Allrun（网格 + interFoam + 重构）→ resistance.py；最后 hull_compare.py 汇总到 hull_study/results/
# 可重复执行：已跑完（log.interFoam 有 End）的算例自动跳过。
# 环境变量：SPEEDS="0.2 0.1"  HULLS="M1_current M2_roundbottom M3_streamlined"  ENDTIME=1500  PAUSE_OVERSET=1（暂停正在跑的整机 overset 算例，跑完再恢复）
#           T=-0.010（水线高度，可改；如 M2 同排水量校核 T=-0.0078）  SUFFIX=""（算例名后缀，如 _eqdisp，避免覆盖已跑完的算例，也不进入三船体汇总）
# 例：只补跑 M2 在同排水量水线（下沉 2.2 mm）的 U=0.2 算例：
#   T=-0.0078 SUFFIX=_eqdisp SPEEDS=0.2 HULLS=M2_roundbottom nohup bash run_hulls.sh > ~/run/run_hulls_eqdisp.out 2>&1 &
# 每个算例跑完后把 force.dat 与网格/用时摘要拷到 hull_study/results/（force_<算例名>.dat、run_summary.txt），Windows 侧可直接读。
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
OF=/mnt/c/Users/L/Desktop/openfoam
SETUP="$OF/robot_simple/setup_robot.sh"
SPEEDS="${SPEEDS:-0.2 0.1}"; HULLS="${HULLS:-M1_current M2_roundbottom M3_streamlined}"
T="${T:--0.010}"; SUFFIX="${SUFFIX:-}"; ENDTIME="${ENDTIME:-1500}"; PAUSE_OVERSET="${PAUSE_OVERSET:-1}"
RES="$HERE/results"; mkdir -p "$RES"
RUN=~/run; mkdir -p "$RUN"
say(){ printf "\n===== [%s] %s =====\n" "$(date '+%H:%M:%S')" "$*"; }
if [ -z "${WM_PROJECT_DIR:-}" ]; then set +u; source /usr/lib/openfoam/openfoam2606/etc/bashrc 2>/dev/null; set -u; fi
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "!! OpenFOAM 环境没加载：先 source /usr/lib/openfoam/openfoam2606/etc/bashrc"; exit 1; }
[ -f "$SETUP" ] || { echo "!! 找不到 $SETUP"; exit 1; }
[ -f "$HERE/hull_M1_current.stl" ] || cp "$HERE/../hull_body_only.stl" "$HERE/hull_M1_current.stl"      # M1 = 原船体
for H in $HULLS; do [ -f "$HERE/hull_$H.stl" ] || { echo "!! 找不到 $HERE/hull_$H.stl"; exit 1; }; done

# 让正在跑的整机 overset 算例先停一停（SIGSTOP），本脚本结束时恢复（SIGCONT）——否则两个 8 核任务抢 CPU，都慢一倍
PAUSED=""
if [ "$PAUSE_OVERSET" = "1" ]; then
    PAUSED=$(pgrep -f "overInterDyMFoa[m] -parallel" | tr '\n' ' ')     # 进程名 comm 只有 15 字符，pgrep -x 匹配不到 16 字符的 overInterDyMFoam，须用 -f 匹配整条命令行
    if [ -n "$PAUSED" ]; then say "暂停整机算例进程 $PAUSED（跑完船体算例后自动恢复）"; kill -STOP $PAUSED; fi
fi
resume(){ if [ -n "$PAUSED" ]; then say "恢复整机算例进程 $PAUSED"; kill -CONT $PAUSED 2>/dev/null || true; fi; }
trap resume EXIT

for U in $SPEEDS; do
  for H in $HULLS; do
    TAG=$(python3 -c "print('U%02d' % round($U*10))")
    CASE="$RUN/hull_${H}_${TAG}${SUFFIX}"; LOG="$CASE/log.interFoam"; NAME="$(basename "$CASE")"
    if [ -f "$LOG" ] && grep -q "^End" "$LOG"; then say "$H  U=$U 已跑完，跳过"; continue; fi
    say "生成 $H  U=$U  水线 T=$T → $CASE"
    rm -rf "$CASE"
    RUNDIR="$CASE" GEOM="$HERE/hull_$H.stl" T=$T U=$U FULL=1 NLEVELS=5 SNAPLEVEL=1 LAYERS=0 COARSE=1 PARMESH=0 ENDTIME="$ENDTIME" NFRAMES=30 \
        bash "$SETUP" < /dev/null > "$CASE.setup.out" 2>&1 || { echo "!! setup 失败，看 $CASE.setup.out"; tail -15 "$CASE.setup.out"; continue; }
    grep -E "blockMesh 通过|closed|Surface is" "$CASE.setup.out" | sed 's/^/    /'
    say "跑 $CASE（网格 ~10 min，求解 $ENDTIME 步 LTS ~20 min）"
    ( cd "$CASE" && ./Allrun > allrun.out 2>&1 )
    if ! grep -q "^End" "$LOG" 2>/dev/null; then echo "!! $H U=$U 没有正常结束：看 $LOG 末尾与 $CASE/allrun.out"; tail -5 "$CASE/allrun.out"; continue; fi
    echo "    网格：$(grep -a 'Snapped mesh' "$CASE/log.snappyHexMesh" 2>/dev/null | tail -1)"
    python3 "$OF/resistance.py" "$CASE" > "$CASE/resistance.out" 2>&1; tail -12 "$CASE/resistance.out" | sed 's/^/    /'
    echo "    用时：$(grep -a ClockTime "$LOG" | tail -1)"
    # 把结果拷到 Windows 可见的 results/：force.dat + 网格/用时摘要
    FD=$(ls -t "$CASE"/postProcessing/*/*/force.dat 2>/dev/null | head -1); [ -n "$FD" ] && cp "$FD" "$RES/force_$NAME.dat"
    { echo "== $NAME  U=$U  T=$T  $(date '+%F %H:%M')"; grep -a 'Snapped mesh' "$CASE/log.snappyHexMesh" 2>/dev/null | tail -1
      grep -a 'Layer mesh' "$CASE/log.snappyHexMesh" 2>/dev/null | tail -1; grep -a ClockTime "$LOG" | tail -1; tail -6 "$CASE/resistance.out"; echo; } >> "$RES/run_summary.txt"
  done
done

say "汇总三个船体的阻力 → $HERE/results/"
python3 "$HERE/hull_compare.py" --run "$RUN" --out "$HERE/results" --hydro "$HERE/hull_hydrostatics.json"
say "全部完成：结果在 $HERE/results/（hull_resistance.md 表、fig_hull_resistance.png、fig_hull_convergence.png、各算例 force_*.dat）"
