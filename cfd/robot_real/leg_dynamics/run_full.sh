#!/usr/bin/env bash
# 整机一键流程：生成算例 → 几何核对 → 跑 closed → 跑 open → 受力（x/y/z）→ 四视角视频 + 四合一
#   bash run_full.sh [U] [cycles]        默认 U = 0.1 m/s、2 个周期；每个算例 8 核约 3–5 h
#   nohup bash run_full.sh 0.1 > ~/run/run_full_U01.out 2>&1 &     # 挂后台跑，之后看 ~/run/run_full_U01.out
# 可重复执行：已生成 / 已跑完的步骤自动跳过（判据：gait.json、log.overInterDyMFoam 里的 End、mp4 是否存在）。
# 算例名：~/run/full_<closed|open>_U01（U = 0.1）、U00（U = 0）、U02（U = 0.2）
set -e
U="${1:-0.1}"; CYC="${2:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"
TAG=$(python3 -c "u=$U*10; print('U%02d'%round(u) if abs(u-round(u))<1e-6 else 'U%03d'%round(u*10))")   # 0.1 → U01, 0.15 → U015
RUN=~/run; mkdir -p "$RUN"
say(){ printf "\n===== [%s] %s =====\n" "$(date '+%H:%M:%S')" "$*"; }

for fan in closed open; do
    CASE="$RUN/full_${fan}_${TAG}"
    if [ ! -f "$CASE/gait.json" ]; then
        say "生成 $fan 算例 $CASE"
        python3 "$HERE/make_full_cfd.py" --case "$CASE" --fan "$fan" --U "$U" --cycles "$CYC" $GEN_ARGS     # GEN_ARGS="--overset trackingInverseDistance" 等附加参数
        python3 "$HERE/full_geometry_check.py" "$CASE" || true         # 几何核对图 → results_full/
    fi
    LOG="$CASE/log.overInterDyMFoam"
    busy(){ [ -n "$(find "$CASE" -maxdepth 1 -name 'log.*' -mmin -3 2>/dev/null)" ]; }   # 3 分钟内有 log 在写 = 有别的进程在跑这个算例
    if busy && ! { [ -f "$LOG" ] && grep -q "^End" "$LOG"; }; then
        say "$CASE 里有正在跑的 Allrun（手动启动的？），等它结束，不重复启动"
        while busy && ! { [ -f "$LOG" ] && grep -q "^End" "$LOG"; }; do sleep 120; done
    fi
    if [ -f "$LOG" ] && grep -q "^End" "$LOG"; then
        while busy; do echo "    求解已 End，等重构 (redistributePar -reconstruct) 写完 ..."; sleep 60; done
        say "$fan 已跑完，跳过"
    else
        if [ -f "$LOG" ]; then           # 上次求解中断：Allrun 会跳过已有 log 的步骤，所以删掉求解/重构 log → 从最近写出的时刻续跑（controlDict 是 latestTime）
            say "$fan 上次求解没跑完，从最近写出的时刻续跑"
            rm -f "$LOG" "$CASE/log.redistributePar.reconstruct"
        elif ls "$CASE"/log.* >/dev/null 2>&1; then     # 网格阶段就断了：清掉从头做（网格 15–20 min）
            say "$fan 网格阶段没做完，./Allclean 后从头做"
            ( cd "$CASE" && ./Allclean >/dev/null 2>&1 ) || true
        fi
        say "跑 $fan（网格 + 求解 + 重构；进度：awk 那条命令看 $LOG）"
        ( cd "$CASE" && ./Allrun > allrun.out 2>&1 ) || true
        echo "--- 网格：$(grep -a -m1 'cells:' "$CASE/log.checkMesh.merged" 2>/dev/null || echo '没有 log.checkMesh.merged')"
        echo "--- 洞切（第一步）："; grep -a -m1 -A3 "calculated" "$LOG" 2>/dev/null | sed 's/^/    /' || true
        grep -q "^End" "$LOG" 2>/dev/null || { echo "!! $fan 求解没有正常结束，看 $LOG 末尾（tail -50）和 $CASE/allrun.out"; exit 1; }
        echo "--- 步数：$(grep -ac '^Time = ' "$LOG")，末时刻：$(grep -a '^Time = ' "$LOG" | tail -1)"
    fi
done

say "受力后处理（两态拼接，含连杆条带法阻力，x/y/z 三向力与力矩）"
python3 "$HERE/full_forces.py" --open "$RUN/full_open_${TAG}" --closed "$RUN/full_closed_${TAG}"

for fan in closed open; do
    CASE="$RUN/full_${fan}_${TAG}"
    if [ -f "$CASE/full_quad.mp4" ]; then say "$fan 视频已有，跳过"; continue; fi
    say "$fan 视频：3d / side / top / rear + 四合一"
    bash "$HERE/full_video.sh" "$CASE" quad
done
say "全部完成。结果在 /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results_full/"
