#!/usr/bin/env bash
# 看 overset 算例进度（整机 overInterDyMFoam 或单腿 overPimpleDyMFoam；网格阶段 / 求解阶段都能看）
#   WSL 里：   bash progress.sh [算例目录] [-w]        默认 ~/run/full_closed_U01；-w = 每 60 s 刷新一次（Ctrl+C 退出）
#              bash progress.sh ~/run/leg_V1_open
#   PowerShell：wsl bash /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/progress.sh
WATCH=0; CASE=""
for x in "$@"; do case "$x" in -w) WATCH=1 ;; *) CASE="$x" ;; esac; done
CASE="${CASE:-$HOME/run/full_closed_U01}"; CASE="${CASE/#\~/$HOME}"
show() {
    [ -d "$CASE" ] || { echo "没有 $CASE（算例还没生成？）"; return; }
    local LOG; LOG=$(ls "$CASE"/log.overInterDyMFoam "$CASE"/log.overPimpleDyMFoam 2>/dev/null | head -1); [ -n "$LOG" ] || LOG="$CASE/log.overInterDyMFoam"
    if [ ! -f "$LOG" ]; then
        echo "[$CASE] 还在网格 / 准备阶段（求解器还没开始写 log）"
        local last; last=$(ls -t "$CASE"/log.* "$CASE"/paddleMesh_*/log.* 2>/dev/null | head -1)
        [ -n "$last" ] && echo "  正在写：${last#$CASE/}   $(date -r "$last" '+%H:%M:%S')"
        [ -f "$CASE/allrun.out" ] && { echo "  allrun.out 末尾："; tail -3 "$CASE/allrun.out" | sed 's/^/    /'; } || echo "  还没有 allrun.out（Allrun 没启动？）"
        return
    fi
    local END; END=$(python3 -c "import json;g=json.load(open('$CASE/gait.json'));print(g['cycles']*g['T'])" 2>/dev/null || echo 2.5)
    awk -v E="$END" '/^Time = /{t=$3; n++; T[n]=t} /ClockTime =/{c=$7; C[n]=c} /^Courant Number/{co=$NF} /^deltaT = /{dt=$NF}
        END{ if(n==0){print "求解器刚启动，还没走时间步"; exit}
             printf "t = %.4f / %.2f s  (%.1f%%)   步数 %d   用时 %.1f 分   Δt = %s   Co_max = %s\n", t, E, 100*t/E, n, c/60, dt, co
             rem1 = (E-t)*c/t/60; k = (n>200) ? n-200 : 1
             if (n>1 && C[n]>C[k]) { rate=(T[n]-T[k])/(C[n]-C[k]); rem2=(E-t)/rate/60; printf "  约剩：按平均速度 %.0f 分 (%.1f h)，按最近 %d 步的速度 %.0f 分 (%.1f h)  ← 起步阶段 Δt 很小，以后者为准\n", rem1, rem1/60, n-k, rem2, rem2/60 }
             else printf "  约剩：按平均速度 %.0f 分 (%.1f h)\n", rem1, rem1/60 }' "$LOG"
    grep -a -m1 -A3 "calculated" "$LOG" | awk 'NR==1{printf "  计算单元 %s", $NF} /hole/{printf "   洞 %s\n", $NF}'
    grep -q "^End" "$LOG" && echo "  求解已结束 (End)$( [ -f "$CASE/log.redistributePar.reconstruct" ] && echo '，重构已开始/完成' )"
    grep -a -q "FOAM FATAL" "$LOG" && echo "  !! log 里有 FOAM FATAL ERROR，看 tail -40 $LOG"
    return 0
}
if [ "$WATCH" = 1 ]; then while true; do clear; date '+%H:%M:%S'; show; sleep 60; done; else show; fi
