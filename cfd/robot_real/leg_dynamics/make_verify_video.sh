#!/usr/bin/env bash
# 步态验证算例对比视频：现步态 / V1 / V2 / V3 各一块，同一视角（y=0 切面 |U| + 桨板 + 连杆骨架）、同一色标、相位锁定（每块每周期 20 帧）；
# 每块都是"两态合成"：划水相用 open 算例的帧、回收相用 closed 算例的帧（与力的合成规则一致）。
#   cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics
#   bash make_verify_video.sh                      # → results_verify/verify_video.mp4（4 个方案 2×2；只有 3 个时上下排列）+ 末帧 png
# 环境变量：SCHEMES="current V1 V2 V3"  RUN=~/run  FPS=8（每周期 20 帧 → 2.5 s/周期）  CYCLES_KEEP=2  UMAX=0.6  HOLD=1.5  FORCE=1 重渲染帧
#           现步态算例默认在 ~/run/leg_open 与 ~/run/leg_closed（没有就只出 V1–V3）
set -o pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; RUN="${RUN:-$HOME/run}"; FPS="${FPS:-8}"; CYC="${CYCLES_KEEP:-2}"; UMAX="${UMAX:-0.6}"; HOLD="${HOLD:-1.5}"
SCHEMES="${SCHEMES:-current V1 V2 V3}"; OUTDIR="$HERE/results_verify"; mkdir -p "$OUTDIR"; WORK="$RUN/verify_video"; mkdir -p "$WORK"
command -v pvbatch >/dev/null || { echo "没有 pvbatch：sudo apt install -y paraview"; exit 1; }
command -v ffmpeg  >/dev/null || { echo "没有 ffmpeg：sudo apt install -y ffmpeg"; exit 1; }
FONT=""; for f in /mnt/c/Windows/Fonts/msyh.ttc /mnt/c/Windows/Fonts/msyhbd.ttc /mnt/c/Windows/Fonts/simhei.ttf; do [ -f "$f" ] && { FONT="$f"; break; }; done
[ -n "$FONT" ] || FONT=$(fc-list :lang=zh -f '%{file}\n' 2>/dev/null | head -1); echo "字体：${FONT:-（无中文字体）}"
case_of(){ if [ "$1" = current ]; then echo "$RUN/leg_$2"; else echo "$RUN/leg_$1_$2"; fi; }
render(){   # $1 算例目录：渲染 side 视角帧（复用 paddle_video.py）
    local case="$1" fr="$1/frames_side"
    if [ "${FORCE:-0}" != "1" ] && [ -n "$(ls "$fr"/frame.*.png 2>/dev/null | head -1)" ]; then echo "    帧已存在，跳过 $(basename "$case")"; return 0; fi
    rm -rf "$fr"; mkdir -p "$fr"; echo "    渲染 $(basename "$case") …"
    if [ -n "${DISPLAY:-}" ]; then VIEW=side UMAX="$UMAX" LEG=1 pvbatch "$HERE/paddle_video.py" "$case" "$fr"
    else VIEW=side UMAX="$UMAX" LEG=1 xvfb-run -a pvbatch "$HERE/paddle_video.py" "$case" "$fr"; fi 2>&1 | grep -vE "vtkOpenFOAMReader|Unsupported|Unexpected token|^\s*$" | tail -2
    [ -n "$(ls "$fr"/frame.*.png 2>/dev/null | head -1)" ] || { echo "!! 没有生成帧：$case"; return 1; }
}
label_of(){   # 从 verify_analysis.json 取数字做标签
    python3 - "$OUTDIR/verify_analysis.json" "$1" <<'PY'
import json, sys
fn, k = sys.argv[1:]; nm = {"current": "现步态", "V1": "V1 前沿点", "V2": "V2 狗式占空比 0.35", "V3": "V3 低中心角"}.get(k, k)
try:
    r = json.load(open(fn))[k]
    s = f"{nm}   T {r['T']:.2f} s  DUTY {r['DUTY']:.2f}  ψ {r['psi_range'][1]:.0f}→{r['psi_range'][0]:.0f}°   F̄x {r['Fx_mean']:.1f} mN  F̄z {r['Fz_mean']:.1f} mN"
except Exception: s = nm
sys.stdout.write(s)
PY
}
STRIPS=(); NAMES=()
for S in $SCHEMES; do
    CO=$(case_of "$S" open); CC=$(case_of "$S" closed)
    [ -d "$CO" ] && [ -d "$CC" ] || { echo "-- $S：缺算例（$CO / $CC），跳过"; continue; }
    ls -d "$CO"/[1-9]* "$CO"/0.[0-9]* >/dev/null 2>&1 || { echo "-- $S：open 算例没有时间目录（没 reconstruct？），跳过"; continue; }
    echo "==> $S"
    render "$CO" || continue; render "$CC" || continue
    python3 "$HERE/compose_verify_frames.py" "$CO" "$CC" "$WORK/frames_$S" --cycles-keep "$CYC" || continue
    printf '%s' "$(label_of "$S")" > "$WORK/label_$S.txt"; echo "    标签：$(cat "$WORK/label_$S.txt")"
    DT=""; [ -n "$FONT" ] && DT="fontfile=$FONT:"
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$WORK/frames_$S/frame.%04d.png" \
        -vf "drawtext=${DT}textfile=$WORK/label_$S.txt:x=16:y=h-48:fontsize=26:fontcolor=0x0f2129:box=1:boxcolor=0xffffff@0.65:boxborderw=10,tpad=stop_mode=clone:stop_duration=$HOLD,scale=trunc(iw/2)*2:trunc(ih/2)*2" \
        -c:v libx264 -pix_fmt yuv420p -crf 18 "$WORK/strip_$S.mp4" && { STRIPS+=("$WORK/strip_$S.mp4"); NAMES+=("$S"); } || echo "!! $S 条带合成失败"
done
[ "${#STRIPS[@]}" -gt 0 ] || { echo "!! 没有任何条带"; exit 1; }
OUT="$OUTDIR/verify_video.mp4"; IN=(); for s in "${STRIPS[@]}"; do IN+=(-i "$s"); done
if [ "${#STRIPS[@]}" -eq 4 ]; then
    ffmpeg -y -loglevel error "${IN[@]}" -filter_complex "[0:v]scale=960:-2[a];[1:v]scale=960:-2[b];[2:v]scale=960:-2[c];[3:v]scale=960:-2[d];[a][b][c][d]xstack=inputs=4:layout=0_0|w0_0|0_h0|w0_h0[v]" -map "[v]" -c:v libx264 -pix_fmt yuv420p -crf 18 "$OUT"
elif [ "${#STRIPS[@]}" -eq 1 ]; then cp "${STRIPS[0]}" "$OUT"
else
    FC=""; i=0; for s in "${STRIPS[@]}"; do FC+="[$i:v]scale=1120:-2[s$i];"; i=$((i+1)); done; j=0; CH=""; while [ $j -lt $i ]; do CH+="[s$j]"; j=$((j+1)); done
    ffmpeg -y -loglevel error "${IN[@]}" -filter_complex "${FC}${CH}vstack=inputs=$i[v]" -map "[v]" -c:v libx264 -pix_fmt yuv420p -crf 18 "$OUT"
fi || { echo "!! 拼接失败"; exit 1; }
ffmpeg -y -loglevel error -sseof "-$(python3 -c "print(max(0.3, 1.5/$FPS))")" -i "$OUT" -frames:v 1 -update 1 "${OUT%.mp4}_last.png"
echo "完成：$OUT（${NAMES[*]}；每块每周期 20 帧、${FPS} fps，相位锁定，两态合成）   末帧：${OUT%.mp4}_last.png"
ffprobe -v error -select_streams v:0 -show_entries stream=width,height,nb_frames -of csv=p=0 "$OUT"
