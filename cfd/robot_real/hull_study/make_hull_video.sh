#!/usr/bin/env bash
# 三个船体（同一航速）的 CFD 动画：每个船体一条 1920×380 的横条（左：中纵剖面 |U| + 自由液面线 + 船体；右：自由液面波高俯视），
# 三条从上到下拼成一个视频（1920×1158，条带之间 6 px 白线）。在 WSL 里执行：
#   cd /mnt/c/Users/L/Desktop/openfoam/robot_real/hull_study
#   bash make_hull_video.sh                 # U=0.2 → results/hull_video_U02.mp4（另存最后一帧 results/hull_video_U02_last.png）
#   U=0.1 bash make_hull_video.sh           # 0.1 m/s 的三个算例
# 环境变量：RUN=~/run  HULLS="M1_current M2_roundbottom M3_streamlined"  FPS=5（30 帧 ≈ 6 s）  HOLD=2（末帧停留秒数）
#           FORCE=1 重新渲染帧（默认已有帧目录就跳过）  UMAX/BAND 透传给 make_hull_video.py
set -o pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; RUN="${RUN:-$HOME/run}"; U="${U:-0.2}"; FPS="${FPS:-5}"; HOLD="${HOLD:-2}"
HULLS="${HULLS:-M1_current M2_roundbottom M3_streamlined}"
TAG=$(python3 -c "print('U%02d' % round($U*10))"); OUT="${OUT:-$HERE/results/hull_video_$TAG.mp4}"; mkdir -p "$(dirname "$OUT")"
command -v pvbatch >/dev/null || { echo "没有 pvbatch：sudo apt install -y paraview"; exit 1; }
command -v ffmpeg  >/dev/null || { echo "没有 ffmpeg：sudo apt install -y ffmpeg"; exit 1; }
ffmpeg -hide_banner -filters 2>/dev/null | grep drawtext >/dev/null || echo "ffmpeg 没有 drawtext（libfreetype）滤镜，标签会缺失"
# 中文字体：优先 Windows 的微软雅黑，其次系统里的 CJK 字体
FONT=""
for f in /mnt/c/Windows/Fonts/msyh.ttc /mnt/c/Windows/Fonts/msyhbd.ttc /mnt/c/Windows/Fonts/simhei.ttf /mnt/c/Windows/Fonts/simsun.ttc; do [ -f "$f" ] && { FONT="$f"; break; }; done
[ -n "$FONT" ] || FONT=$(fc-match -f '%{file}' 'Noto Sans CJK SC' 2>/dev/null); [ -n "$FONT" ] && [ -f "$FONT" ] || FONT=$(fc-list :lang=zh -f '%{file}\n' 2>/dev/null | head -1)
echo "字体：${FONT:-（没找到中文字体，标签用默认字体，中文可能显示为方块）}"
NAME_M1_current="M1 现船体（平底）"; NAME_M2_roundbottom="M2 圆滑底"; NAME_M3_streamlined="M3 流线型圆舭"
label_of() {   # 船体名 + 阻力（从 results/hull_resistance.json 取）
    local h="$1" nm; nm=$(eval echo "\$NAME_$h"); [ -n "$nm" ] || nm="$h"
    python3 - "$HERE/results/hull_resistance.json" "$h" "$U" "$nm" <<'PY'
import json, sys
fn, h, U, nm = sys.argv[1:]; s = f"{nm}    U = {float(U):g} m/s"
try:
    r = json.load(open(fn))["results"][h][f"{float(U):.1f}"]
    s += f"    阻力 {r['Fx_mN']:.1f} mN（压差 {r['Px_mN']:.1f} / 摩擦 {r['Vx_mN']:.1f}）"
except Exception: pass
sys.stdout.write(s)
PY
}
render() {   # $1 算例目录  $2 视角
    local case="$1" v="$2" fr="$1/frames_$2"
    if [ "${FORCE:-0}" != "1" ] && [ -n "$(ls "$fr"/frame.*.png 2>/dev/null | head -1)" ]; then echo "    [$v] 已有帧，跳过（FORCE=1 重渲染）"; return 0; fi
    rm -rf "$fr"; mkdir -p "$fr"; echo "    [$v] 渲染帧 …"
    if [ -z "${DISPLAY:-}" ] && command -v xvfb-run >/dev/null; then
        VIEW="$v" U="$U" DECOMPOSED="$DEC" xvfb-run -a -s "-screen 0 1920x1080x24" pvbatch "$HERE/make_hull_video.py" "$case" "$fr"
    else
        VIEW="$v" U="$U" DECOMPOSED="$DEC" pvbatch "$HERE/make_hull_video.py" "$case" "$fr"
    fi 2>&1 | grep -vE "vtkOpenFOAMReader|Unsupported directive|Unexpected token|^\s*$" | tail -3
    [ -n "$(ls "$fr"/frame.*.png 2>/dev/null | head -1)" ] || { echo "!! [$v] 没有生成帧（display 报错的话：sudo apt install -y xvfb）"; return 1; }
}
STRIPS=(); TMP="$(mktemp -d)"
for H in $HULLS; do
    CASE="$RUN/hull_${H}_$TAG"
    [ -d "$CASE" ] || { echo "!! 没有算例 $CASE，跳过"; continue; }
    DEC=0; ls -d "$CASE"/[1-9]* >/dev/null 2>&1 || { [ -d "$CASE/processor0" ] && DEC=1 || { echo "!! $CASE 没有时间目录，跳过"; continue; }; }
    echo "==> $H  ($CASE)$([ "$DEC" = 1 ] && echo '  （未 reconstruct，直接读 processor* 分块结果）')"
    render "$CASE" side || continue; render "$CASE" top || continue
    ns=$(ls "$CASE"/frames_side/frame.*.png | wc -l); nt=$(ls "$CASE"/frames_top/frame.*.png | wc -l); [ "$ns" = "$nt" ] || echo "    注意：side $ns 帧 / top $nt 帧 不一致"
    printf '%s' "$(label_of "$H")" > "$TMP/label_$H.txt"; echo "    标签：$(cat "$TMP/label_$H.txt")"
    DT=""; [ -n "$FONT" ] && DT="fontfile=$FONT:"
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$CASE/frames_side/frame.%04d.png" -framerate "$FPS" -i "$CASE/frames_top/frame.%04d.png" \
        -filter_complex "[0:v][1:v]hstack=inputs=2[s];[s]drawtext=${DT}textfile=$TMP/label_$H.txt:x=16:y=36:fontsize=26:fontcolor=0x0f2129:box=1:boxcolor=0xffffff@0.6:boxborderw=8[t];[t]tpad=stop_mode=clone:stop_duration=$HOLD[u];[u]pad=w=iw:h=ih+6:x=0:y=0:color=0xffffff[v]" \
        -map "[v]" -r "$FPS" -c:v libx264 -pix_fmt yuv420p -crf 18 "$TMP/strip_$H.mp4" && STRIPS+=("$TMP/strip_$H.mp4") || echo "!! $H 条带合成失败"
done
[ "${#STRIPS[@]}" -gt 0 ] || { echo "!! 没有任何条带"; exit 1; }
if [ "${#STRIPS[@]}" -eq 1 ]; then cp "${STRIPS[0]}" "$OUT"; else
    IN=(); FC=""; i=0; for s in "${STRIPS[@]}"; do IN+=(-i "$s"); FC+="[$i:v]"; i=$((i+1)); done
    ffmpeg -y -loglevel error "${IN[@]}" -filter_complex "${FC}vstack=inputs=${#STRIPS[@]}[v]" -map "[v]" -c:v libx264 -pix_fmt yuv420p -crf 18 "$OUT" || { echo "!! 拼接失败"; exit 1; }
fi
ffmpeg -y -loglevel error -sseof "-$(python3 -c "print(max(0.3, 1.5/$FPS))")" -i "$OUT" -frames:v 1 -update 1 "${OUT%.mp4}_last.png"
rm -rf "$TMP"
echo "完成：$OUT   （最后一帧：${OUT%.mp4}_last.png）"; ffprobe -v error -select_streams v:0 -show_entries stream=width,height,nb_frames,r_frame_rate -of csv=p=0 "$OUT"
