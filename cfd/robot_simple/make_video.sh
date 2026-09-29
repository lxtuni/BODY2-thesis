#!/usr/bin/env bash
# 一条命令出视频：pvbatch 渲染帧 → ffmpeg 合成
#   bash make_video.sh ~/run/robot_simple               # VIEW=flow（水下总览）→ robot_flow.mp4
#   VIEW=slice bash make_video.sh ~/run/robot_simple    # slice | topuw | pressure | surface
#   VIEW=grid  bash make_video.sh ~/run/robot_simple    # 四合一：flow / slice / topuw / pressure
#   VIEW=all   bash make_video.sh ~/run/robot_simple    # 五个单独视频 + 四合一
#   SKIP=4 FPS=8 ... 其余环境变量见 make_video.py 开头
set -o pipefail
CASE="${1:-$PWD}"; CASE="$(cd "$CASE" && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
FPS="${FPS:-8}"
VIEW="${VIEW:-flow}"

command -v pvbatch >/dev/null || { echo "没有 pvbatch：先  sudo apt install -y paraview"; exit 1; }
command -v ffmpeg  >/dev/null || { echo "没有 ffmpeg：先   sudo apt install -y ffmpeg"; exit 1; }

render_one() {   # $1 = 视角
    local v="$1" frames="$CASE/frames_$1" out="$CASE/robot_$1.mp4"
    rm -rf "$frames"; mkdir -p "$frames"
    echo "==> [$v] 渲染帧"
    if [ -z "$DISPLAY" ] && command -v xvfb-run >/dev/null; then
        VIEW="$v" xvfb-run -a -s "-screen 0 1920x1080x24" pvbatch "$HERE/make_video.py" "$CASE" "$frames"
    else
        VIEW="$v" pvbatch "$HERE/make_video.py" "$CASE" "$frames"
    fi 2>&1 | grep -vE "vtkOpenFOAMReader|Unsupported directive|Unexpected token|^\s*$"
    local n; n=$(ls "$frames"/frame.*.png 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || { echo "[$v] 没有生成帧（display 相关报错的话：sudo apt install -y xvfb 再试）"; return 1; }
    echo "==> [$v] 合成 $n 帧 → $out"
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$frames/frame.%04d.png" -c:v libx264 -pix_fmt yuv420p -crf 18 "$out" && echo "完成：$out"
}

make_grid() {
    local out="$CASE/robot_grid.mp4"
    for v in flow slice topuw pressure; do [ -f "$CASE/robot_$v.mp4" ] || render_one "$v" || return 1; done
    echo "==> 四合一 → $out"
    ffmpeg -y -loglevel error -i "$CASE/robot_flow.mp4" -i "$CASE/robot_slice.mp4" -i "$CASE/robot_topuw.mp4" -i "$CASE/robot_pressure.mp4" \
        -filter_complex "[0:v]scale=960:540[a];[1:v]scale=960:540[b];[2:v]scale=960:540[c];[3:v]scale=960:540[d];[a][b]hstack[t];[c][d]hstack[u];[t][u]vstack[v]" \
        -map "[v]" -c:v libx264 -pix_fmt yuv420p -crf 18 "$out" && echo "完成：$out"
}

case "$VIEW" in
    grid) make_grid ;;
    all)  for v in flow slice topuw pressure surface; do render_one "$v" || exit 1; done; make_grid ;;
    *)    render_one "$VIEW" ;;
esac
