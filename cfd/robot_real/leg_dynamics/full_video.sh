#!/usr/bin/env bash
# 整机 overset 算例出视频：bash full_video.sh <case> [3d|side|top|rear|both|all|quad]   quad = 四视角 2×2 合成
#   VIEW 默认 3d；FPS 默认 12；输出 <case>/full_<view>.mp4
set -e
CASE="${1:?用法: bash full_video.sh <case目录> [3d|side|top|rear|both|all|quad]}"
WHAT="${2:-3d}"; FPS="${FPS:-12}"
HERE="$(cd "$(dirname "$0")" && pwd)"
render() {
    local v="$1"; local fr="$CASE/frames_$v"
    rm -rf "$fr"; mkdir -p "$fr"
    echo "==> [$v] 渲染帧 ..."
    if [ -n "$DISPLAY" ]; then VIEW=$v pvbatch "$HERE/full_video.py" "$CASE" "$fr"
    else VIEW=$v xvfb-run -a pvbatch "$HERE/full_video.py" "$CASE" "$fr"; fi
    echo "==> [$v] 合成 mp4 ..."
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$fr/frame.%04d.png" -c:v libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" "$CASE/full_$v.mp4"
    echo "完成：$CASE/full_$v.mp4"
    # 顺手拷一份到共享文件夹（Windows 能直接打开）；COPY_TO= 空 可关闭
    local dest="${COPY_TO-/mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results_full}"
    if [ -n "$dest" ] && [ -d "$(dirname "$dest")" ]; then
        mkdir -p "$dest"; cp "$CASE/full_$v.mp4" "$dest/full_${v}_$(basename "$CASE").mp4"
        echo "已拷到：$dest/full_${v}_$(basename "$CASE").mp4"
    fi
}
quad() {
    # 四合一：3d | side 上排，top | rear 下排（同一算例四个视角帧数相同）
    for v in 3d side top rear; do [ -f "$CASE/full_$v.mp4" ] || render "$v"; done
    echo "==> 四合一 ..."
    ffmpeg -y -loglevel error -i "$CASE/full_3d.mp4" -i "$CASE/full_side.mp4" -i "$CASE/full_top.mp4" -i "$CASE/full_rear.mp4" \
        -filter_complex "[0:v][1:v]hstack=inputs=2[t];[2:v][3:v]hstack=inputs=2[b];[t][b]vstack=inputs=2,scale=1920:-2[v]" -map "[v]" -c:v libx264 -pix_fmt yuv420p "$CASE/full_quad.mp4"
    echo "完成：$CASE/full_quad.mp4"
    local dest="${COPY_TO-/mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results_full}"
    if [ -n "$dest" ] && [ -d "$(dirname "$dest")" ]; then mkdir -p "$dest"; cp "$CASE/full_quad.mp4" "$dest/full_quad_$(basename "$CASE").mp4"; echo "已拷到：$dest/full_quad_$(basename "$CASE").mp4"; fi
}
case "$WHAT" in
    side|3d|top|rear) render "$WHAT" ;;
    both) render 3d; render side ;;  all) render 3d; render side; render top; render rear ;;
    quad) quad ;;
    *) echo "视角只能是 3d / side / top / rear / both / all / quad（四合一）"; exit 1 ;;
esac
