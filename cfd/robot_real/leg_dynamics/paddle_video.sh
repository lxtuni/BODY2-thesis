#!/usr/bin/env bash
# 单腿 overset 算例出视频：bash paddle_video.sh <case> [side|3d|both]
#   VIEW 默认 side；FPS 默认 12；输出 <case>/paddle_<view>.mp4
set -e
CASE="${1:?用法: bash paddle_video.sh <case目录> [side|3d|both]}"
WHAT="${2:-side}"; FPS="${FPS:-12}"
HERE="$(cd "$(dirname "$0")" && pwd)"
render() {
    local v="$1"; local fr="$CASE/frames_$v"
    rm -rf "$fr"; mkdir -p "$fr"
    echo "==> [$v] 渲染帧 ..."
    if [ -n "$DISPLAY" ]; then VIEW=$v pvbatch "$HERE/paddle_video.py" "$CASE" "$fr"
    else VIEW=$v xvfb-run -a pvbatch "$HERE/paddle_video.py" "$CASE" "$fr"; fi
    echo "==> [$v] 合成 mp4 ..."
    ffmpeg -y -loglevel error -framerate "$FPS" -i "$fr/frame.%04d.png" -c:v libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" "$CASE/paddle_$v.mp4"
    echo "完成：$CASE/paddle_$v.mp4"
    # 顺手拷一份到共享文件夹（Windows 能直接打开）；COPY_TO= 空 可关闭
    local dest="${COPY_TO-/mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results}"
    if [ -n "$dest" ] && [ -d "$(dirname "$dest")" ]; then
        mkdir -p "$dest"; cp "$CASE/paddle_$v.mp4" "$dest/paddle_${v}_$(basename "$CASE").mp4"
        echo "已拷到：$dest/paddle_${v}_$(basename "$CASE").mp4"
    fi
}
case "$WHAT" in
    side|3d) render "$WHAT" ;;
    both) render side; render 3d ;;
    *) echo "视角只能是 side / 3d / both"; exit 1 ;;
esac
