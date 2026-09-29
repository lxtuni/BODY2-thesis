#!/usr/bin/env bash
# =====================================================================
#  DTCHull 船体阻力算例 —— 自动配置脚本
#
#  用法（在 Ubuntu 终端里）：
#      bash /mnt/c/Users/L/Desktop/openfoam/setup_dtchull.sh
#      NFRAMES=80 bash .../setup_dtchull.sh      # 想要更多视频帧
#
#  它做四件事：定位教程 → 复制到 ~/run/ → 按本机核数配并行 → 按帧数配输出
# =====================================================================
set -o pipefail

NFRAMES=${NFRAMES:-50}
RUNDIR="$HOME/run/dtcHull"

B=$'\033[1;36m'; G=$'\033[0;32m'; Y=$'\033[1;33m'; R=$'\033[0;31m'; N=$'\033[0m'
step(){ printf "\n${B}==> %s${N}\n" "$*"; }
ok(){   printf "${G}    OK  %s${N}\n" "$*"; }
warn(){ printf "${Y}    !!  %s${N}\n" "$*"; }
die(){  printf "\n${R}    XX  %s${N}\n\n" "$*"; exit 1; }

# ---------------------------------------------------------------- 环境
step "检查 OpenFOAM 环境"
[ -n "$WM_PROJECT_DIR" ] || die "OpenFOAM 环境没加载。先跑 source ~/.bashrc"
ok "$(foamVersion 2>/dev/null || echo "$WM_PROJECT_DIR")"

# ---------------------------------------------------------------- 定位
step "定位 DTCHull 教程算例"
SRC=$(find "$FOAM_TUTORIALS" -maxdepth 5 -type d -name 'DTCHull' 2>/dev/null | head -1)
if [ -z "$SRC" ]; then
    warn "没找到名为 DTCHull 的算例。你的 OpenFOAM 里有这些船体相关算例："
    find "$FOAM_TUTORIALS" -maxdepth 5 -type d \
        \( -iname '*hull*' -o -iname '*DTC*' -o -iname '*ship*' -o -iname '*wigley*' \) \
        2>/dev/null | sed "s|$FOAM_TUTORIALS|\$FOAM_TUTORIALS|" | sed 's/^/      /'
    die "把上面的列表发给我，我告诉你用哪个。"
fi
ok "找到：${SRC/$FOAM_TUTORIALS/\$FOAM_TUTORIALS}"

# ---------------------------------------------------------------- 复制
step "复制到 $RUNDIR"
if [ -d "$RUNDIR" ]; then
    warn "$RUNDIR 已存在"
    read -r -p "    删掉重来？[y/N] " ans
    case "$ans" in
        y|Y) rm -rf "$RUNDIR" ;;
        *)   die "已取消。想保留就改个名字，或手动删除后重跑。" ;;
    esac
fi
mkdir -p "$(dirname "$RUNDIR")"
cp -r "$SRC" "$RUNDIR" || die "复制失败"
cd "$RUNDIR" || die "进不去 $RUNDIR"
ok "算例已就位"

# ---------------------------------------------------------------- 并行
step "配置并行分解"
PHYS=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l)
LOGI=$(nproc)
[ "$PHYS" -ge 1 ] 2>/dev/null || PHYS=$LOGI
NSUB=$PHYS
[ "$NSUB" -lt 1 ] && NSUB=1
[ "$NSUB" -gt 16 ] && NSUB=16
ok "物理核 $PHYS / 逻辑核 $LOGI  →  分解为 $NSUB 块"
[ "$PHYS" -lt "$LOGI" ] && warn "有超线程。CFD 是内存带宽瓶颈，用物理核数通常比用满逻辑核更快。"

cat > system/decomposeParDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains $NSUB;
method            scotch;
EOF
ok "system/decomposeParDict 已写入（scotch 自动分区，不用手工指定 x/y/z）"

# ---------------------------------------------------------------- 输出
step "配置结果输出（为出视频准备 $NFRAMES 帧）"
ET=$(foamDictionary -entry endTime   -value system/controlDict 2>/dev/null)
ST=$(foamDictionary -entry startTime -value system/controlDict 2>/dev/null || echo 0)
[ -n "$ET" ] || die "读不到 controlDict 的 endTime"

if grep -qE '^\s*default\s+localEuler' system/fvSchemes 2>/dev/null; then
    # LTS 伪瞬态：deltaT=1 代表一次迭代，按步数控制输出最干净
    LTS=1
    WI=$(python3 -c "print(max(1, round(($ET-$ST)/$NFRAMES)))")
    foamDictionary -entry writeControl  -set timeStep -disableFunctionEntries system/controlDict >/dev/null
    foamDictionary -entry writeInterval -set "$WI"    -disableFunctionEntries system/controlDict >/dev/null
    ok "检测到 localEuler（LTS 伪瞬态）：endTime=$ET 步，每 $WI 步输出一次 → 约 $NFRAMES 帧"
    warn "LTS 下的“时间”是迭代次数，不是物理秒。视频展示的是尾流波形从静水逐步建立到稳态的过程。"
else
    LTS=0
    WI=$(python3 -c "print(f'{max(($ET-$ST)/$NFRAMES, 1e-9):.6g}')")
    foamDictionary -entry writeControl  -set adjustableRunTime -disableFunctionEntries system/controlDict >/dev/null
    foamDictionary -entry writeInterval -set "$WI"             -disableFunctionEntries system/controlDict >/dev/null
    ok "真实瞬态：endTime=$ET s，每 $WI s 输出一次 → 约 $NFRAMES 帧"
fi
foamDictionary -entry purgeWrite       -set 0      -disableFunctionEntries system/controlDict >/dev/null
foamDictionary -entry writeFormat      -set binary -disableFunctionEntries system/controlDict >/dev/null
foamDictionary -entry writeCompression -set on     -disableFunctionEntries system/controlDict >/dev/null
ok "写二进制 + 压缩（ASCII 格式会占几十 GB）"

# ---------------------------------------------------------------- 估算
step "资源估算"
EST=$(python3 -c "print(f'{$NFRAMES*1.33e6*10*8/1e9*0.45:.1f}')")
AVAIL=$(df -Pk "$HOME" | awk 'NR==2{printf "%.1f", $4/1024/1024}')
echo "    网格量      ~1.33M 单元（半船模型，对称面简化）"
echo "    结果占用    ~${EST} GB（压缩后估算）"
echo "    可用空间    ${AVAIL} GB"
python3 -c "import sys; sys.exit(0 if $AVAIL > $EST*1.5 else 1)" \
    || warn "空间偏紧。把 NFRAMES 调小，或清理磁盘。"

MEMGB=$(awk '/MemTotal/{printf "%.1f", $2/1024/1024}' /proc/meminfo)
echo "    WSL 可用内存 ${MEMGB} GB"
python3 -c "import sys; sys.exit(0 if $MEMGB >= 7.5 else 1)" \
    || warn "内存可能不够（建议 8 GB+）。在 Windows 建 C:\\Users\\L\\.wslconfig 写 [wsl2] 和 memory=16GB，然后 wsl --shutdown。"

# ---------------------------------------------------------------- 完成
cat <<DONE

${G}======================================================
  配置完成
======================================================${N}

算例位置：$RUNDIR

$( if [ -f Allmesh ]; then cat <<'STAGE'
分两阶段跑，中间有检查点：

  阶段 A —— 只建网格（约 10 分钟）：
    cd ~/run/dtcHull
    ./Allmesh
    checkMesh | tail -30          # 看到 "Mesh OK" 再往下

  阶段 B —— 求解（几小时，建议开 tmux 或单独标签）：
    ./Allrun
    （网格已建好，Allrun 会直接复用，接着做 setFields → decomposePar → interFoam → reconstructPar）
STAGE
else cat <<'STAGE'
这个版本的算例没有单独的 Allmesh，建网格包含在 Allrun 里：

    cd ~/run/dtcHull
    ./Allrun

  前 10 分钟左右是建网格（blockMesh → snappyHexMesh），日志在 log.snappyHexMesh；
  之后自动进入求解，日志切到 log.interFoam。想在求解前先检查网格，可以等
  log.snappyHexMesh 结束后另开终端跑：  checkMesh | tail -30
STAGE
fi )

盯进度（另开一个终端）：
    tail -f $RUNDIR/log.interFoam

跑完提取阻力：
    python3 /mnt/c/Users/L/Desktop/openfoam/resistance.py $RUNDIR --full-ship

DONE
