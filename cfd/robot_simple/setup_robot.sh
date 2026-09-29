#!/usr/bin/env bash
# =====================================================================
#  robot_simple —— 自定义几何算例配置脚本（以 ESI v2606 的 DTCHull 为模板）
#
#  用法（Ubuntu 终端）：
#      bash /mnt/c/Users/L/Desktop/openfoam/robot_simple/setup_robot.sh
#
#  可调参数（环境变量）：
#      U=0.2        巡航速度 m/s（来流沿 -x，前进方向 +x）
#      T=0.012      水线高度 m（几何的 z 坐标；身体底面 z=0）
#      ENDTIME=2000 LTS 伪时间步数
#      NFRAMES=40   视频帧数
#      COARSE=1     1 = 粗网格（快，几十万单元）；0 = DTC 同等密度
#      LAYERS=1     1 = 在壁面加 3 层边界层（正式算阻力用）；0 = 不加（验证流程时快很多）
#      PARMESH=0    1 = snappyHexMesh 并行建网格（快 4-6 倍）；0 = 和 DTC 模板一样串行
#      FULL=0       1 = 全模型（几何左右不对称时必须用；计算域 y 两侧对称，阻力不用 ×2）
#      NLEVELS=6    refineMesh 逐级加密的级数（真几何建议 5，配合 SNAPLEVEL=1）
#      SNAPLEVEL=0  snappyHexMesh 在物体表面再加密的级数（真几何建议 1：薄桨板才有足够单元）
#      GEOM=...     几何 STL 路径（二进制 STL，米制，封闭）
#
#  思路：几何相关的字典（计算域 / 加密盒 / 水面 / hRef / locationInMesh）全部重写，
#        数值格式、边界条件、湍流设置、Allrun 流程原样沿用 DTCHull（已验证）。
# =====================================================================
set -o pipefail

U=${U:-0.2}
T=${T:-0.012}
ENDTIME=${ENDTIME:-2000}
NFRAMES=${NFRAMES:-40}
COARSE=${COARSE:-1}
LAYERS=${LAYERS:-1}
PARMESH=${PARMESH:-0}
FULL=${FULL:-0}
NLEVELS=${NLEVELS:-6}
SNAPLEVEL=${SNAPLEVEL:-0}
GEOM=${GEOM:-/mnt/c/Users/L/Desktop/openfoam/robot_simple/hull.stl}
RUNDIR="${RUNDIR:-$HOME/run/robot_simple}"

B=$'\033[1;36m'; G=$'\033[0;32m'; Y=$'\033[1;33m'; R=$'\033[0;31m'; N=$'\033[0m'
step(){ printf "\n${B}==> %s${N}\n" "$*"; }
ok(){   printf "${G}    OK  %s${N}\n" "$*"; }
warn(){ printf "${Y}    !!  %s${N}\n" "$*"; }
die(){  printf "\n${R}    XX  %s${N}\n\n" "$*"; exit 1; }

# ------------------------------------------------------------ 0 环境
step "检查环境"
[ -n "$WM_PROJECT_DIR" ] || die "OpenFOAM 环境没加载。先 source ~/.bashrc"
[ -f "$GEOM" ] || die "找不到几何文件 $GEOM"
ok "$(foamVersion 2>/dev/null || echo "$WM_PROJECT_DIR")"
ok "几何：$GEOM"

# ------------------------------------------------------------ 1 模板
step "复制 DTCHull 作为模板"
SRC=$(find "$FOAM_TUTORIALS" -maxdepth 5 -type d -name 'DTCHull' 2>/dev/null | head -1)
[ -n "$SRC" ] || die "找不到 DTCHull 教程"
if [ -d "$RUNDIR" ]; then
    warn "$RUNDIR 已存在"
    read -r -p "    删掉重来？[y/N] " ans
    case "$ans" in y|Y) rm -rf "$RUNDIR";; *) die "已取消";; esac
fi
mkdir -p "$(dirname "$RUNDIR")"
cp -r "$SRC" "$RUNDIR" || die "复制失败"
cd "$RUNDIR" || die "进不去 $RUNDIR"
ok "模板：${SRC/$FOAM_TUTORIALS/\$FOAM_TUTORIALS}"

# 结构核对：blockMesh 的边界名、0.orig 里显式写的边界、分级加密字典
for p in atmosphere inlet outlet bottom side midPlane; do
    grep -q "^\s*$p\s*$" system/blockMeshDict || die "模板 blockMeshDict 里没有边界 '$p'，结构不同，把 system/blockMeshDict 发我"
done
for p in inlet outlet atmosphere hull; do
    grep -q "^\s*$p\s*$" 0.orig/U || die "模板 0.orig/U 里没有边界 '$p'，把 0.orig/U 发我"
done
[ -f system/topoSetDict.1 ] && [ -f system/topoSetDict.6 ] || die "模板没有 system/topoSetDict.1…6，结构不同，把 system/ 目录列表发我"
[ -f system/refineMeshDict ] || die "模板没有 system/refineMeshDict"
grep -q "DTC-scaled" system/snappyHexMeshDict || die "snappyHexMeshDict 里找不到 DTC-scaled，结构不同"
ok "模板结构核对通过（ESI 版 DTCHull）"

# ------------------------------------------------------------ 2 几何
step "放入几何"
mkdir -p constant/triSurface
rm -f constant/triSurface/DTC-scaled.*
cp "$GEOM" constant/triSurface/hull.stl
ok "constant/triSurface/hull.stl"

# 字典里所有对 DTC-scaled 的引用改成 hull（snappyHexMeshDict / surfaceFeatureExtractDict）
grep -rl "DTC-scaled" system 2>/dev/null | xargs -r sed -i 's/DTC-scaled/hull/g'
grep -q 'hull\.stl' system/snappyHexMeshDict && grep -q 'hull\.eMesh' system/snappyHexMeshDict \
    || die "snappyHexMeshDict 里的几何引用没改成功"
grep -q '^hull\.stl' system/surfaceFeatureExtractDict || die "surfaceFeatureExtractDict 的键没改成功"
ok "snappyHexMeshDict / surfaceFeatureExtractDict 已指向 hull.stl / hull.eMesh"

# Allrun：照抄 ESI 的流程，只去掉"从教程资源目录复制 DTC 几何"那一段
rm -f Allrun
if [ "$PARMESH" = "1" ]; then
cat > Allrun <<'EOF'
#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
# robot_simple（并行建网格版）：几何已由 setup_robot.sh 放在 constant/triSurface/hull.stl
# 流程参考 ESI 的 motorBike 教程：背景网格串行 → 分块 → snappyHexMesh 并行 → 求解并行 → 重组

runApplication surfaceFeatureExtract

runApplication blockMesh

for i in $(seq 1 NLEVELS_PLACEHOLDER)
do
    runApplication -s "$i" \
        topoSet -dict system/topoSetDict.${i}

    runApplication -s "$i" \
        refineMesh -dict system/refineMeshDict -overwrite
done

runApplication decomposePar

runParallel snappyHexMesh -overwrite

runParallel checkMesh -constant

restore0Dir -processor

runParallel setFields

runParallel renumberMesh -overwrite

runParallel $(getApplication)

# 合并网格和结果（ParaView 选 "Decomposed Case" 的话这两步可以省掉）
runApplication reconstructParMesh -constant

runApplication reconstructPar

#------------------------------------------------------------------------------
EOF
ok "Allrun 已写入（PARMESH=1：snappyHexMesh 并行）"
else
cat > Allrun <<'EOF'
#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
# robot_simple：几何已由 setup_robot.sh 放在 constant/triSurface/hull.stl

runApplication surfaceFeatureExtract

runApplication blockMesh

for i in $(seq 1 NLEVELS_PLACEHOLDER)
do
    runApplication -s "$i" \
        topoSet -dict system/topoSetDict.${i}

    runApplication -s "$i" \
        refineMesh -dict system/refineMeshDict -overwrite
done

runApplication snappyHexMesh -overwrite

restore0Dir

runApplication setFields

runParallel -s decompose redistributePar -decompose

runParallel renumberMesh -overwrite

runParallel $(getApplication)

runParallel -s reconstruct redistributePar -reconstruct

#------------------------------------------------------------------------------
EOF
ok "Allrun 已重写（流程与 ESI 模板一致，去掉了复制 DTC 几何的步骤）"
fi
sed -i "s/NLEVELS_PLACEHOLDER/$NLEVELS/" Allrun
chmod +x Allrun
ok "refineMesh 逐级加密 $NLEVELS 级"

# snappy 表面加密
if [ "$SNAPLEVEL" != "0" ]; then
    sed -i -E "/refinementSurfaces/,/^\s*}/ s/level \(0 0\);/level ($SNAPLEVEL $SNAPLEVEL);/" system/snappyHexMeshDict
    sed -i -E "/features/,/\);/ s/level\s+0;/level $SNAPLEVEL;/" system/snappyHexMeshDict
    grep -q "level ($SNAPLEVEL $SNAPLEVEL)" system/snappyHexMeshDict || die "snappyHexMeshDict 的表面加密级数没改成功"
    ok "SNAPLEVEL=$SNAPLEVEL：snappyHexMesh 在表面和特征边再加密 $SNAPLEVEL 级"
fi

# 边界层开关
if [ "$LAYERS" = "0" ]; then
    foamDictionary -entry addLayers -set false -disableFunctionEntries system/snappyHexMeshDict >/dev/null \
        || die "改 snappyHexMeshDict 的 addLayers 失败"
    ok "LAYERS=0：不加边界层（验证流程用；正式算阻力请用 LAYERS=1）"
else
    ok "LAYERS=1：壁面加 3 层边界层（沿用 DTC 设置）"
fi

# ------------------------------------------------------------ 3 尺寸
step "计算域尺寸（前进 +x，来流 -x，$([ "$FULL" = "1" ] && echo 全模型 || echo 半模型 y ≤ 0)）"
# 纯标准库读 STL 包围盒（二进制或 ASCII 都行，不依赖 numpy）
read -r XMN YMN ZMN XMX YMX ZMX < <(python3 - "$GEOM" <<'PY'
import sys, struct, os, re
p = sys.argv[1]
size = os.path.getsize(p)
with open(p, "rb") as f:
    head = f.read(84)
n = struct.unpack("<I", head[80:84])[0] if size >= 84 else -1
mn = [float("inf")] * 3; mx = [float("-inf")] * 3
if size == 84 + 50 * n:                      # 二进制
    with open(p, "rb") as f:
        f.seek(84); data = f.read(50 * n)
    for rec in struct.iter_unpack("<12fH", data):
        for k in range(3):
            for j in range(3):
                v = rec[3 + 3 * k + j]
                if v < mn[j]: mn[j] = v
                if v > mx[j]: mx[j] = v
else:                                        # ASCII
    txt = open(p, "r", errors="replace").read()
    for m in re.finditer(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt):
        for j in range(3):
            v = float(m.group(j + 1))
            if v < mn[j]: mn[j] = v
            if v > mx[j]: mx[j] = v
if mn[0] == float("inf"):
    sys.exit(1)
print(*[f"{x:.5f}" for x in (*mn, *mx)])
PY
)
[ -n "$ZMX" ] || die "读不出 STL 包围盒（文件不是合法的 STL？）"
LX=$(python3 -c "print(f'{$XMX-$XMN:.4f}')")
echo "    几何包围盒  x[$XMN, $XMX]  y[$YMN, $YMX]  z[$ZMN, $ZMX]   L = $LX m"
python3 -c "import sys; sys.exit(0 if $ZMN < $T < $ZMX else 1)" \
    || warn "水线 T=$T 不在几何的 z 范围 [$ZMN, $ZMX] 内 —— 物体要么全在水下要么全在水上，确认一下 T"

# 计算域：比例取自 DTC（上游 1.6L，下游 4.2L，侧 3.1L，下 2.6L，上 0.65L）
X0=$(python3 -c "print(f'{$XMN-4.2*$LX:.4f}')");  X1=$(python3 -c "print(f'{$XMX+1.6*$LX:.4f}')")
Y0=$(python3 -c "print(f'{-3.1*$LX:.4f}')")
if [ "$FULL" = "1" ]; then Y1=$(python3 -c "print(f'{3.1*$LX:.4f}')"); else Y1=0; fi
Z0=$(python3 -c "print(f'{$T-2.6*$LX:.4f}')");     Z6=$(python3 -c "print(f'{$T+0.65*$LX:.4f}')")
# 竖向分层：底 → 几何最低点下方 3 cm → 水线 ±1.8 mm 细层 → 几何顶上 2 cm → 顶
Z1=$(python3 -c "print(f'{$ZMN-0.03:.4f}')")
Z2=$(python3 -c "print(f'{$T-0.0018:.4f}')");  Z3=$T;  Z4=$(python3 -c "print(f'{$T+0.0018:.4f}')")
Z5=$(python3 -c "print(f'{$ZMX+0.02:.4f}')")
if [ "$COARSE" = "1" ]; then NX=21; NY=10; NZA=25; NZB=25; NZC=20; NZD=10; else NX=42; NY=19; NZA=50; NZB=50; NZC=40; NZD=20; fi
echo "    计算域      x[$X0, $X1]  y[$Y0, $Y1]  z[$Z0, $Z6]"
echo "    竖向分层    $Z0 | $Z1 | $Z2 | $Z3(水线) | $Z4 | $Z5 | $Z6"
echo "    背景网格    ${NX} × ${NY} × (${NZA}+${NZB}+4+4+${NZC}+${NZD})   COARSE=$COARSE"

cat > system/blockMeshDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// robot_simple：分块思路照搬 DTCHull，尺寸由 setup_robot.sh 生成
scale   1;

vertices
(
    ($X0 $Y0 $Z0) ($X1 $Y0 $Z0) ($X1 $Y1 $Z0) ($X0 $Y1 $Z0)
    ($X0 $Y0 $Z1) ($X1 $Y0 $Z1) ($X1 $Y1 $Z1) ($X0 $Y1 $Z1)
    ($X0 $Y0 $Z2) ($X1 $Y0 $Z2) ($X1 $Y1 $Z2) ($X0 $Y1 $Z2)
    ($X0 $Y0 $Z3) ($X1 $Y0 $Z3) ($X1 $Y1 $Z3) ($X0 $Y1 $Z3)
    ($X0 $Y0 $Z4) ($X1 $Y0 $Z4) ($X1 $Y1 $Z4) ($X0 $Y1 $Z4)
    ($X0 $Y0 $Z5) ($X1 $Y0 $Z5) ($X1 $Y1 $Z5) ($X0 $Y1 $Z5)
    ($X0 $Y0 $Z6) ($X1 $Y0 $Z6) ($X1 $Y1 $Z6) ($X0 $Y1 $Z6)
);

blocks
(
    hex (0 1 2 3 4 5 6 7)         ($NX $NY $NZA) simpleGrading (1 1 0.05)
    hex (4 5 6 7 8 9 10 11)       ($NX $NY $NZB) simpleGrading (1 1 1)
    hex (8 9 10 11 12 13 14 15)   ($NX $NY 4)    simpleGrading (1 1 1)
    hex (12 13 14 15 16 17 18 19) ($NX $NY 4)    simpleGrading (1 1 1)
    hex (16 17 18 19 20 21 22 23) ($NX $NY $NZC) simpleGrading (1 1 1)
    hex (20 21 22 23 24 25 26 27) ($NX $NY $NZD) simpleGrading (1 1 5)
);

edges
(
);

boundary
(
    atmosphere
    {
        type patch;
        faces ( (24 25 26 27) );
    }
    inlet
    {
        type patch;
        faces ( (1 2 6 5) (5 6 10 9) (9 10 14 13) (13 14 18 17) (17 18 22 21) (21 22 26 25) );
    }
    outlet
    {
        type patch;
        faces ( (0 4 7 3) (4 8 11 7) (8 12 15 11) (12 16 19 15) (16 20 23 19) (20 24 27 23) );
    }
    bottom
    {
        type symmetryPlane;
        faces ( (0 3 2 1) );
    }
    side
    {
        type symmetryPlane;
        faces ( (0 1 5 4) (4 5 9 8) (8 9 13 12) (12 13 17 16) (16 17 21 20) (20 21 25 24) );
    }
    midPlane
    {
        type symmetryPlane;
        faces ( (3 7 6 2) (7 11 10 6) (11 15 14 10) (15 19 18 14) (19 23 22 18) (23 27 26 22) );
    }
);

mergePatchPairs
(
);
EOF
ok "system/blockMeshDict 已写入"

# 六级嵌套加密盒：ESI 版每级一个文件 topoSetDict.i，集合名固定 c0（refineMeshDict 里 set c0）
python3 - "$XMN" "$XMX" "$ZMN" "$ZMX" "$Y0" "$YMX" "$NLEVELS" "$FULL" <<'PY'
import sys
xmn,xmx,zmn,zmx,y0,ymx=map(float,sys.argv[1:7]); nlev=int(sys.argv[7]); full=sys.argv[8]=="1"
# (下游, 上游, 侧向, 下方, 上方) 边距，逐级收紧
M=[(0.32,0.12,0.20,0.10,0.10),(0.16,0.06,0.10,0.05,0.05),(0.09,0.04,0.06,0.03,0.03),
   (0.05,0.025,0.035,0.02,0.02),(0.03,0.015,0.02,0.012,0.012),(0.016,0.008,0.012,0.007,0.007)]
for i,(d,u,s,b,a) in enumerate(M[:nlev],1):
    ylo=max(y0,-(ymx+s)); yhi=(ymx+s) if full else 0
    txt=f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}}
// robot_simple 第 {i} 级加密盒，由 setup_robot.sh 生成
actions
(
    {{
        name    c0;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     ({xmn-d:.4f} {ylo:.4f} {zmn-b:.4f}) ({xmx+u:.4f} {yhi:.4f} {zmx+a:.4f});
    }}
);
"""
    open(f"system/topoSetDict.{i}","w").write(txt)
    print(f"    c{i}: box ({xmn-d:.4f} {ylo:.4f} {zmn-b:.4f}) ({xmx+u:.4f} {yhi:.4f} {zmx+a:.4f})")
PY
ok "system/topoSetDict.1 … .$NLEVELS 已写入"

cat > system/setFieldsDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      setFieldsDict;
}
defaultFieldValues
(
    volScalarFieldValue alpha.water 0
);
regions
(
    boxToCell
    {
        box (-999 -999 -999) (999 999 $T);
        fieldValues ( volScalarFieldValue alpha.water 1 );
    }
    boxToFace
    {
        box (-999 -999 -999) (999 999 $T);
        fieldValues ( volScalarFieldValue alpha.water 1 );
    }
);
EOF
ok "system/setFieldsDict 已写入（水线 z = $T）"

foamDictionary -entry value -set "$T" -disableFunctionEntries constant/hRef >/dev/null \
    || die "改 constant/hRef 失败"
ok "constant/hRef = $T（p_rgh 的参考高度，必须和水线一致）"

LIM=$(python3 -c "print(f'({$XMX+0.77*$LX:.4f} {-0.107*$LX:.4f} {$T-0.118*$LX:.4f})')")
sed -i -E "s/locationInMesh[[:space:]]*\([^)]*\);/locationInMesh $LIM;/" system/snappyHexMeshDict
grep -q "locationInMesh $LIM" system/snappyHexMeshDict || die "snappyHexMeshDict 里没改到 locationInMesh"
ok "snappyHexMeshDict locationInMesh = $LIM（物体前方水下）"

# ------------------------------------------------------------ 4 来流 / 湍流初值 / 输出
step "来流速度、湍流初值、力的参考点、输出控制"
FILES=$(grep -rl -- "1\.668" 0.orig system constant 2>/dev/null)
[ -n "$FILES" ] || die "没找到 DTC 的来流速度 1.668，模板结构不同"
echo "$FILES" | xargs sed -i "s/-1\.668/-$U/g; s/1\.668/$U/g"
ok "来流 |U| = $U m/s（改动：$(echo $FILES | tr '\n' ' ')）"

# k、omega 初值按新速度和尺度粗略缩放（湍流强度 2%，长度尺度 0.07 L）
KIN=$(python3 -c "print(f'{1.5*(0.02*$U)**2:.3e}')")
OMG=$(python3 -c "import math; k=1.5*(0.02*$U)**2; print(f'{math.sqrt(k)/(0.09**0.25*0.07*$LX):.3g}')")
foamDictionary -entry internalField -set "uniform $KIN" -disableFunctionEntries 0.orig/k     >/dev/null
foamDictionary -entry internalField -set "uniform $OMG" -disableFunctionEntries 0.orig/omega >/dev/null
ok "k = $KIN m²/s²，omega = $OMG 1/s"

ZC=$(python3 -c "print(f'{($ZMN+$ZMX)/2:.4f}')")
sed -i -E "s/(CofR[[:space:]]+)\([^)]*\)/\1(0 0 $ZC)/" system/controlDict
ok "forces CofR = (0 0 $ZC)（只影响力矩）"

foamDictionary -entry endTime -set "$ENDTIME" -disableFunctionEntries system/controlDict >/dev/null
if grep -qE '^\s*default\s+localEuler' system/fvSchemes; then
    WI=$(python3 -c "print(max(1, round($ENDTIME/$NFRAMES)))")
    foamDictionary -entry writeControl  -set timeStep -disableFunctionEntries system/controlDict >/dev/null
    foamDictionary -entry writeInterval -set "$WI"    -disableFunctionEntries system/controlDict >/dev/null
    ok "LTS：endTime=$ENDTIME 步，每 $WI 步写一帧 → ~$NFRAMES 帧"
else
    warn "fvSchemes 不是 localEuler，输出控制沿用模板"
fi
foamDictionary -entry purgeWrite       -set 0      -disableFunctionEntries system/controlDict >/dev/null
foamDictionary -entry writeFormat      -set binary -disableFunctionEntries system/controlDict >/dev/null
foamDictionary -entry writeCompression -set on     -disableFunctionEntries system/controlDict >/dev/null

PHYS=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l); [ "$PHYS" -ge 1 ] 2>/dev/null || PHYS=$(nproc)
NSUB=$PHYS; [ "$NSUB" -gt 16 ] && NSUB=16
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
ok "并行分解 $NSUB 块（scotch）"

# ------------------------------------------------------------ 5 立即自检
step "立即自检（几何 + 背景网格，几秒钟）"
mkdir -p tmpCheck && ( cd tmpCheck && surfaceCheck ../constant/triSurface/hull.stl 2>&1 \
    | grep -E "Bounding Box|closed|Surface has|unconnected|orientation|Triangles" | sed 's/^/    /' )
rm -rf tmpCheck
blockMesh > log.blockMesh.setup 2>&1 || { tail -20 log.blockMesh.setup; die "blockMesh 失败，把上面的输出发我"; }
NCELL=$(grep -m1 "nCells" log.blockMesh.setup | awk '{print $NF}')
ok "blockMesh 通过：背景网格 $NCELL 个单元"
checkMesh 2>&1 | grep -E "Mesh OK|Failed|cells:" | sed 's/^/    /'
rm -f log.blockMesh.setup

# ------------------------------------------------------------ 完成
cat <<DONE

${G}======================================================
  robot_simple 配置完成
======================================================${N}
算例：$RUNDIR
设定：U = $U m/s，水线 z = $T m，endTime = $ENDTIME（LTS），$NSUB 核
      COARSE=$COARSE LAYERS=$LAYERS PARMESH=$PARMESH FULL=$FULL NLEVELS=$NLEVELS SNAPLEVEL=$SNAPLEVEL

跑：
    cd $RUNDIR
    ./Allrun
    tail -f log.interFoam        # 另开终端

建网格阶段（surfaceFeatureExtract → blockMesh → topoSet/refineMesh ×6 → snappyHexMesh）
串行加边界层约 10–20 分钟；PARMESH=1 + LAYERS=0 约 2–4 分钟。看到 log.interFoam 出现就是在求解了。

下次想快：
    PARMESH=1 LAYERS=0 bash /mnt/c/Users/L/Desktop/openfoam/robot_simple/setup_robot.sh

跑完：
    python3 /mnt/c/Users/L/Desktop/openfoam/resistance.py $RUNDIR$([ "$FULL" = "1" ] && echo "        # 全模型：不要加 --full-ship" || echo " --full-ship")

如果哪一步报错，把对应的 log.xxx 最后 30 行发我。
DONE
