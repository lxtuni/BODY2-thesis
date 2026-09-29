#!/usr/bin/env bash
# =====================================================================
#  OpenFOAM-v2012 (windows10 / WSL 预编译包) 一键安装脚本
#
#  在 WSL 的 Ubuntu 终端里运行：
#      bash /mnt/c/Users/L/Desktop/openfoam/install_openfoam.sh
#
#  开头会问一次 sudo 密码（你的 Linux 密码，不是 Windows 密码），
#  之后全程无人值守。整个过程约 10-20 分钟。
# =====================================================================

set -o pipefail   # 故意不用 set -e：OpenFOAM 的 etc/bashrc 在严格模式下会误退出

TGZ="/mnt/c/Users/L/Desktop/openfoam/OpenFOAM-v2012-windows10.tgz"
FOAM_DIR="/opt/OpenFOAM/OpenFOAM-v2012"
NEED_KB=$((7 * 1024 * 1024))   # 解压需要约 7 GB

R=$'\033[0;31m'; G=$'\033[0;32m'; Y=$'\033[1;33m'; B=$'\033[1;36m'; N=$'\033[0m'
step(){ printf "\n${B}==> [%s] %s${N}\n" "$1" "$2"; }
ok(){   printf "${G}    OK  %s${N}\n" "$*"; }
warn(){ printf "${Y}    !!  %s${N}\n" "$*"; }
die(){  printf "\n${R}    XX  %s${N}\n\n" "$*"; exit 1; }

printf "${B}"
echo "======================================================"
echo "  OpenFOAM-v2012  WSL 安装脚本"
echo "======================================================"
printf "${N}"

# ---------------------------------------------------------------- [0]
step 0/6 "环境自检"

[ -n "$BASH_VERSION" ] || die "请用 bash 运行：bash $0"

if grep -qi microsoft /proc/version 2>/dev/null; then
    ok "运行在 WSL 里（distro: ${WSL_DISTRO_NAME:-未知}）"
else
    warn "没检测到 WSL —— 如果这是原生 Linux 也没问题，继续。"
fi

if [ -r /etc/os-release ]; then
    . /etc/os-release
    ok "系统：$PRETTY_NAME"
    case "$VERSION_ID" in
        20.04|18.04) ;;
        *) warn "这个包是 GCC 9.2 / Ubuntu 20.04 时代编译的。"
           warn "在 $VERSION_ID 上可能缺库 —— 第 5 步会自动检测并告诉你怎么办。" ;;
    esac
fi

[ -f "$TGZ" ] || die "找不到安装包：$TGZ"
SZ=$(stat -c%s "$TGZ")
[ "$SZ" -gt 600000000 ] || die "安装包大小异常（$SZ 字节），可能没下完，请重新下载。"
ok "安装包就位（$((SZ/1024/1024)) MB）"

AVAIL=$(df -Pk /opt | awk 'NR==2{print $4}')
if [ "$AVAIL" -lt "$NEED_KB" ]; then
    die "磁盘空间不足：/opt 可用 $((AVAIL/1024/1024)) GB，需要约 7 GB。"
fi
ok "磁盘空间够（可用 $((AVAIL/1024/1024)) GB）"

echo
echo "接下来需要管理员权限。请输入你的 ${B}Linux 密码${N}（输入时不显示字符）："
sudo -v || die "sudo 认证失败。"
# 后台续期，避免中途再次询问密码
( while true; do sleep 50; sudo -n true 2>/dev/null; kill -0 "$$" 2>/dev/null || exit; done ) &
KEEPALIVE=$!
trap 'kill $KEEPALIVE 2>/dev/null' EXIT
ok "已获取管理员权限"

# ---------------------------------------------------------------- [1]
step 1/6 "安装编译依赖 (bison / flex / m4 / build-essential)"
sudo apt-get update -qq          || warn "apt update 有警告，继续。"
sudo apt-get install -y -qq bison flex m4 build-essential \
    || die "依赖安装失败，检查一下 WSL 的网络（试试 ping archive.ubuntu.com）。"
ok "依赖安装完成"

# ---------------------------------------------------------------- [2]
step 2/6 "解压到 /opt（最慢的一步，5-15 分钟）"
if [ -d "$FOAM_DIR" ]; then
    ok "$FOAM_DIR 已存在，跳过解压"
else
    echo "    解压中，每个点代表 2000 个文件："
    printf "    "
    sudo tar -xzvf "$TGZ" -C /opt 2>/dev/null \
        | awk 'NR%2000==0{printf "."; fflush()} END{print ""}'
    [ "${PIPESTATUS[0]}" -eq 0 ] || die "解压失败，安装包可能损坏，请重新下载。"
    ok "解压完成"
fi

# ---------------------------------------------------------------- [3]
step 3/6 "移交目录所有权"
sudo chown -R "$USER" /opt/OpenFOAM || die "chown 失败。"
ok "/opt/OpenFOAM 现在归 $USER 所有"

# ---------------------------------------------------------------- [4]
step 4/6 "配置 ~/.bashrc 环境变量"
FOAM_BASHRC=$(find /opt/OpenFOAM -maxdepth 3 -path '*/etc/bashrc' -type f 2>/dev/null | head -1)
[ -n "$FOAM_BASHRC" ] || die "在 /opt/OpenFOAM 下找不到 etc/bashrc，解压结构不对。"
ok "找到环境脚本：$FOAM_BASHRC"

if grep -qF "$FOAM_BASHRC" "$HOME/.bashrc" 2>/dev/null; then
    ok "~/.bashrc 里已有这一行，不重复添加"
else
    cp "$HOME/.bashrc" "$HOME/.bashrc.bak.$(date +%Y%m%d%H%M%S)" 2>/dev/null
    {
        echo ""
        echo "# --- OpenFOAM-v2012 (由安装脚本添加) ---"
        echo "source $FOAM_BASHRC"
    } >> "$HOME/.bashrc"
    ok "已写入 ~/.bashrc（原文件已备份为 ~/.bashrc.bak.*）"
fi

# ---------------------------------------------------------------- [5]
step 5/6 "验证可执行文件与动态库"
BM=$(bash -c "source '$FOAM_BASHRC' >/dev/null 2>&1; command -v blockMesh" 2>/dev/null)
[ -n "$BM" ] || die "环境加载后仍找不到 blockMesh，安装不完整。"
ok "blockMesh 路径：$BM"

MISSING=$(bash -c "source '$FOAM_BASHRC' >/dev/null 2>&1; ldd '$BM'" 2>/dev/null \
          | grep 'not found' | awk '{print $1}' | sort -u)
if [ -n "$MISSING" ]; then
    printf "${R}\n    检测到缺失的动态库：${N}\n"
    echo "$MISSING" | sed 's/^/      - /'
    cat <<'HINT'

    这就是 v2012 预编译包在新版 Ubuntu 上的典型问题。两个解决办法：

      A) 装一个 Ubuntu 20.04 的 WSL 发行版专门跑 OpenFOAM（推荐）
         在 Windows PowerShell 里：  wsl --install -d Ubuntu-20.04
         然后在新的 Ubuntu 20.04 里重新运行本脚本。

      B) 改用 openfoam.com 的 apt 源装新版本（v2412），不需要这个 tgz。

    把上面这几行库名发给我，我帮你判断哪种更省事。
HINT
    exit 1
fi
ok "动态库检查通过，没有缺失"

# ---------------------------------------------------------------- [6]
step 6/6 "跑通 cavity 教程算例（icoFoam）"
TESTROOT="$HOME/OpenFOAM-test"
bash -c "
    source '$FOAM_BASHRC' >/dev/null 2>&1  || exit 10
    rm -rf '$TESTROOT' && mkdir -p '$TESTROOT' || exit 11
    cp -r \"\$FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity\" '$TESTROOT/' || exit 12
    cd '$TESTROOT/cavity' || exit 13
    blockMesh > log.blockMesh 2>&1 || exit 14
    icoFoam   > log.icoFoam   2>&1 || exit 15
    touch cavity.foam
"
RC=$?
case $RC in
    0)  ok "cavity 算例跑通了" ;;
    12) die "找不到教程文件，FOAM_TUTORIALS 变量没设好。" ;;
    14) die "blockMesh 失败，看日志：$TESTROOT/cavity/log.blockMesh" ;;
    15) die "icoFoam 失败，看日志：$TESTROOT/cavity/log.icoFoam" ;;
    *)  die "验证失败（错误码 $RC）。" ;;
esac

grep -q '^End' "$TESTROOT/cavity/log.icoFoam" && ok "求解器正常收敛结束（End）"

# ---------------------------------------------------------------- 完成
DISTRO="${WSL_DISTRO_NAME:-Ubuntu}"
# 把 Linux 路径转成 Windows 资源管理器能用的 UNC 路径
WSLPATH='\\wsl$\'"$DISTRO"$(printf '%s' "$TESTROOT/cavity" | tr '/' '\\')
printf "${G}"
cat <<DONE

======================================================
  安装完成
======================================================
DONE
printf "${N}"
cat <<DONE
下一步：

1. 关掉这个终端，重新打开一个 Ubuntu 终端（让环境变量生效）。
   或者在当前终端执行：  source ~/.bashrc

2. 确认环境：      echo \$WM_PROJECT_DIR
   应该输出：      $FOAM_DIR

3. 测试算例在：    $TESTROOT/cavity
   在 Windows 资源管理器地址栏输入下面这行，可以直接用 ParaView 打开
   里面的 cavity.foam：

       $WSLPATH

   （ParaView 从 https://www.paraview.org/download/ 下 Windows 版）

4. 以后所有算例都放在 Linux 家目录（~/ 下面），不要放 /mnt/c/，
   否则 I/O 会慢一个数量级。

DONE
