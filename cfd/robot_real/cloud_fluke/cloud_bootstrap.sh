#!/usr/bin/env bash
# 云端机器第一次运行：检查环境 → 后台启动 run_fluke.sh（尾鳍 + 阻力型航速对照，默认 10 个算例 / 13 次计算）
#   Windows PowerShell（两条，各输一次密码）：
#     scp -r C:\Users\L\Desktop\openfoam\robot_real\cloud_fluke <user>@<server>:~/
#     ssh <user>@<server> "bash ~/cloud_fluke/cloud_bootstrap.sh"
#   只跑 L0：    ssh <user>@<server> "ORDER=L0 bash ~/cloud_fluke/cloud_bootstrap.sh"
#   看进度：     ssh <user>@<server> "grep -a '=====\|背景网格\|洞切设置\|第一步\|!!' ~/run/run_fluke.out; bash ~/cloud_fluke/progress.sh ~/run/fluke_L0"
#   拿结果回来： scp -r <user>@<server>:~/cloud_fluke/results_fluke C:\Users\L\Desktop\openfoam\robot_real\leg_dynamics\
#   两路同时跑：性能核 ORDER="..." bash cloud_bootstrap.sh；能效核 ORDER="..." PIN=16-23 OUT=run_fluke_E bash cloud_bootstrap.sh
#   环境变量：OUT（日志名，默认 run_fluke）  ORDER（默认全部）  NP（默认 8）  PIN（绑定的逻辑 CPU，默认 0-15 = 性能核）  WRITES（默认 4）  KEEP_PROC（默认 0）
HERE="$(cd "$(dirname "$0")" && pwd)"
echo "== 机器：$(hostname)  物理核 $(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l)（逻辑 $(nproc)）  内存 $(free -g | awk '/Mem/{print $2}') GB  磁盘可用 $(df -h ~ | awk 'NR==2{print $4}')  $(lsb_release -ds 2>/dev/null || head -1 /etc/os-release)"
if [ -z "${WM_PROJECT_DIR:-}" ]; then
  for rc in /usr/lib/openfoam/openfoam2606/etc/bashrc /usr/lib/openfoam/openfoam2512/etc/bashrc /usr/lib/openfoam/openfoam2506/etc/bashrc /usr/lib/openfoam/openfoam2412/etc/bashrc /usr/lib/openfoam/openfoam2406/etc/bashrc /opt/openfoam*/etc/bashrc "$HOME"/OpenFOAM/OpenFOAM-v*/etc/bashrc; do
    [ -f "$rc" ] && { source "$rc" 2>/dev/null; break; }; done; fi
if [ -n "${WM_PROJECT_DIR:-}" ] && command -v overPimpleDyMFoam >/dev/null; then
  echo "== OpenFOAM：$WM_PROJECT_VERSION  ($WM_PROJECT_DIR)  overPimpleDyMFoam ✓"
else
  echo "!! 没找到 OpenFOAM。装法（Ubuntu，需 sudo）："
  echo "   curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash && sudo apt-get install -y openfoam2606-default"
  echo "   （Ubuntu 20.04 若没有 2606 的包：apt-cache search openfoam 看有哪个 openfoam24xx/25xx-default，装最新的即可）"
  echo "   装完再跑一次：bash ~/cloud_fluke/cloud_bootstrap.sh"; ls /usr/lib/openfoam /opt 2>/dev/null; exit 2
fi
python3 -c "import numpy, matplotlib" 2>/dev/null || { echo "== 装 python 依赖（numpy matplotlib）"; python3 -m pip install --user numpy matplotlib >/dev/null 2>&1; }
python3 -c "import numpy, matplotlib; print('== python3', __import__('sys').version.split()[0], 'numpy', numpy.__version__)" || { echo "!! python3 缺 numpy/matplotlib（sudo apt install python3-numpy python3-matplotlib）"; exit 3; }
command -v mpirun >/dev/null || { echo "!! 没有 mpirun（OpenFOAM 的 openmpi 没装全）"; exit 4; }
mkdir -p ~/run
# 硬盘：这台只剩十几 GB → 每周期只写 4 帧流场（力每步都记，不受影响），跑完删 processor*
export WRITES="${WRITES:-4}" KEEP_PROC="${KEEP_PROC:-0}"
FREE=$(df -BG ~ | awk 'NR==2{gsub("G","",$4); print $4}'); echo "== 硬盘可用 ${FREE} GB（每个算例约 1–2 GB）"; [ "$FREE" -lt 8 ] && { echo "!! 硬盘不到 8 GB，先清理"; exit 5; }
# i9-12900K = 8 性能核（逻辑 CPU 0–15）+ 8 能效核（16–23）。内核 5.4 不识别大小核 → 把求解限定在性能核上，8 进程
export NP="${NP:-8}" OMPI_MCA_hwloc_base_binding_policy=none
# 这台装了 Xilinx XRT：Open MPI 2.1 的 hwloc 枚举 OpenCL 设备时会加载 XRT 的 ICD → 找不到 libxrt_core.so.2 → mpirun 直接崩（连 mpirun hostname 都崩）
mkdir -p "$HOME/.empty_ocl_icd"; export HWLOC_COMPONENTS=-opencl OCL_ICD_VENDORS="$HOME/.empty_ocl_icd"
mpirun -np 2 hostname >/dev/null 2>&1 || { echo "!! mpirun -np 2 hostname 仍然失败："; mpirun -np 2 hostname 2>&1 | tail -5; exit 6; }
PIN="${PIN:-0-15}"
OUT="${OUT:-run_fluke}"      # 第二路（能效核）用 OUT=run_fluke_E，两路各写各的日志
if [ -f ~/run/$OUT.pid ] && kill -0 "$(cat ~/run/$OUT.pid)" 2>/dev/null; then echo "== $OUT 这一路已在跑（PID $(cat ~/run/$OUT.pid)）"; tail -5 ~/run/$OUT.out; exit 0; fi
chmod +x "$HERE"/*.sh 2>/dev/null
nohup taskset -c "$PIN" bash "$HERE/run_fluke.sh" > ~/run/$OUT.out 2>&1 &
echo $! > ~/run/$OUT.pid
sleep 20; echo "== 已启动 $OUT（PID $!，CPU $PIN，$NP 进程）ORDER=${ORDER:-默认全部}；日志 ~/run/$OUT.out 前 20 行："; head -20 ~/run/$OUT.out
echo "== 退出 ssh 不影响计算。看进度见本文件开头的命令。"
