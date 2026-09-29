#!/usr/bin/env bash
# 云端机器上第一次运行：检查环境 → 启动 run_cloud.sh 批处理（最小决定集 V3c12 V3c06 V3trap）
#   Windows PowerShell 里（两条，各输一次密码）：
#     scp -r C:\Users\L\Desktop\openfoam\robot_real\leg_dynamics <user>@<server>:~/
#     ssh <user>@<server> "bash ~/leg_dynamics/cloud_bootstrap.sh"
#   之后看进度： ssh <user>@<server> "tail -20 ~/run/run_cloud.out; bash ~/leg_dynamics/progress.sh ~/run/leg_V3c12_closed"
#   环境变量：ORDER（默认 "V3c12 V3c06 V3trap"）NP（默认 min(nproc,16)）
HERE="$(cd "$(dirname "$0")" && pwd)"; export ORDER="${ORDER:-V3c12 V3c06 V3trap}"
echo "== 机器：$(hostname)  $(nproc) 核  内存 $(free -g | awk '/Mem/{print $2}') GB  磁盘可用 $(df -h ~ | awk 'NR==2{print $4}')  $(lsb_release -ds 2>/dev/null || cat /etc/os-release | head -1)"
if [ -z "${WM_PROJECT_DIR:-}" ]; then
  for rc in /usr/lib/openfoam/openfoam2606/etc/bashrc /opt/openfoam2606/etc/bashrc "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc" /usr/lib/openfoam/openfoam2512/etc/bashrc /usr/lib/openfoam/openfoam2506/etc/bashrc /usr/lib/openfoam/openfoam2412/etc/bashrc; do
    [ -f "$rc" ] && { source "$rc" 2>/dev/null; break; }; done; fi
if [ -n "${WM_PROJECT_DIR:-}" ] && command -v overPimpleDyMFoam >/dev/null; then
  echo "== OpenFOAM：$WM_PROJECT_VERSION  ($WM_PROJECT_DIR)  overPimpleDyMFoam ✓"
else
  echo "!! 没找到 OpenFOAM（v2606/2512/2506/2412 的 /usr/lib/openfoam 或 /opt 或 ~/OpenFOAM）。装法（Ubuntu，需 sudo）："
  echo "   curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash && sudo apt-get install -y openfoam2606-default"
  echo "   装完再跑一次：bash ~/leg_dynamics/cloud_bootstrap.sh"; ls /usr/lib/openfoam /opt 2>/dev/null; exit 2
fi
python3 -c "import numpy, matplotlib" 2>/dev/null || { echo "== 装 python 依赖（numpy matplotlib）"; pip3 install --user numpy matplotlib >/dev/null 2>&1 || python3 -m pip install --user numpy matplotlib >/dev/null 2>&1; }
python3 -c "import numpy, matplotlib; print('== python3', __import__('sys').version.split()[0], 'numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)" || { echo "!! python3 缺 numpy/matplotlib，装不上（没网？）"; exit 3; }
command -v mpirun >/dev/null || { echo "!! 没有 mpirun（OpenFOAM 的 openmpi 没装全）"; exit 4; }
mkdir -p ~/run
if pgrep -f run_cloud.sh >/dev/null; then echo "== run_cloud.sh 已在跑：$(pgrep -f run_cloud.sh | head -1)"; tail -5 ~/run/run_cloud.out; exit 0; fi
chmod +x "$HERE"/*.sh 2>/dev/null
nohup bash "$HERE/run_cloud.sh" > ~/run/run_cloud.out 2>&1 &
sleep 20; echo "== 已启动（PID $!）ORDER=$ORDER；前 20 行："; head -20 ~/run/run_cloud.out
echo "== 看进度：ssh 上来后  tail -f ~/run/run_cloud.out   或  bash ~/leg_dynamics/progress.sh ~/run/leg_V3c12_closed"
