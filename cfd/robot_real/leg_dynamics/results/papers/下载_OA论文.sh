#!/usr/bin/env bash
# 把开放获取（OA）的必读论文下载到当前文件夹。
#   cd /mnt/c/Users/L/Desktop/openfoam/robot_real/leg_dynamics/results/papers
#   bash 下载_OA论文.sh
# 已存在的文件跳过；下不来的会列在最后（多半是出版社挡了脚本，点 论文索引.html 里的链接手动存即可）。
# 🔒 付费的四篇（Qu 2025、Chen 2019、Kim & Gharib 2011、Herrera-Amaya 2024）不在这里，走 eaccess.ub.tum.de。
set -u
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126 Safari/537.36"
FAIL=()

get() {   # get <输出文件名> <URL>
    local out="$1" url="$2"
    if [ -s "$out" ]; then echo "  已有，跳过  $out"; return 0; fi
    echo "  下载 $out …"
    if curl -fsSL --connect-timeout 20 --max-time 180 -A "$UA" -o "$out.part" "$url" \
       && [ -s "$out.part" ] && head -c 5 "$out.part" | grep -q '%PDF'; then
        mv "$out.part" "$out"; echo "        ok  $(du -h "$out" | cut -f1)"
    else
        rm -f "$out.part"; FAIL+=("$out  ←  $url")
    fi
}

echo "=== 第一梯队 ==="
get 01_Fish1984_muskrat.pdf                     "https://www.wcupa.edu/sciences-mathematics/biology/fFish/documents/1984JEBMuskratMechanics.pdf"
get 02_Wang2025_quadruped_paddling_CFD.pdf      "https://www.mdpi.com/2313-7673/10/3/148/pdf"

echo "=== 第二梯队 ==="
get 04_Lin2024_tendon_coupling.pdf              "https://arxiv.org/pdf/2409.14707"
get 05_Csillag2025_duck_paddling.pdf            "https://journals.biologists.com/jeb/article-pdf/doi/10.1242/jeb.249274/3717449/jeb249274.pdf"
get 06_Ribak2023_duck_feet_oars.pdf             "https://www.nature.com/articles/s41598-023-42784-w.pdf"
get 07_Han2025_learn_to_swim.pdf                "https://arxiv.org/pdf/2505.03146"

echo "=== 第三梯队（方法章） ==="
get 08_Chen2022_beaver_leg.pdf                  "https://ms.copernicus.org/articles/13/831/2022/ms-13-831-2022.pdf"
get 09_Hu2024_passive_joint_paddle.pdf          "https://www.mdpi.com/2313-7673/9/1/56/pdf"
get 10a_ITTC_resistance_test.pdf                "https://ittc.info/media/1217/75-02-02-01.pdf"
get 10b_ITTC_CFD_guidelines.pdf                 "https://www.ittc.info/media/11960/75-03-02-04.pdf"
get 11_Windt2020_overset_evaluation.pdf         "https://link.springer.com/content/pdf/10.1007/s40722-019-00156-5.pdf"

echo "=== 第四梯队（背景） ==="
get 12_Picardi2023_underwater_legged_review.pdf "https://iopscience.iop.org/article/10.1088/1748-3190/acc0bb/pdf"
get 13_Walker2000_rowing_vs_flying.pdf          "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1690750/pdf/11052540.pdf"

echo "=== 附加（有用但不急） ==="
get A1_Ding2025_diving_beetle_robot.pdf         "https://www.mdpi.com/2313-7673/10/3/182/pdf"
get A2_Qi2021_diving_beetle_kinematics.pdf      "https://www.nature.com/articles/s41598-021-96158-1.pdf"
get A3_Terziev2022_scale_effects.pdf            "https://momchil-terziev.github.io/files/Draft_1.pdf"
get A4_Stern2001_VV_methodology.pdf             "http://servidor.demec.ufpr.br/CFD/bibliografia/erros_numericos/Stern_et_al_2001.pdf"
get A5_Blake1980_recovery_stroke.pdf            "https://journals.biologists.com/jeb/article-pdf/85/1/337/3193575/jexbio_85_1_337.pdf"

echo
echo "完成：$(ls -1 *.pdf 2>/dev/null | wc -l) 个 PDF"
if [ "${#FAIL[@]}" -gt 0 ]; then
    echo
    echo "以下没下来（出版社挡脚本，用浏览器打开 论文索引.html 点链接手动存）："
    printf '  %s\n' "${FAIL[@]}"
fi
