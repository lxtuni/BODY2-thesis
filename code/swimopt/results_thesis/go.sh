#!/bin/bash
# run up to 2 jobs in parallel for one sandbox call; each job resumes from its checkpoint
cd /sessions/zealous-happy-mayer/mnt/SA/swimopt
W=${W:-135}
for job in "$@"; do
  IFS='|' read -r tag kind mode seed budget ov <<< "$job"
  MUJOCO_GL=disable timeout 170 python3 results_thesis/cma_run.py $tag $kind $mode $seed $budget $W $ov >> results_thesis/log_$tag.txt 2>&1 &
done
wait
for job in "$@"; do IFS='|' read -r tag rest <<< "$job"; tail -1 results_thesis/log_$tag.txt; done
