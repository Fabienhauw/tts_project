#!/bin/bash
 cmd=$( printf "j%02d_dcm_estimate" ${SLURM_ARRAY_TASK_ID})
 bash /network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/scripts/dcm/dcm_estimate/$cmd


 echo seff -d ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} >> do_seff
