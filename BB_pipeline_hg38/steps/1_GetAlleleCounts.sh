#!/usr/bin/env bash
CONFIG=$1
source ${CONFIG} || exit 99

if [[ $SLURM_JOB_ID != "" ]]; then
  JOB_ID=$SLURM_JOB_ID
  TASK_ID=$SLURM_ARRAY_TASK_ID
elif [[ $LSB_JOBID != "" ]]; then
  JOB_ID=$LSB_JOBID
  TASK_ID=$LSB_JOBINDEX
else
  echo "Unknown scheduler" && exit 99
fi

if [[ $TASK_ID -gt 23 ]]; then
  SAMPLE_NAME=$NORMALNAME
  BAM=$NORMALBAM
  CHR="$(($TASK_ID-23))"
else
  SAMPLE_NAME=$TUMOURNAME
  BAM=$TUMOURBAM
  CHR=$TASK_ID
fi

CHR_NAME=$CHR
if [[ $CHR -eq 23 ]]; then CHR_NAME='X'; fi

if [[ "${BAM#*.}" == "cram" ]]; then
  ADD="-r ${REFERENCE_FASTA}"
else
  ADD=""
fi

${ALLELECOUNTER} \
  -b ${BAM} \
  -l ${G1000_PREFIX_AC}${CHR_NAME}.txt \
  -o ${OUTPUT_DIR}/A-GetAlleleCounts/${SAMPLE_NAME}_alleleFrequencies_chr${CHR_NAME}.txt \
  ${ADD} -m 20 -q 35 --dense-snps || exit 1

sed -i "s+chr${CHR_NAME}+${CHR_NAME}+g" ${OUTPUT_DIR}/A-GetAlleleCounts/${SAMPLE_NAME}_alleleFrequencies_chr${CHR_NAME}.txt || exit 2

exit 0
