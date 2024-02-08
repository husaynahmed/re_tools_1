#!/usr/bin/env bash
CONFIG=$1
source ${CONFIG} || exit 99

find ${OUTPUT_DIR}/A-GetAlleleCounts/ -type f -name "*_alleleFrequencies_chr*.txt" -exec gzip -9 -f {} \; || exit 1
if [[ -e ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_alleleCounts.tab ]]; then gzip -9 ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_alleleCounts.tab || exit 2; fi
if [[ -e ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_normalLogR.tab ]]; then gzip -9 ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_normalLogR.tab || exit 3; fi
if [[ -e ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_normalBAF.tab ]]; then gzip -9 ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_normalBAF.tab || exit 4; fi
if [[ -e ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_mutantLogR.tab ]]; then gzip -9 ${OUTPUT_DIR}/B-RunBAFLogR/${TUMOURNAME}_mutantLogR.tab || exit 5; fi
find ${OUTPUT_DIR}/D-GenerateImputeInputFromAlleleFrequencies/ -type f -name "*_impute_input_chr*.txt" -exec gzip -9 -f {} \; || exit 6
find ${OUTPUT_DIR}/E-RunBeagle5/ -type f -name "*_beagle5_input_chr*.vcf" -exec gzip -9 -f {} \; || exit 7
find ${OUTPUT_DIR}/E-RunBeagle5/ -type f -name "*_impute_output_chr*_allHaplotypeInfo.txt" -exec gzip -9 -f {} \; || exit 8
find ${OUTPUT_DIR}/G-GetHaplotypedBAFs/ -type f -name "*_chr*_heterozygousMutBAFs_haplotyped.txt" -exec gzip -9 -f {} \; || exit 9

exit 0
