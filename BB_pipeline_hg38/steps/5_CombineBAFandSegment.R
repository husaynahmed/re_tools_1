#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Battenberg))

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
CONFIG=args[1]
stopifnot(file.exists(CONFIG))
CONFIG=strsplit(readLines(CONFIG),'=')
names(CONFIG)=sapply(CONFIG,function(x) x[1])
CONFIG=lapply(CONFIG,function(x) x[2])

setwd(paste0(CONFIG$OUTPUT_DIR,'/K-SegmentBAFphased'))

print('Running: combine.baf.files')
combine.baf.files(inputfile.prefix=paste0(CONFIG$OUTPUT_DIR,"/G-GetHaplotypedBAFs/",CONFIG$TUMOURNAME, "_chr"),
  inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
  outputfile=paste0(CONFIG$OUTPUT_DIR,"/J-CombineBAFfiles/",CONFIG$TUMOURNAME, "_heterozygousMutBAFs_haplotyped.txt"),
  chr_names=get.chrom.names(CONFIG$IMPUTEINFOFILE, as.logical(CONFIG$IS_MALE)))

print('Running: segment.baf.phased')
segment.baf.phased(samplename=CONFIG$TUMOURNAME,
  inputfile=paste0(CONFIG$OUTPUT_DIR,"/J-CombineBAFfiles/",CONFIG$TUMOURNAME, "_heterozygousMutBAFs_haplotyped.txt"),
  outputfile=paste0(CONFIG$OUTPUT_DIR,"/K-SegmentBAFphased/",CONFIG$TUMOURNAME, ".BAFsegmented.txt"),
  prior_breakpoints_file=if (file.exists(paste0(CONFIG$OUTPUT_DIR,"/Additional_data/",CONFIG$TUMOURNAME,"_SV.txt"))) paste0(CONFIG$OUTPUT_DIR,"/Additional_data/",CONFIG$TUMOURNAME,"_SV.txt") else NULL,
  gamma=as.numeric(CONFIG$SEGMENTATION_GAMMA),
  phasegamma=as.numeric(CONFIG$PHASING_GAMMA))
  
print('All done.')
