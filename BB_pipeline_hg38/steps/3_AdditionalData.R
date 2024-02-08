#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(Battenberg))

Generate_SV_file=function(SV_BREAKPOINTS_FILE,OUTPUT_DIR,TUMOURNAME,NORMALNAME) {
  print('Generating SV file...')
  if (!file.exists(SV_BREAKPOINTS_FILE)) {
    warning('### SV file does not exist! ###')
    return(F)
  }
  # stopifnot(file.exists(SV_BREAKPOINTS_FILE))
  SV=suppressWarnings(readVcf(SV_BREAKPOINTS_FILE))
  SV=SV[which(rowRanges(SV)$FILTER=='PASS'),]
  SV=SV[which(qual(SV)>=1000 & info(SV)$AS>0 & info(SV)$RAS>0),] # high quality calls from https://github.com/PapenfussLab/gridss#quality-score
  mySV=data.frame(chr=as.character(seqnames(SV)),
                  pos=start(SV),
                  LOC2=as.character(rowRanges(SV)$ALT),
                  stringsAsFactors=F)
  stopifnot(length(grep('^.*(\\[|\\])(.*)(\\[|\\]).*$',mySV$LOC2))==nrow(mySV))
  LOC2=strsplit(gsub('^.*(\\[|\\])(.*)(\\[|\\]).*$','\\2',mySV$LOC2),':')
  mySV$chr2=sapply(LOC2,function(x) x[1])
  mySV$pos2=as.numeric(sapply(LOC2,function(x) x[2]))
  mySV$LOC2=NULL
  rm(LOC2)
  mySV$chr=gsub('^chr','',mySV$chr)
  mySV$chr2=gsub('^chr','',mySV$chr2)
  
  mySV=mySV[which(mySV$chr %in% c(1:22,'X') & mySV$chr2 %in% c(1:22,'X')),]
  SV=data.frame(chromosome=factor(c(mySV$chr,mySV$chr2),levels=c(1:22,'X')),
                position=c(mySV$pos,mySV$pos2),
                stringsAsFactors=F)
  SV=SV[order(SV$chromosome,SV$position),]
  SV=unique(SV)
  rm(mySV)
  
  LINE=1
  while(LINE<nrow(SV)) {
    if (SV$chromosome[LINE]==SV$chromosome[LINE+1] && (SV$position[LINE+1]-SV$position[LINE])<=1e3) {
      SV=SV[-(LINE+1),]
    } else {
      LINE=LINE+1
    }
  }; rm(LINE)
  SV$chromosome=as.character(SV$chromosome)
  write.table(SV,file=paste0(OUTPUT_DIR,"/Additional_data/",TUMOURNAME,"_SV.txt"),
              sep='\t',row.names=F,col.names=T,quote=F)
}

Generate_SNV_file=function(VCFFILEPATH,OUTPUT_DIR,TUMOURNAME) {
  print('Generating SNV file...')
  LOCI_FILE=paste0(OUTPUT_DIR,'/Additional_data/',TUMOURNAME,'_loci.txt')
  AF_FILE=paste0(OUTPUT_DIR,'/Additional_data/',TUMOURNAME,'_AF.txt')
  VCF=readVcf(VCFFILEPATH)
  VCF=VCF[which(rowRanges(VCF)$FILTER=='PASS')]
  VCF=VCF[which(as.character(seqnames(VCF)) %in% paste0('chr',1:22)),]
  VCF=VCF[which(as.character(ref(VCF)) %in% c('A','C','G','T')),]
  ALT=as.character(unlist(alt(VCF)))
  stopifnot(identical(length(ALT),length(VCF)))
  VCF=VCF[which(ALT %in% c('A','C','G','T')),]
  rm(ALT)
  IDs=paste0(as.character(seqnames(VCF)),'_',start(VCF))
  MULTI_ALLELIC_SNV=IDs[duplicated(IDs)]
  if (length(MULTI_ALLELIC_SNV)>=1) {
    VCF=VCF[-which(IDs %in% MULTI_ALLELIC_SNV)]
  }
  rm(MULTI_ALLELIC_SNV,IDs)
  AD=do.call(rbind,geno(VCF)$AD[,2])
  VCF=data.frame(Chromosome=gsub('^chr','',as.character(seqnames(VCF))),
                 Position=start(VCF),
                 Ref=ref(VCF),
                 Alt=as.character(unlist(alt(VCF))),
                 Ref_C=AD[,1],
                 Alt_C=AD[,2],
                 stringsAsFactors=F)
  VCF=VCF[which(VCF$Alt_C>=2),]
  AF=VCF[,1:2]
  AF$count_A=0
  AF$count_C=0
  AF$count_G=0
  AF$count_T=0
  AF$total_depth=apply(AD[AD[,2]>=2,],1,sum)
  MAPPING=c(3:6)
  names(MAPPING)=c('A','C','G','T')
  AF[cbind(1:nrow(AF),MAPPING[VCF$Ref])]=VCF$Ref_C
  AF[cbind(1:nrow(AF),MAPPING[VCF$Alt])]=VCF$Alt_C
  INDEX=which(AF$total_depth>=10)
  AF=AF[INDEX,]
  VCF=VCF[INDEX,]
  write.table(AF,file=AF_FILE,sep='\t',quote=F,row.names=F)
  write.table(VCF[,1:4],file=LOCI_FILE,sep='\t',col.names=F,row.names=F,quote=F)
}

Compute_coverage=function(OUTPUT_DIR,TUMOURNAME,NORMALNAME,PROBLEMLOCI) {
  print('Computing coverage...')
  AC=as.data.frame(Battenberg::read_table_generic(paste0(OUTPUT_DIR,'/B-RunBAFLogR/',TUMOURNAME,'_alleleCounts.tab')))
  coverage_plot(TUMOURNAME,AC,paste0(OUTPUT_DIR,'/B-RunBAFLogR/',TUMOURNAME,'_alleleCounts.png'))
  rownames(AC)=paste0(AC$Chromosome,'_',AC$Position)
  PROBLOCI=as.data.frame(Battenberg::read_table_generic(CONFIG$PROBLEMLOCI))
  PROBLOCI=paste0(PROBLOCI$Chr,'_',PROBLOCI$Pos)
  AC=AC[is.na(match(rownames(AC),PROBLOCI)),]
  writeLines(as.character(round(mean(AC$mutCountT1+AC$mutCountT2,na.rm=T),4)),con=paste0(OUTPUT_DIR,'/Additional_data/',TUMOURNAME,'_coverage.txt'))
  writeLines(as.character(round(mean(AC$mutCountN1+AC$mutCountN2,na.rm=T),4)),con=paste0(OUTPUT_DIR,'/Additional_data/',NORMALNAME,'_coverage.txt'))
}

Compute_MAPD=function(OUTPUT_DIR,TUMOURNAME) {
  logR=as.data.frame(Battenberg:::read_logr(paste0(OUTPUT_DIR,'/C-RunGCcorrect/',TUMOURNAME,'_mutantLogR_gcCorrected.tab')))[,3]
  MAPD=median(abs(diff(logR)))
  writeLines(as.character(round(MAPD,6)),con=paste0(OUTPUT_DIR,'/Additional_data/',TUMOURNAME,'_MAPD.txt'))
}

ARGS=commandArgs(trailingOnly=TRUE)
CONFIG=ARGS[1]
stopifnot(file.exists(CONFIG))
CONFIG=strsplit(readLines(CONFIG),'=')
names(CONFIG)=sapply(CONFIG,function(x) x[1])
CONFIG=lapply(CONFIG,function(x) x[2])

Generate_SV_file(SV_BREAKPOINTS_FILE = CONFIG$SV_BREAKPOINTS_FILE,
                 OUTPUT_DIR = CONFIG$OUTPUT_DIR,
                 TUMOURNAME = CONFIG$TUMOURNAME,
                 NORMALNAME = CONFIG$NORMALNAME)

Generate_SNV_file(VCFFILEPATH = CONFIG$VCFFILEPATH,
                  OUTPUT_DIR = CONFIG$OUTPUT_DIR,
                  TUMOURNAME = CONFIG$TUMOURNAME)

Compute_coverage(OUTPUT_DIR = CONFIG$OUTPUT_DIR,
                 TUMOURNAME = CONFIG$TUMOURNAME,
                 NORMALNAME = CONFIG$NORMALNAME,
                 PROBLEMLOCI = CONFIG$PROBLEMLOCI)

Compute_MAPD(OUTPUT_DIR = CONFIG$OUTPUT_DIR,
             TUMOURNAME = CONFIG$TUMOURNAME)