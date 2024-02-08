#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Battenberg))
suppressPackageStartupMessages(library(ASCAT))
suppressPackageStartupMessages(library(data.table))

myCuratedPlot=function(TumourLogR,TumourBAF,NormalLogR,NormalBAF,probloci,SampleName) {
  CheckGZ=function(FILE) {
    if (!file.exists(FILE)) {return(paste0(FILE,'.gz'))} else {return(FILE)}
  }
  read_logR_BAF_probloci=function(PATH,isprobloci=F) {
    print(paste0('   Reading: ',PATH))
    stopifnot(file.exists(PATH) & file.info(PATH)$size>0)
    tmp=fread(file=PATH,data.table=F,sep='\t',showProgress=F,header=T)
    print(paste0('   nrows=',nrow(tmp)))
    if (isprobloci) {
      return(paste0(tmp[,1],'_',tmp[,2]))
    } else {
      rownames(tmp)=paste0(tmp[,1],'_',tmp[,2])
      return(tmp)
    }
  }
  TumourLogR=CheckGZ(TumourLogR)
  TumourBAF=CheckGZ(TumourBAF)
  NormalLogR=CheckGZ(NormalLogR)
  NormalBAF=CheckGZ(NormalBAF)
  stopifnot(all(file.exists(c(TumourLogR,TumourBAF,NormalLogR,NormalBAF,probloci))))
  TumourLogR=read_logR_BAF_probloci(TumourLogR)
  TumourBAF=read_logR_BAF_probloci(TumourBAF)
  NormalLogR=read_logR_BAF_probloci(NormalLogR)
  NormalBAF=read_logR_BAF_probloci(NormalBAF)
  probloci=read_logR_BAF_probloci(probloci,isprobloci=T)
  TO_KEEP=setdiff(Reduce(intersect,list(rownames(TumourLogR),rownames(TumourBAF),rownames(NormalLogR),rownames(NormalBAF))),probloci)
  DATA=cbind(TumourLogR[TO_KEEP,],
             TumourBAF[TO_KEEP,3],
             NormalLogR[TO_KEEP,3],
             NormalBAF[TO_KEEP,3])
  colnames(DATA)=c('Chromosome','Position','Tumour_LogR','Tumour_BAF','Normal_LogR','Normal_BAF')
  rm(TumourLogR,TumourBAF,NormalLogR,NormalBAF,probloci,TO_KEEP)
  save(DATA,file='SNPpos.Rdata')

  ch=list()
  chr_names=c(1:22,'X')
  for (i in 1:length(chr_names)) {
    temp=which(DATA$Chromosome==chr_names[i])
    if (length(temp)==0) {
      ch[[i]]=0
    } else {
      ch[[i]]=temp[1]:temp[length(temp)]
    }
  }; rm(i)
  ascat.bc=list(Tumor_LogR=as.data.frame(DATA[,'Tumour_LogR']),
                Tumor_BAF=as.data.frame(DATA[,'Tumour_BAF']),
                Germline_LogR=as.data.frame(DATA[,'Normal_LogR']),
                Germline_BAF=as.data.frame(DATA[,'Normal_BAF']),
                Tumor_LogR_segmented=NULL,
                Tumor_BAF_segmented=NULL,
                SNPpos=DATA[,1:2],
                chrs=chr_names,
                samples=SampleName,
                chrom=Battenberg:::split_genome(DATA[,1:2]),
                ch=ch)
  ascat.plotRawData(ascat.bc)
}

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
CONFIG=args[1]
stopifnot(file.exists(CONFIG))
CONFIG=strsplit(readLines(CONFIG),'=')
names(CONFIG)=sapply(CONFIG,function(x) x[1])
CONFIG=lapply(CONFIG,function(x) x[2])

setwd(paste0(CONFIG$OUTPUT_DIR,'/B-RunBAFLogR'))

print('Running: getBAFsAndLogRs')
getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste0(CONFIG$OUTPUT_DIR,"/A-GetAlleleCounts/",CONFIG$TUMOURNAME,"_alleleFrequencies_chr"),
  normalAlleleCountsFile.prefix=paste0(CONFIG$OUTPUT_DIR,"/A-GetAlleleCounts/",CONFIG$NORMALNAME,"_alleleFrequencies_chr"),
  figuresFile.prefix=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME, "_"),
  BAFnormalFile=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_normalBAF.tab"),
  BAFmutantFile=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_mutantBAF.tab"),
  logRnormalFile=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_normalLogR.tab"),
  logRmutantFile=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_mutantLogR.tab"),
  combinedAlleleCountsFile=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_alleleCounts.tab"),
  chr_names=get.chrom.names(CONFIG$IMPUTEINFOFILE, as.logical(CONFIG$IS_MALE)),
  g1000file.prefix=CONFIG$G1000_PREFIX,
  minCounts=as.numeric(CONFIG$MIN_NORMAL_DEPTH),
  samplename=CONFIG$TUMOURNAME,
  seed=as.numeric(CONFIG$SEED))

print('Running: gc.correct.wgs')
gc.correct.wgs(Tumour_LogR_file=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_mutantLogR.tab"),
  outfile=paste0(CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/",CONFIG$TUMOURNAME,"_mutantLogR_gcCorrected.tab"),
  correlations_outfile=paste0(CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/",CONFIG$TUMOURNAME, "_GCwindowCorrelations.txt"),
  gc_content_file_prefix=CONFIG$GCCORRECTPREFIX,
  replic_timing_file_prefix=CONFIG$RTCORRECTPREFIX,
  chrom_names=get.chrom.names(CONFIG$IMPUTEINFOFILE, as.logical(CONFIG$IS_MALE)),
  recalc_corr_afterwards=T)

print('myCuratedPlot')
setwd(paste0(CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/"))
myCuratedPlot(TumourLogR=paste0(CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/",CONFIG$TUMOURNAME,"_mutantLogR_gcCorrected.tab"),
              TumourBAF=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_mutantBAF.tab"),
              NormalLogR=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_normalLogR.tab"),
              NormalBAF=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_normalBAF.tab"),
              probloci=CONFIG$PROBLEMLOCI,
              SampleName=CONFIG$TUMOURNAME)

print('All done.')