#!/usr/bin/env Rscript

MYcallChrXsubclones = function(tumourname,X_gamma=1000,X_kmin=100,genomebuild,AR=TRUE,prior_breakpoints_file=NULL,chrom_names){
  
  print(tumourname)
  
  PCFinput=data.frame(read_table_generic(paste0(tumourname,"_mutantLogR_gcCorrected.tab")),stringsAsFactors=F)
  ChrNotation=unique(PCFinput[which(!is.na(match(PCFinput$Chromosome,c("X","chrX")))),]$Chromosome) # find the chromosome notation
  PCFinput=PCFinput[which(PCFinput$Chromosome==ChrNotation & PCFinput$Position>2.6e6 & PCFinput$Position<156e6),] # get nonPAR
  colnames(PCFinput)[3]=tumourname
  print(paste("Number of chrX nonPAR SNPs =",nrow(PCFinput)))
  
  if (!is.null(prior_breakpoints_file)) {
    sv=read.table(prior_breakpoints_file, header=T, stringsAsFactors=F)
    sv=sv[which(!is.na(match(sv$chr,c("X","chrX")))),]
    breakpoints=c(min(PCFinput$Position),sv$pos,max(PCFinput$Position))
    PCF=data.frame()
    for (j in 1:(length(breakpoints)-1)) {
      PCFinput_sv=PCFinput[which(PCFinput$Position>=breakpoints[j] & PCFinput$Position<breakpoints[j+1]),]
      if (nrow(PCFinput_sv)==0) next # Fix 1 - case where there is no SNP between two SVs
      PCF_sv=copynumber::pcf(PCFinput_sv,gamma=X_gamma,kmin=X_kmin)
      PCF=rbind(PCF,PCF_sv)
    }
  } else {
    PCF=copynumber::pcf(PCFinput,gamma=X_gamma,kmin=X_kmin)
  }
  write.table(PCF,paste0(tumourname,"_PCF_gamma_",X_gamma,"_chrX.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  print("PCF segmentation done")
  
  if (genomebuild=="hg19"){
    x_centromere=c(58632012,61632012) # hg19
    ar=data.frame(startpos=66763874,endpos=66950461)
  } else {
    x_centromere=c(58605580,62412542) #hg38
    ar=data.frame(startpos=67544021,endpos=67730619)
  }
  
  # INPUT for copy number inference
  SAMPLEsegs=data.frame(PCF,stringsAsFactors=F)
  pupl=read.table(paste0(tumourname,"_purity_ploidy.txt"),header=T,stringsAsFactors=F)
  SAMPLEpurity=pupl[,1] # SAMPLEpurity=pupl$cellularity in previous Battenberg version; change from pupl$purity to pupl[,1] for universality
  #SAMPLEwgd=ifelse(round(pupl$ploidy/2)*2==4,T,F)
  SAMPLEn=pupl$ploidy
  print(paste(SAMPLEpurity,SAMPLEn))
  
  # Estimating LogR deviation in diploid and gained regions (AUTOSOMAL)
  BB=read.table(paste0(tumourname,"_copynumber_extended.txt"),header=T,stringsAsFactors = F)
  
  BBdip=BB[which(BB$nMaj1_A==1 & BB$nMin1_A==1 & BB$frac1_A==1),]
  # correction for LogR values
  BBcorr=-mean(BBdip$LogR) #diploid only
  if (nrow(BBdip)<=1){
    print("likely WGD sample")
    BBdip=BB[which(BB$nMaj1_A==2 & BB$nMin1_A==2 & BB$frac1_A==1),]
    cnloh=BB[which(BB$nMaj1_A==2 & BB$nMin1_A==0 & BB$frac1_A==1),]
    if (nrow(cnloh)>0){
      BBcorr=-mean(cnloh$LogR)
    } else if (nrow(cnloh)==0){
      print("CRUDE estimation of BBcorr based on assumption of 2 copies vs ploidy")
      BBcorr=-log2(2/SAMPLEn)
    }
  }
  BBg1=BB[which(BB$nMaj1_A==2 & BB$nMin1_A==1 & BB$frac1_A==1),]
  BBg2=BB[which(BB$nMaj1_A==3 & BB$nMin1_A==1 & BB$frac1_A==1),]
  BBg3=BB[which(BB$nMaj1_A==4 & BB$nMin1_A==1 & BB$frac1_A==1),]
  BBg4=BB[which(BB$nMaj1_A==3 & BB$nMin1_A==2 & BB$frac1_A==1),] # likely observed in WGD samples
  
  # get max gain N:
  BBcomb=rbind(BBdip,BBg1,BBg2,BBg3,BBg4)
  maxNMaj=max(BBcomb$nMaj1_A)
  
  # SD for LogR values - diploid and gain regions
  BBsd=c(sd(BBdip$LogR),sd(BBg1$LogR),sd(BBg2$LogR),sd(BBg3$LogR))
  #BBsd_mean=mean(BBsd,na.rm=T)
  BBsd_max=max(BBsd, na.rm=T)
  BBsd_max=max(BBsd_max,0.05) # accept a minimum of 5% sd in LogR variation
  
  # BB LOH - estimating sd for LOH/loss events
  BBloh=BB[which(BB$nMaj1_A==1 & BB$nMin1_A==0 & BB$frac1_A==1),]
  if (nrow(BBloh)<=1){ #sd would be NA
    print("likely WGD sample or no clonal LOH event")
    BBloh=BB[which(BB$nMin1_A==0 & BB$frac1_A==1),] # all LOH events with varying nMaj1_A including 2:0 events
  }
  
  # expected ChrX logR values
  explogrgainX=function(x){log2((SAMPLEpurity*x+(1-SAMPLEpurity)*1)/1)}
  explogrGain=sapply(2:10000,explogrgainX) # up to 10000 copies!
  
  explogrLoss=max(log2(0+(1-SAMPLEpurity)*1),log2(0.01)) # if purity ~ 1, then purity of 0.99 is assumed for a realistic explogR estimate
  
  # assign CN
  SEG=data.frame()
  for (j in 1:nrow(SAMPLEsegs)){
    seg=SAMPLEsegs[j,]
    seg$type=ifelse(seg$mean<0,"loss","gain")
    
    # is segment different from zero?
    seg$mean=seg$mean+BBcorr
    
    if (seg$type=="gain"){
      seg$CNA=ifelse(seg$mean>(0+1.96*BBsd_max),"yes","no")
    } else {
      seg$CNA=ifelse(seg$mean<(0-1.96*BBsd_max),"yes","no")
    }
    # copy number
    if (seg$CNA=="yes"){
      if (seg$type=="gain"){
        rank=which(sort(c(explogrGain,seg$mean))==seg$mean) # rank of observed logR mean for segment among the expected logR values
        seg$CN=rank+1
        # clonality test
        if (rank==1){
          seg$clonal=ifelse(round(explogrGain[rank]-seg$mean,digits=2)<=round((BBsd_max/explogrGain[rank]),digits=2),"yes","no") # CV
          
        } else if (rank>=5){ # STOPS calling 'subclonal' events when copy number is >=5
          if (abs(seg$mean-explogrGain[rank-1])<abs(seg$mean-explogrGain[rank])){
            
            seg$clonal="yes"
            seg$CN=seg$CN-1
          } else {
            seg$clonal="yes"
          }
        } else {
          if (abs(seg$mean-explogrGain[rank-1])<abs(seg$mean-explogrGain[rank])){
            
            seg$clonal=ifelse(round(seg$mean-explogrGain[rank-1],digits = 2)<=round((BBsd_max/explogrGain[rank-1]),digits=2),"yes","no")
            if (seg$clonal=="yes"){
              seg$CN=seg$CN-1
            }
          }
          else {
            seg$clonal=ifelse(round(explogrGain[rank]-seg$mean,digits=2)<round((BBsd_max/explogrGain[rank]),digits = 2),"yes","no")
          }
        }
      } else if (seg$type=="loss"){
        seg$CN=0
        if (nrow(BBloh)>0){
          seg$clonal=ifelse(round(abs(explogrLoss-seg$mean),digits=2)<round(abs(sd(BBloh$LogR)/explogrLoss),digits=2),"yes","no")
        } else if (nrow(BBloh)==0){
          seg$clonal=ifelse(round(abs(explogrLoss-seg$mean),digits=2)<round(abs(BBsd_max/explogrLoss),digits=2),"yes","no")
        }
      }
    }
    else {
      seg$CN=1
      seg$clonal=NA
      print(paste("no CNA for segment",j))
    }
    if (seg$arm=="p" & seg$end.pos>x_centromere[1]-1e6 & seg$CNA=="yes" & seg$end.pos<seg$start.pos+1e6){
      print("segment is p-arm centromere noise")
      print(seg)
    } else if (seg$arm=="q" & seg$end.pos<x_centromere[2]+1e6 & seg$CNA=="yes" & seg$end.pos<seg$start.pos+1e6){
      print("segment is q-arm centromere noise")
      print(seg)
    } else {
      SEG=rbind(SEG,seg)
    }
  }
  
  # CALCULATE CCF 
  CCF=data.frame()
  for (j in 1:nrow(SEG)){
    seg=SEG[j,]
    if (seg$CNA=="yes"){
      if (seg$type=="gain"){
        if (seg$clonal=="no"){
          seg$CCF=(2^seg$mean-(SAMPLEpurity*(seg$CN-1)+(1-SAMPLEpurity)*1))/SAMPLEpurity 
        } else {
          seg$CCF=1
        }
      } else if (seg$type=="loss"){
        if (seg$clonal=="no"){
          seg$CCF=(1-2^(seg$mean))/SAMPLEpurity # for Loss (assuming one chrX in all cells prior to Loss)
          if (seg$CCF>=0.95){
            seg$CCF=1
            seg$clonal="yes"
          }
        } else {
          seg$CCF=1
        }
      }
    } else {
      seg$CCF=1
    }
    CCF=rbind(CCF,seg)
  }
  
  # GENERATE FINAL OUTPUT
  SUBCLONES=data.frame()
  for (j in 1:nrow(CCF)){
    subclones=CCF[j,]
    if (subclones$CNA=="no"){
      subclones=data.frame(subclones,nMaj1=1,nMin1=0,frac1=1,nMaj2=0,nMin2=0,frac2=0)
    } else {
      if (subclones$type=="gain" & subclones$clonal=="yes"){
        subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=1,nMaj2=0,nMin2=0,frac2=0)  
      }
      else if (subclones$type=="gain" & subclones$clonal=="no"){
        if(subclones$CCF>0.5){ # switch nMaj/nMin so that the first nMaj/nMin represent the MAJOR CLONE
          subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=subclones$CCF,nMaj2=subclones$CN-1,nMin2=0,frac2=1-subclones$CCF)
        } else {
          subclones=data.frame(subclones,nMaj1=subclones$CN-1,nMin1=0,frac1=1-subclones$CCF,nMaj2=subclones$CN,nMin2=0,frac2=subclones$CCF)  
        }
      }
      else if (subclones$type=="loss" & subclones$clonal=="yes"){
        subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=1,nMaj2=0,nMin2=0,frac2=0) # very unlikely scenario; no sequencing reads should be present!
      }
      else if (subclones$type=="loss" & subclones$clonal=="no"){
        if(subclones$CCF>0.5){ # switch nMaj/nMin so that the first nMaj/nMin represent the MAJOR CLONE
          subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=subclones$CCF,nMaj2=1,nMin2=0,frac2=1-subclones$CCF)
        } else {
          subclones=data.frame(subclones,nMaj1=1,nMin1=0,frac1=1-subclones$CCF,nMaj2=subclones$CN,nMin2=0,frac2=subclones$CCF)
        }
      }
    }
    print(j)
    SUBCLONES=rbind(SUBCLONES,subclones)
  }
  
  SUBCLONES$average=(SUBCLONES$nMaj1+SUBCLONES$nMin1)*SUBCLONES$frac1+(SUBCLONES$nMaj2+SUBCLONES$nMin2)*SUBCLONES$frac2
  
  SUBCLONESout=data.frame(SUBCLONES[,c("chrom","arm")],startpos=SUBCLONES$start.pos,endpos=SUBCLONES$end.pos,nSNPs=SUBCLONES$n.probes,
                          LogR=SUBCLONES$mean,SUBCLONES[,c("type","CNA","CN","clonal","nMaj1","nMin1","frac1","nMaj2","nMin2","frac2")],
                          subclonalCN=SUBCLONES$average,stringsAsFactors = F)
  SUBCLONESout$type[SUBCLONESout$type=="gain"]="+ve"
  SUBCLONESout$type[SUBCLONESout$type=="loss"]="-ve"
  
  # merge adjacent segments with same copy number  
  SUBCLONESout$rank=1:nrow(SUBCLONESout)
  SUBCLONESout=SUBCLONESout[order(SUBCLONESout$subclonalCN),]
  
  SPLIT=split(SUBCLONESout$rank, cumsum(c(1, diff(SUBCLONESout$rank) != 1))) # find consecutive segments with same subclonalCN
  outputDF=data.frame()
  for (j in 1:length(SPLIT)){
    if (length(SPLIT[[j]])>1){
      print(length(SPLIT[[j]]))
      SUBsplit=SUBCLONESout[which(!is.na(match(SUBCLONESout$rank,SPLIT[[j]]))),]
      if (length(unique(SUBsplit$arm))==1){
        if (sd(SUBsplit$subclonalCN)<=0.01){
          mergedseg=SUBsplit[1,]
          mergedseg$endpos=SUBsplit[length(SPLIT[[j]]),"endpos"]
          mergedseg$nSNPs=sum(SUBsplit$nSNPs)
          mergedseg$LogR=weighted.mean(SUBsplit$LogR,SUBsplit$nSNPs)
          outputDF=rbind(outputDF,mergedseg)
        } else {
          outputDF=rbind(outputDF,SUBsplit)
          print("adjacent not same subclonalCN in SPLIT")
        }
      } else if (length(SPLIT[[j]])==2){
        outputDF=rbind(outputDF,SUBsplit)
      } else{
        # if (length(SUBsplit$arm=="p"))
        pseg=SUBsplit[SUBsplit$arm=="p",]
        if (nrow(pseg)>1){
          if (sd(pseg$subclonalCN)<=0.01){
            mergedseg=pseg[1,]
            mergedseg$endpos=pseg[nrow(pseg),"endpos"]
            mergedseg$nSNPs=sum(pseg$nSNPs)
            mergedseg$LogR=weighted.mean(pseg$LogR,pseg$nSNPs)
            outputDF=rbind(outputDF,mergedseg)
          } else {
            outputDF=rbind(outputDF,pseg)
            print("adjacent not same subclonalCN in pseg")
          }
        } else {outputDF=rbind(outputDF,pseg)}
        qseg=SUBsplit[SUBsplit$arm=="q",]
        if (nrow(qseg)>1){
          if (sd(qseg$subclonalCN)<=0.01){
            mergedseg=qseg[1,]
            mergedseg$endpos=qseg[nrow(qseg),"endpos"]
            mergedseg$nSNPs=sum(qseg$nSNPs)
            mergedseg$LogR=weighted.mean(qseg$LogR,qseg$nSNPs)
            outputDF=rbind(outputDF,mergedseg)
          } else {
            outputDF=rbind(outputDF,qseg)
            print("adjacent not same subclonalCN in qseg")
          }
        } else {outputDF=rbind(outputDF,qseg)}
      }
    }  else {
      SUBsplit=SUBCLONESout[which(SUBCLONESout$rank==SPLIT[[j]]),]
      outputDF=rbind(outputDF,SUBsplit)
    }
  }
  outputDF=outputDF[order(outputDF$startpos),]
  outputDF$rank=NULL
  
  print(paste("Number of rows merged =",nrow(SUBCLONESout)-nrow(outputDF)))
  
  BBnew=BB[which(is.na(match(BB$chr,c("X","chrX")))),c(1:3,8:13)] # copynumber.txt columns to be populated with chrX calls
  
  outputDF_for_merge=data.frame(chr=outputDF$chrom,startpos=outputDF$startpos,endpos=outputDF$endpos,
                                nMaj1_A=outputDF$nMaj1,nMin1_A=outputDF$nMin1,frac1_A=outputDF$frac1,
                                nMaj2_A=outputDF$nMaj2,nMin2_A=outputDF$nMin2,frac2_A=outputDF$frac2,
                                stringsAsFactors = F)
  
  BBnew=rbind(BBnew,outputDF_for_merge)
  write.table(BBnew,paste0(tumourname,"_copynumber.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
  BBnew_extended=BB[which(is.na(match(BB$chr,c("X","chrX")))),] # copynumber_extended.txt columns for chrX
  
  outputDF_for_merge_extended=data.frame(chr=outputDF$chrom,startpos=outputDF$startpos,endpos=outputDF$endpos,BAF=NA,pval=NA,LogR=outputDF$LogR,ntot=outputDF$CN,
                                         nMaj1_A=outputDF$nMaj1,nMin1_A=outputDF$nMin1,frac1_A=outputDF$frac1,nMaj2_A=outputDF$nMaj2,nMin2_A=outputDF$nMin2,
                                         frac2_A=outputDF$frac2)
  BtoFsolutions=data.frame(matrix(nrow= nrow(outputDF),ncol = ncol(BB)-ncol(outputDF_for_merge_extended)))
  names(BtoFsolutions)=names(BB)[(ncol(outputDF_for_merge_extended)+1):ncol(BB)]
  
  BBnew_extended=rbind(BBnew_extended,cbind(outputDF_for_merge_extended,BtoFsolutions))
  write.table(BBnew_extended,paste0(tumourname,"_copynumber_extended.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  write.table(outputDF,paste0(tumourname,"_chrX_copynumber.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
  # PLOT
  outputDF$diff=outputDF$endpos-outputDF$startpos
  PGAclonal=sum(outputDF[which(outputDF$clonal=="yes"),]$diff)/sum(outputDF[which(!is.na(outputDF$clonal)),]$diff)
  
  
  plot_BB=ggplot()+geom_hline(yintercept = 0:ceiling(max(outputDF$subclonalCN)),linetype="longdash",col="grey",size=0.2)+
    geom_rect(data=outputDF,aes(xmin=startpos,xmax=endpos,ymin=subclonalCN-0.02,ymax=subclonalCN+0.02))+
    geom_vline(xintercept = x_centromere,linetype="longdash",col="green")+
    #geom_hline(yintercept = nonpar,linetype="dotted",col="blue")+
    ylim(-0.2,ceiling(max(outputDF$subclonalCN))+0.2)+labs(x="ChrX coordinate (bp)",y="Average Ploidy")+
    theme(plot.title = element_text(hjust = 0.5,size=12),panel.background = element_blank())+
    ggtitle(paste0(tumourname," , PLOIDY: ",round(SAMPLEn,digits = 3)," , PURITY: ",round(SAMPLEpurity*100,digits = 0),
                   "%, PGAclonal: ",round(PGAclonal*100,digits = 1),"%"))
  
  # ANDROGEN RECEPTOR LOCUS
  if (AR){
    data.table::setDT(ar)
    data.table::setkey(ar,"startpos","endpos")
    data.table::setDT(outputDF)
    data.table::setkey(outputDF,"startpos","endpos")
    segAR=data.table::foverlaps(ar,outputDF,type="any",nomatch = 0)
    segAR$subclonalCN=(segAR$nMaj1+segAR$nMin1)*segAR$frac1+(segAR$nMaj2+segAR$nMin2)*segAR$frac2
    plot_BB=plot_BB+geom_rect(data=segAR,aes(xmin=startpos,xmax=endpos,ymin=subclonalCN-0.02,ymax=subclonalCN+0.02),fill="red")
  }
  
  pdf(paste0(tumourname,"_chrX_average_ploidy.pdf"))
  print(plot_BB)
  dev.off()
  
  # Update the genomewide Battenberg plots
  # goodness from rho_psi file (i.e. column named 'distance')
  goodness=read.table(paste0(tumourname,"_rho_and_psi.txt"),header=T,stringsAsFactors = F,sep="\t")
  goodness=goodness[which(goodness$is.best=="TRUE"),"distance"]
  # rho and ploidy from purity_ploidy file
  rho_psi=read.table(paste0(tumourname,"_purity_ploidy.txt"),header=T,stringsAsFactors = F,sep="\t")
  rho=rho_psi$purity # FIX 2 - this was rho_psi$cellularity
  ploidy=rho_psi$ploidy
  # Need BAFsegment file
  BAFvals=as.data.frame(Battenberg:::read_bafsegmented(paste0(tumourname,".BAFsegmented.txt")))
  BAFvals=rbind(BAFvals[which(is.na(match(BAFvals$Chromosome,c("X","chrX")))),],
                data.frame(Chromosome="X",Position=sort(sample(1:155e6,90000,replace=F)), # 155e6: approximate length of chrX
                           BAF=sample(c(0,1),90000,replace=T),BAFphased=1,BAFseg=1)) # 90000 = typical no. of het SNPs expected based on chrX length (roughly around chr 7 and 8 average hetSNP counts) 
  
  Battenberg:::plot.gw.subclonal.cn(subclones=BBnew, 
                                    BAFvals=BAFvals, 
                                    rho=rho, 
                                    ploidy=ploidy, 
                                    goodness=goodness, 
                                    output.gw.figures.prefix=paste(tumourname,"_BattenbergProfile", sep=""), 
                                    chr.names=chrom_names, 
                                    tumourname=tumourname)
}


suppressPackageStartupMessages(library(Battenberg))
suppressPackageStartupMessages(library(ggplot2))
options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
CONFIG=args[1]
stopifnot(file.exists(CONFIG))
CONFIG=strsplit(readLines(CONFIG),'=')
names(CONFIG)=sapply(CONFIG,function(x) x[1])
CONFIG=lapply(CONFIG,function(x) x[2])

setwd(paste0(CONFIG$OUTPUT_DIR,'/L-FitCopyNumber_',CONFIG$RUN))
file.copy(from=paste0(CONFIG$OUTPUT_DIR,'/K-SegmentBAFphased/',CONFIG$TUMOURNAME,'.BAFsegmented.txt'),
  to=paste0(CONFIG$OUTPUT_DIR,'/L-FitCopyNumber_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'.BAFsegmented.txt'))

print('Running: fit.copy.number')
fit.copy.number(samplename=CONFIG$TUMOURNAME,
  outputfile.prefix=paste0(CONFIG$OUTPUT_DIR,"/L-FitCopyNumber_",CONFIG$RUN,"/",CONFIG$TUMOURNAME, "_"),
  inputfile.baf.segmented=paste0(CONFIG$OUTPUT_DIR,"/L-FitCopyNumber_",CONFIG$RUN,"/",CONFIG$TUMOURNAME, ".BAFsegmented.txt"),
  inputfile.baf=paste0(CONFIG$OUTPUT_DIR,"/B-RunBAFLogR/",CONFIG$TUMOURNAME,"_mutantBAF.tab"),
  inputfile.logr=paste0(CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/",CONFIG$TUMOURNAME,"_mutantLogR_gcCorrected.tab"),
  dist_choice=as.numeric(CONFIG$CLONALITY_DIST_METRIC),
  ascat_dist_choice=as.numeric(CONFIG$ASCAT_DIST_METRIC),
  min.ploidy=as.numeric(CONFIG$MIN_PLOIDY),
  max.ploidy=as.numeric(CONFIG$MAX_PLOIDY),
  min.rho=as.numeric(CONFIG$MIN_RHO),
  max.rho=as.numeric(CONFIG$MAX_RHO),
  min.goodness=as.numeric(CONFIG$MIN_GOODNESS_OF_FIT),
  uninformative_BAF_threshold=as.numeric(CONFIG$BALANCED_THRESHOLD),
  gamma_param=as.numeric(CONFIG$PLATFORM_GAMMA),
  use_preset_rho_psi=as.logical(CONFIG$USE_PRESET_RHO_PSI),
  preset_rho=as.numeric(CONFIG$PRESET_RHO),
  preset_psi=as.numeric(CONFIG$PRESET_PSI),
  read_depth=30)
  
print('Running: callSubclones')
callSubclones(sample.name=CONFIG$TUMOURNAME,
  baf.segmented.file=paste0(CONFIG$OUTPUT_DIR,"/L-FitCopyNumber_",CONFIG$RUN,"/",CONFIG$TUMOURNAME, ".BAFsegmented.txt"),
  logr.file=paste0(CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/",CONFIG$TUMOURNAME,"_mutantLogR_gcCorrected.tab"),
  rho.psi.file=paste0(CONFIG$OUTPUT_DIR,"/L-FitCopyNumber_",CONFIG$RUN,"/",CONFIG$TUMOURNAME, "_rho_and_psi.txt"),
  output.file=paste0(CONFIG$OUTPUT_DIR,"/M-CallSubclones_",CONFIG$RUN,"/",CONFIG$TUMOURNAME,"_copynumber.txt"), # must be called *_copynumber.txt
  output.figures.prefix=paste0(CONFIG$OUTPUT_DIR,"/M-CallSubclones_",CONFIG$RUN,"/",CONFIG$TUMOURNAME,"_subclones_chr"),
  output.gw.figures.prefix=paste0(CONFIG$OUTPUT_DIR,"/M-CallSubclones_",CONFIG$RUN,"/",CONFIG$TUMOURNAME,"_BattenbergProfile"),
  chr_names=get.chrom.names(CONFIG$IMPUTEINFOFILE, as.logical(CONFIG$IS_MALE)),
  masking_output_file=paste0(CONFIG$OUTPUT_DIR,"/M-CallSubclones_",CONFIG$RUN,"/",CONFIG$TUMOURNAME, "_segment_masking_details.txt"),
  max_allowed_state=as.numeric(CONFIG$MAX_CN_STATE),
  prior_breakpoints_file=if (file.exists(paste0(CONFIG$OUTPUT_DIR,"/Additional_data/",CONFIG$TUMOURNAME,"_SV.txt"))) paste0(CONFIG$OUTPUT_DIR,"/Additional_data/",CONFIG$TUMOURNAME,"_SV.txt") else NULL,
  gamma=as.numeric(CONFIG$PLATFORM_GAMMA),
  segmentation.gamma=as.numeric(CONFIG$SEGMENTATION_GAMMA),
  siglevel=0.05,
  maxdist=0.01,
  noperms=1000,
  seed=as.numeric(CONFIG$SEED))

setwd(paste0(CONFIG$OUTPUT_DIR,'/M-CallSubclones_',CONFIG$RUN))

print('make_posthoc_plots')
make_posthoc_plots(samplename=CONFIG$TUMOURNAME,
                   logr_file=paste0(CONFIG$OUTPUT_DIR,'/C-RunGCcorrect/',CONFIG$TUMOURNAME,'_mutantLogR_gcCorrected.tab'),
                   bafsegmented_file=paste0(CONFIG$OUTPUT_DIR,'/L-FitCopyNumber_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'.BAFsegmented.txt'),
                   logrsegmented_file=paste0(CONFIG$OUTPUT_DIR,'/L-FitCopyNumber_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'.logRsegmented.txt'),
                   allelecounts_file=NULL)

print('cnfit_to_refit_suggestions')
cnfit_to_refit_suggestions(samplename=CONFIG$TUMOURNAME,
                           subclones_file=paste0(CONFIG$OUTPUT_DIR,'/M-CallSubclones_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_copynumber_extended.txt'),
                           rho_psi_file=paste0(CONFIG$OUTPUT_DIR,'/L-FitCopyNumber_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_rho_and_psi.txt'),
                           gamma_param=1)

if (as.logical(CONFIG$IS_MALE)) {
  file.copy(from=paste0(CONFIG$TUMOURNAME,"_copynumber.txt"),to=paste0(CONFIG$TUMOURNAME,"_copynumber_beforeX.txt")) # Save
  file.copy(from=paste0(CONFIG$TUMOURNAME,"_copynumber_extended.txt"),to=paste0(CONFIG$TUMOURNAME,"_copynumber_extended_beforeX.txt")) # Save
  system(paste0('ln -s ',CONFIG$OUTPUT_DIR,"/C-RunGCcorrect/",CONFIG$TUMOURNAME,"_mutantLogR_gcCorrected.tab .")) # Symlink, remove at the end
  system(paste0('ln -s ',CONFIG$OUTPUT_DIR,"/L-FitCopyNumber_",CONFIG$RUN,"/",CONFIG$TUMOURNAME,".BAFsegmented.txt .")) # Symlink, remove at the end
  system(paste0('ln -s ',CONFIG$OUTPUT_DIR,"/L-FitCopyNumber_",CONFIG$RUN,"/",CONFIG$TUMOURNAME, "_rho_and_psi.txt .")) # Symlink, remove at the end
  print('callChrXsubclones')
  MYcallChrXsubclones(tumourname=CONFIG$TUMOURNAME,
                      X_gamma=1000,
                      X_kmin=100,
                      genomebuild='hg38',
                      AR=FALSE,
                      prior_breakpoints_file=if (file.exists(paste0(CONFIG$OUTPUT_DIR,"/Additional_data/",CONFIG$TUMOURNAME,"_SV.txt"))) paste0(CONFIG$OUTPUT_DIR,"/Additional_data/",CONFIG$TUMOURNAME,"_SV.txt") else NULL,
                      chrom_names=c(1:22,'X'))
  system(paste0('rm ',CONFIG$TUMOURNAME,"_mutantLogR_gcCorrected.tab"))
  system(paste0('rm ',CONFIG$TUMOURNAME,".BAFsegmented.txt"))
  system(paste0('rm ',CONFIG$TUMOURNAME, "_rho_and_psi.txt"))
}

print('All done.')