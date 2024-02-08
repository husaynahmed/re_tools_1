#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dpclust3p))
suppressPackageStartupMessages(library(DPClust))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(CNAqc))

GenerateDPinputt2BamAC=function(OUTPUT_DIR,TUMOURNAME,IS_MALE,RUN) {
  LOCI_FILE=paste0(OUTPUT_DIR,'/Additional_data/',TUMOURNAME,'_loci.txt')
  AF_FILE=paste0(OUTPUT_DIR,'/Additional_data/',TUMOURNAME,'_AF.txt')
  cellularity_file=paste0(OUTPUT_DIR,'/L-FitCopyNumber_',RUN,'/',TUMOURNAME,'_rho_and_psi.txt')
  subclone_file=paste0(OUTPUT_DIR,'/M-CallSubclones_',RUN,'/',TUMOURNAME,'_copynumber_extended.txt')
  output_file=paste0(OUTPUT_DIR,'/Refit_',RUN,'/',TUMOURNAME,'_DPinput.txt')
  for (i in c(LOCI_FILE,AF_FILE,cellularity_file,subclone_file)) {
    stopifnot(file.exists(i))
  }; rm(i)
  
  if (file.exists(output_file) && file.exists(gsub('_DPinput.txt','_DPinput_original.txt',output_file))) {
    print('Using pre-existing DPinput file')
  } else {
    runGetDirichletProcessInfo(loci_file=LOCI_FILE,
                               allele_frequencies_file=AF_FILE, 
                               cellularity_file=cellularity_file,
                               subclone_file=subclone_file,
                               gender=if (IS_MALE=='TRUE' || IS_MALE==T) 'male' else 'female',
                               SNP.phase.file="NA", 
                               mut.phase.file="NA",
                               output_file=output_file)
    
    file.copy(output_file,gsub('_DPinput.txt','_DPinput_original.txt',output_file))
    DP.info=read.table(output_file,sep="\t",header=T,stringsAsFactors = F,check.names=F)
    DP.info=DP.info[!is.na(DP.info$subclonal.fraction) & DP.info$subclonal.fraction>0, ]
    set.seed(1234)
    if (nrow(DP.info)>75e3) DP.info=DP.info[sort(sample(1:nrow(DP.info),75e3,replace=F)),]
    write.table(DP.info,output_file,sep='\t',row.names=F,col.names=T,quote=F,append=F)
  }
}

RunDPClust=function(OUTPUT_DIR,PARTICIPANTID,TUMOURNAME,NORMALNAME,RUN) {
  OUT_DIR=paste0('/Refit_',RUN,'/')
  dpinfofilepath=paste0(OUTPUT_DIR,OUT_DIR,TUMOURNAME,'_DPinput.txt')
  rhoandpsifilepath=paste0(OUTPUT_DIR,'/L-FitCopyNumber_',RUN,'/',TUMOURNAME,'_rho_and_psi.txt')

  cellularity=read.table(rhoandpsifilepath,header=T,stringsAsFactors=F,sep="\t")[2,1]
  DP.info=read.table(dpinfofilepath,sep="\t",header=T,stringsAsFactors = F)
  mutCount=array(DP.info$mut.count,c(nrow(DP.info),1))
  WTCount=array(DP.info$WT.count,c(nrow(DP.info),1))
  totalCopyNumber=array(DP.info$subclonal.CN,c(nrow(DP.info),1))
  copyNumberAdjustment=array(DP.info$no.chrs.bearing.mut,c(nrow(DP.info),1))
  mutation.copy.number=array(DP.info$mutation.copy.number,c(nrow(DP.info),1))
  no.iters=1000
  DPClust:::DirichletProcessClustering(mutCount=mutCount,
                             WTCount=WTCount,
                             totalCopyNumber=totalCopyNumber,
                             copyNumberAdjustment=copyNumberAdjustment,
                             mutation.copy.number=mutation.copy.number,
                             cellularity=cellularity,
                             output_folder=paste0(OUTPUT_DIR,OUT_DIR),
                             no.iters=no.iters,
                             no.iters.burn.in=ceiling(no.iters/5),
                             subsamplesrun=TUMOURNAME,
                             samplename=TUMOURNAME,
                             conc_param=1,
                             cluster_conc=5,
                             mut.assignment.type=1,
                             most.similar.mut=NA,
                             mutationTypes=NA,
                             max.considered.clusters=30)
}

AssessBBDPCRun=function (sampledir,participantid,tumourplatekey,normalplatekey,RUN,lowerclonal,upperclonal) {
  dpcdir=paste0('/Refit_',RUN,'/')
  subfilename = paste0(sampledir, "/M-CallSubclones_",RUN,"/", tumourplatekey, "_copynumber_extended.txt")
  purityfilename = paste0(sampledir, "/M-CallSubclones_",RUN,"/", tumourplatekey, "_purity_ploidy.txt")

  tumournormal = tumourplatekey
  dpcfilename = paste0(sampledir, dpcdir, tumournormal, "_optimaInfo.txt")
  if (file.exists(dpcfilename)) {
    optima = read.table(dpcfilename, header = T, stringsAsFactors = F)
    dpcmetrics = matrix(NA, nrow = 0, ncol = 18)
    number.mutations = sum(as.numeric(optima[, 3]))
    number.clusters = nrow(optima)
    toremove = which(optima[, 3] <= number.mutations/100)
    if (length(toremove) > 0) {
      optima = optima[-toremove, ]
    }
    number.mutations.post1pcremoval = sum(as.numeric(optima[, 3]))
    number.clusters.post1pcremoval = nrow(optima)
    dists = abs(optima[, 2] - (upperclonal + lowerclonal)/2)
    position.clonalcluster = optima[which(dists == min(dists)), 2]
    numbermutations.clonalcluster = optima[which(dists == min(dists)), 3]
    position.topccfcluster = optima[nrow(optima), 2]
    numbermutations.topccfcluster = optima[nrow(optima), 3]
    most.muts = which(optima[, 3] == max(optima[, 3]))
    if (length(most.muts) > 1) {
      position.clustermostmutations = optima[most.muts[length(most.muts)], 2]
    } else {
      position.clustermostmutations = optima[which(optima[, 3] == max(optima[, 3])), 2]
    }
    if (length(most.muts) > 1) {
      numbermutations.clustermostmutations = optima[most.muts[length(most.muts)], 3]
    } else {
      numbermutations.clustermostmutations = optima[which(optima[, 3] == max(optima[, 3])), 3]
    }
    superclonal = which(optima[, 2] > upperclonal)
    if (length(superclonal) > 0) {
      position.superclonalcluster = optima[max(superclonal), 2]
      numbermutations.superclonalcluster = optima[max(superclonal), 3]
    } else {
      position.superclonalcluster = 0
      numbermutations.superclonalcluster = 0
    }
    subclonal = which(optima[, 2] < lowerclonal)
    if (length(subclonal) == 0) {
      position.biggestsubclonalcluster = 0
      numbermutations.biggestsubclonalcluster = 0
    } else if (length(subclonal) == 1) {
      position.biggestsubclonalcluster = optima[which(optima[subclonal, 3] == max(optima[subclonal, 3])), 2]
      numbermutations.biggestsubclonalcluster = optima[which(optima[subclonal, 3] == max(optima[subclonal, 3])), 3]
    } else if (length(subclonal) > 1) {
      subclonal.clusters = which(optima[subclonal, 3] == max(optima[subclonal, 3]))
      position.biggestsubclonalcluster = optima[subclonal.clusters[length(subclonal.clusters)], 2]
      numbermutations.biggestsubclonalcluster = optima[subclonal.clusters[length(subclonal.clusters)], 3]
    }
    if ((optima[which(dists == min(dists)), 3] == max(optima[, 3]) | optima[which(dists == min(dists)), 3] == optima[nrow(optima), 3]) & (position.clonalcluster >= lowerclonal & position.clonalcluster <= upperclonal)) {
      passeddpc = "Yes"
    } else {
      passeddpc = "No"
    }
    dpcprepfilename = paste0(sampledir, dpcdir, tumourplatekey, "_DPinput.txt")
    dpcprep = read.table(dpcprepfilename, header = T, stringsAsFactors = F)
    print(paste("Read dpcprepfile", dpcprepfilename))
    no.muts.chr1 = length(which(dpcprep[, 16] == 1))
    no.muts.chr2 = length(which(dpcprep[, 16] == 2))
    no.muts.chrmore1 = length(which(dpcprep[, 16] > 1))
    upperbound = (position.clonalcluster/2) + 0.05
    lowerbound = (position.clonalcluster/2) - 0.05
    if (length(which(optima$location < upperbound & optima$location > lowerbound)) == 0) {
      cluster.at.50 = FALSE
      position.cluster.at.50 = FALSE
      number.muts.cluster.at.50 = FALSE
    } else {
      cluster.at.50 = TRUE
      position.cluster.at.50 = optima[which(optima$location < upperbound & optima$location > lowerbound), 2]
      number.muts.cluster.at.50 = optima[which(optima$location < upperbound & optima$location > lowerbound), 3]
      # FIX tlesluyes - one sample had two possible solutions here
      if (length(position.cluster.at.50>1)) {
        position.cluster.at.50=position.cluster.at.50[which.max(number.muts.cluster.at.50)]
        number.muts.cluster.at.50=number.muts.cluster.at.50[which.max(number.muts.cluster.at.50)]
      }
    }
    dpcmetrics = cbind(number.mutations, number.mutations.post1pcremoval, 
                       number.clusters, number.clusters.post1pcremoval, 
                       position.clonalcluster, numbermutations.clonalcluster, 
                       position.topccfcluster, numbermutations.topccfcluster, 
                       position.clustermostmutations, numbermutations.clustermostmutations, 
                       position.superclonalcluster, numbermutations.superclonalcluster, 
                       position.biggestsubclonalcluster, numbermutations.biggestsubclonalcluster, 
                       no.muts.chr1, no.muts.chr2, no.muts.chrmore1, cluster.at.50, 
                       position.cluster.at.50, number.muts.cluster.at.50, 
                       passeddpc)
    print(paste("Collected DPClust metrics"))
  } else {
    dpcmetricscolnames = c("number.mutations", "number.mutations.post1pcremoval", 
                           "number.clusters", "number.clusters.post1pcremoval", 
                           "position.clonalcluster", "numbermutations.clonalcluster", 
                           "position.topccfcluster", "numbermutations.topccfcluster", 
                           "position.clustermostmutations", "numbermutations.clustermostmutations", 
                           "position.superclonalcluster", "numbermutations.superclonalcluster", 
                           "position.biggestsubclonalcluster", "numbermutations.biggestsubclonalcluster", 
                           "no.muts.chr1", "no.muts.chr2", "no.muts.chrmore1", 
                           "cluster.at.50", "position.cluster.at.50", "number.muts.cluster.at.50", 
                           "passeddpc")
    dpcmetrics = matrix(NA, 1, 21)
    colnames(dpcmetrics) = dpcmetricscolnames
    print(paste("No DPClust runs completed"))
  }
  if (file.exists(subfilename)) {
    sub = read.table(subfilename, header = T, stringsAsFactors = F)
    purityploidy = read.table(purityfilename, header = T, stringsAsFactors = F)
    purity = purityploidy[1, 1]
    ploidy = purityploidy[1, 3]
    segments = nrow(sub)
    genome.length = as.numeric(sub[, 3]) - as.numeric(sub[, 2])
    genome.length = sum(as.numeric(genome.length))
    clonallydiploid = sub[which(sub[, 8] == 1 & sub[, 9] == 1 & sub[, 10] == 1), ]
    segments.clonallydiploid = nrow(clonallydiploid)
    length.clonallydiploid = sum(as.numeric(clonallydiploid[, 3]) - as.numeric(clonallydiploid[, 2]))
    subclonallydiploid = sub[which((sub[, 8] == 1 & sub[, 9] == 1 & sub[, 10] != 1) | sub[, 11] == 1 & sub[, 12] == 1), ]
    segments.subclonallydiploid = nrow(subclonallydiploid)
    length.subclonallydiploid = sum(as.numeric(subclonallydiploid[, 3]) - as.numeric(subclonallydiploid[, 2]))
    clonallytetraploid = sub[which(sub[, 8] == 2 & sub[, 9] == 2 & sub[, 10] == 1), ]
    segments.clonallytetraploid = nrow(clonallytetraploid)
    length.clonallytetraploid = sum(as.numeric(clonallytetraploid[, 3]) - as.numeric(clonallytetraploid[, 2]))
    subclonallytetraploid = sub[which((sub[, 8] == 2 & sub[, 9] == 2 & sub[, 10] != 1) | sub[, 11] == 2 & sub[, 12] == 2), ]
    segments.subclonallytetraploid = nrow(subclonallytetraploid)
    length.subclonallytetraploid = sum(as.numeric(subclonallytetraploid[, 3]) - as.numeric(subclonallytetraploid[, 2]))
    clonalaberration = sub[which(sub$frac1_A == 1 & (sub$nMaj1_A != 1 | sub$nMin1_A != 1)), ]
    segments.clonalaberration = nrow(clonalaberration)
    length.clonalaberration = sum(as.numeric(clonalaberration[, 3] - clonalaberration[, 2]))
    subclonalaberration = sub[which(sub$frac1_A != 1 & !is.na(sub$frac1_A)), ]
    segments.subclonalaberration = nrow(subclonalaberration)
    length.subclonalaberration = sum(as.numeric(subclonalaberration[, 3] - subclonalaberration[, 2]))
    clonalLOH = sub[which(sub$frac1_A == 1 & (sub$nMaj1_A == 0 | sub$nMin1_A == 0)), ]
    segments.clonalLOH = nrow(clonalLOH)
    length.clonalLOH = sum(as.numeric(clonalLOH[, 3] - clonalLOH[, 2]))
    subclonalLOH = sub[which(sub$frac1_A != 1 & (sub$nMaj1_A == 0 | sub$nMin1_A == 0 | sub$nMaj2_A == 0 | sub$nMin2_A == 0)), ]
    segments.subclonalLOH = nrow(subclonalLOH)
    length.subclonalLOH = sum(as.numeric(subclonalLOH[, 3] - subclonalLOH[, 2]))
    clonalodd = sub[which(sub$frac1_A == 1 & (sub$nMaj1_A%%2 != 0 | sub$nMin1_A%%2 != 0)), ]
    segments.clonalodd = nrow(clonalodd)
    length.clonalodd = sum(as.numeric(clonalodd[, 3] - clonalodd[, 2]))
    midpoint = (upperclonal + lowerclonal)/2/2
    fiftylower = midpoint - 0.025
    fiftyupper = midpoint + 0.025
    fiftypcsubclones = sub[which(sub$frac1_A < fiftyupper & sub$frac1_A > fiftylower), ]
    segments.fiftypcsubclones = nrow(fiftypcsubclones)
    length.fiftypcsubclones = sum(as.numeric(fiftypcsubclones[, 3] - fiftypcsubclones[, 2]))
    clonalhomdels = sub[which(sub$frac1_A == 1 & sub$nMaj1_A == 0 & sub$nMin1_A == 0), ]
    segments.clonalhomdels = nrow(clonalhomdels)
    length.clonalhomdels = sum(as.numeric(clonalhomdels[, 3] - clonalhomdels[, 2]))
    if (segments.clonalhomdels > 0) {
      maxlength.clonalhomdels = max(clonalhomdels[, 3] - clonalhomdels[, 2])
    } else {
      maxlength.clonalhomdels = 0
    }
    subclonalhomdels = sub[which(sub$frac1_A != 1 & ((sub$nMaj1_A == 0 & sub$nMin1_A == 0) | (sub$nMaj2_A == 0 & sub$nMin2_A == 0))), ]
    segments.subclonalhomdels = nrow(subclonalhomdels)
    length.subclonalhomdels = sum(as.numeric(subclonalhomdels[, 3] - subclonalhomdels[, 2]))
    if (segments.subclonalhomdels > 0) {
      maxlength.subclonalhomdels = max(subclonalhomdels[, 3] - subclonalhomdels[, 2])
    } else {
      maxlength.subclonalhomdels = 0
    }
    bbmetrics = cbind(purity, ploidy, segments, genome.length, 
                      segments.clonallydiploid, length.clonallydiploid, 
                      segments.subclonallydiploid, length.subclonallydiploid, 
                      segments.clonallytetraploid, length.clonallytetraploid, 
                      segments.subclonallytetraploid, length.subclonallytetraploid, 
                      segments.clonalaberration, length.clonalaberration, 
                      segments.subclonalaberration, length.subclonalaberration, 
                      segments.clonalLOH, length.clonalLOH, segments.subclonalLOH, 
                      length.subclonalLOH, segments.clonalodd, length.clonalodd, 
                      segments.fiftypcsubclones, length.fiftypcsubclones, 
                      segments.clonalhomdels, length.clonalhomdels, maxlength.clonalhomdels, 
                      segments.subclonalhomdels, length.subclonalhomdels, 
                      maxlength.subclonalhomdels)
  } else {
    bbmetricscolnames = c("purity", "ploidy", "segments", 
                          "genome.length", "segments.clonallydiploid", "length.clonallydiploid", 
                          "segments.subclonallydiploid", "length.subclonallydiploid", 
                          "segments.clonallytetraploid", "length.clonallytetraploid", 
                          "segments.subclonallytetraploid", "length.subclonallytetraploid", 
                          "segments.clonalaberration", "length.clonalaberration", 
                          "segments.subclonalaberration", "length.subclonalaberration", 
                          "segments.clonalLOH", "length.clonalLOH", "segments.subclonalLOH", 
                          "length.subclonalLOH", "segments.clonalodd", "length.clonalodd", 
                          "segments.fiftypcsubclones", "length.fiftypcsubclones", 
                          "segments.clonalhomdels", "length.clonalhomdels", 
                          "maxlength.clonalhomdels", "segments.subclonalhomdels", 
                          "length.subclonalhomdels", "maxlength.subclonalhomdels")
    bbmetrics = matrix(NA, 1, 30)
    colnames(bbmetrics) = bbmetricscolnames
    print(paste("No BB runs completed"))
  }
  samplemetrics = as.data.frame(cbind(participantid, tumourplatekey, normalplatekey, bbmetrics, dpcmetrics))
  write.csv(samplemetrics, paste0(sampledir, dpcdir, tumourplatekey, "_metrics_run", RUN), row.names = F, quote = F)
  if (samplemetrics$passeddpc == "Yes") {
    write.csv("PASS", paste0(sampledir, dpcdir, "PASS"), quote = F, row.names = F)
  } else if (samplemetrics$passeddpc == "No") {
    write.csv("FAIL", paste0(sampledir, dpcdir, "FAIL"), quote = F, row.names = F)
  }
  print(paste("Written final metrics and indicator pass/fail file for", tumourplatekey))
  return(samplemetrics)
}

RunCNAqc=function (sampledir,TUMOURNAME,RUN) {
  read.snvs = function (AF_FILE,LOCI_FILE) {
    AF=read.table(AF_FILE,sep='\t',header=T,stringsAsFactors=F)
    loci=read.table(LOCI_FILE,sep='\t',header=F,stringsAsFactors=F)
    stopifnot(identical(paste0(AF[,1],':',AF[,2]),paste0(loci[,1],':',loci[,2])))
    AF$REF=AF[cbind(1:nrow(loci),unlist(list('A'=3,'C'=4,'G'=5,'T'=6)[loci$V3]))]
    AF$ALT=AF[cbind(1:nrow(loci),unlist(list('A'=3,'C'=4,'G'=5,'T'=6)[loci$V4]))]
    AF$VAF=AF$ALT/(AF$REF+AF$ALT)
    IS_NOT_NA=which(!is.na(AF$VAF))
    OUT=data.frame(chr=paste0('chr',AF$Chromosome),
                   from=AF$Position,
                   to=AF$Position,
                   ref=loci$V3,
                   alt=loci$V4,
                   DP=AF$REF+AF$ALT,
                   NV=AF$ALT,
                   VAF=AF$VAF,
                   stringsAsFactors=F)[IS_NOT_NA,]
    return(OUT)
  }
  read.segs = function (filename.segs) {
    cna=read.table(filename.segs,sep='\t',header=T,stringsAsFactors=F)[,1:13]
    cna$chr=paste0('chr',cna$chr)
    cna=cna[which(!is.na(cna$nMaj1_A)),]
    cna$Major=cna$nMaj1_A
    cna$minor=cna$nMin1_A
    cna$CCF=cna$frac1_A
    cna$Major_2=cna$nMaj2_A
    cna$minor_2=cna$nMin2_A
    INDEX=which(!is.na(cna$frac2_A) & cna$frac2_A>cna$frac1_A)
    if (length(INDEX)>=1) {
      cna$Major[INDEX]=cna$nMaj2_A[INDEX]
      cna$minor[INDEX]=cna$nMin2_A[INDEX]
      cna$CCF[INDEX]=cna$frac2_A[INDEX]
      cna$Major_2[INDEX]=cna$nMaj1_A[INDEX]
      cna$minor_2[INDEX]=cna$nMin1_A[INDEX]
    }
    rm(INDEX)
    stopifnot(min(cna$CCF)>=0.5)
    cna=cna[,c(1:3,14:18)]
    colnames(cna)=c('chr','from','to','Major','minor','CCF','Major_2','minor_2')
    return(cna)
  }

  AF_FILE = paste0(sampledir,'/Additional_data/',TUMOURNAME,"_AF.txt")
  LOCI_FILE = paste0(sampledir,'/Additional_data/',TUMOURNAME,"_loci.txt")
  filename.segs = paste0(sampledir,"/M-CallSubclones_",RUN,"/",TUMOURNAME,"_copynumber_extended.txt")
  filename.purity.ploidy = paste0(sampledir,"/M-CallSubclones_",RUN,"/",TUMOURNAME,"_purity_ploidy.txt")
  filename.output = paste0(sampledir,"/Refit_",RUN,"/",TUMOURNAME,"_peaks_output.Rdata")
  filename.output.figure = paste0(sampledir,"/Refit_",RUN,"/",TUMOURNAME,"_peaks_output.png")
  for (i in c(AF_FILE,LOCI_FILE,filename.segs,filename.purity.ploidy)) {
    stopifnot(file.exists(i))
  }; rm(i)
  mutations=read.snvs(AF_FILE,LOCI_FILE)
  cna=read.segs(filename.segs)
  purity=read.table(filename.purity.ploidy, sep = "\t", header = T, stringsAsFactors = F)[1,1]
  data=init(mutations=mutations,cna=cna,purity=purity,ref='GRCh38')
  tryCatch(
    expr={
      peaks=analyze_peaks(data,
                          karyotypes = c('1:0', '1:1', '2:1', '2:0', '2:2'),
                          min_karyotype_size = 0.05,
                          min_absolute_karyotype_mutations = 50,
                          p_binsize_peaks = 0.005,
                          purity_error = 0.025,
                          n_bootstrap = 1,
                          kernel_adjust = 1,
                          matching_strategy = "closest")
      if (is.null(peaks$peaks_analysis)) peaks$peaks_analysis$QC='FLAG'
      
      save(peaks,file=filename.output)
      if (peaks$peaks_analysis$QC=="FLAG") {write.csv("FLAGPEAKS",paste0(sampledir,"/Refit_",RUN,"/FLAGPEAKS"),quote=F,row.names=F)}
      if (peaks$peaks_analysis$QC=="PASS") {write.csv("PASSPEAKS",paste0(sampledir,"/Refit_",RUN,"/PASSPEAKS"),quote=F,row.names=F)}
      if (peaks$peaks_analysis$QC=="FAIL") {write.csv("FAILPEAKS",paste0(sampledir,"/Refit_",RUN,"/FAILPEAKS"),quote=F,row.names=F)}
      
      png(filename.output.figure,width=25,height=15,units='cm',res=300,pointsize=6)
      plot(plot_peaks_analysis(peaks))
      dev.off()
    },
    error=function(cond) {
      message(cond)
      write.csv("FLAGPEAKS",paste0(sampledir,"/Refit_",RUN,"/FLAGPEAKS"),quote=F,row.names=F)
    }
  )
}

reestimate_ploidy=function(rho_old,rho_new,psi_old,WGD) {
  # From Anna's bbdpcrecall package
  COEF=if(WGD) 4 else 2
  return(((rho_old * psi_old) + COEF * (rho_new - rho_old))/rho_new)
}

GenerateIDcard=function(CONFIG) {
  assessWGD=function(ploidy,fracHomo) {
    return(ploidy>fracHomo*-2+2.9)
  }
  ploidy_from_BB=function(subclones) {
    subclones=subclones[which(subclones$chr %in% 1:22),]
    stopifnot(nrow(subclones)>0)
    seg_length = floor((subclones$endpos-subclones$startpos)/1000)
    is_subclonal_maj = abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
    is_subclonal_min = abs(subclones$nMin1_A - subclones$nMin2_A) > 0
    is_subclonal_maj[is.na(is_subclonal_maj)] = F
    is_subclonal_min[is.na(is_subclonal_min)] = F
    segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1)  + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0) 
    segment_states_maj = subclones$nMaj1_A * ifelse(is_subclonal_maj, subclones$frac1_A, 1)  + ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0) 
    ploidy = sum((segment_states_min+segment_states_maj) * seg_length, na.rm=T) / sum(seg_length, na.rm=T)
    return(ploidy)
  }
  getModeAllele=function(cn,col) {
    y=round(cn[,col])
    y[y>=5]=5
    y=tapply(1:nrow(cn),y,function(z) sum((cn[z,"endpos"]-cn[z,"startpos"])/1000000))
    ord=order(y,decreasing=T)
    y=y[ord]
    return(as.numeric(names(y)[which.max(y)]))
  }
  myPlotImage=function(myImage,YLIM=c(0,1),XLIM=c(0,1),TEXT=NULL,COL='black') {
    par_mar=par()$mar
    par(mar=rep(0,4))
    plot(NULL,ylim=c(0,1),xlim=c(0,1),axes=F,xlab='',ylab='')
    rasterImage(myImage,xleft=XLIM[1],ybottom=YLIM[1],xright=XLIM[2],ytop=YLIM[2])
    if (!is.null(TEXT)) text(0.5,0.95,TEXT,col=COL,cex=2)
    par(mar=par_mar)
  }
  plot_CNseg_mutCN=function(BB,DPinput) {
    COLS=c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999') # Set1 from RColorBrewer
    CHRsize=data.frame(chr=1:22,
                       size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566),
                       middle=c(124625310.5,370850307.5,591461209,786049562,972084330,1148099493.5,1313226358.5,1465977701,1609766427.5,1748140516.5,1883411148,2017840353.5,2142351240,2253610949,2358551415,2454994487.5,2540769469,2620405698,2689008813.5,2750086065,2805663772.5,2855381003),
                       cumul=c(249250621,492449994,690472424,881626700,1062541960,1233657027,1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,2199936179,2307285719,2409817111,2500171864,2581367074,2659444322,2718573305,2781598825,2829728720,2881033286),
                       add=c(0,249250621,492449994,690472424,881626700,1062541960,1233657027,1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,2199936179,2307285719,2409817111,2500171864,2581367074,2659444322,2718573305,2781598825,2829728720)) # Hg19 data!
    Clust=table(DPinput$cluster)
    while(length(Clust)>length(COLS)) {
      TO_RM=which.min(Clust)
      DPinput=DPinput[which(DPinput$cluster!=names(Clust[TO_RM])),]
      Clust=Clust[-TO_RM]
      rm(TO_RM)
    }
    for (i in 2:22) {
      INDEX=which(BB$chr==i)
      if (length(INDEX)>0) {
        BB$startpos[INDEX]=BB$startpos[INDEX]+CHRsize$add[i]
        BB$endpos[INDEX]=BB$endpos[INDEX]+CHRsize$add[i]
      }
      INDEX=which(DPinput$chr==i)
      if (length(INDEX)>0) {
        DPinput$start[INDEX]=DPinput$start[INDEX]+CHRsize$add[i]
        DPinput$end[INDEX]=DPinput$end[INDEX]+CHRsize$add[i]
      }
    }; rm(i)
    PAR_MAR=par()$mar
    par(mar=c(2.1,4.2,3.15,1.05))
    plot(NULL,ylim=c(0,5),xlim=c(0,CHRsize$cumul[22]),xlab='',xaxs='i',xaxt='n',ylab='Major CN',cex.lab=1.25)
    points(DPinput$start,DPinput$mutation.copy.number,pch=16,cex=0.5,col=COLS[DPinput$cluster])
    segments(BB$startpos,BB$unrounded_major_cn,BB$endpos,BB$unrounded_major_cn,lwd=2,col='black')
    abline(v=CHRsize$add[2:22],col='grey25')
    axis(1,at=CHRsize$middle,1:22,las=2)
    legend('top',paste0('C',names(Clust),': ',Clust),col=COLS[as.numeric(names(Clust))],xpd=T,horiz=T,pch=16,x.intersp=0.5,inset=c(0,-0.2))
    par(mar=PAR_MAR)
  }
  
  FILES=c(paste0(CONFIG$OUTPUT_DIR,'/Refit_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_DPinput.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/Refit_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_DP_and_cluster_info.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/M-CallSubclones_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_copynumber_extended.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/M-CallSubclones_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_BattenbergProfile_average.png'),
          paste0(CONFIG$OUTPUT_DIR,'/L-FitCopyNumber_1/',CONFIG$TUMOURNAME,'_distance.png'),
          paste0(CONFIG$OUTPUT_DIR,'/Refit_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_DirichletProcessplot_with_cluster_locations.png'),
          paste0(CONFIG$OUTPUT_DIR,'/Additional_data/',CONFIG$TUMOURNAME,'_SV.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/Additional_data/',CONFIG$TUMOURNAME,'_coverage.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/Additional_data/',CONFIG$NORMALNAME,'_coverage.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/Additional_data/',CONFIG$TUMOURNAME,'_MAPD.txt'),
          paste0(CONFIG$OUTPUT_DIR,'/M-CallSubclones_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_purity_ploidy.txt'))
  names(FILES)=c('DPinput','DPclusters','BB','BBprofile','Sunrise','DPplot','SV','CoverageT','CoverageN','MAPD','CelPlo')
  for (i in 1:length(FILES)) {
    if (i!=7) stopifnot(file.exists(FILES[i]))
  }
  CONFIG$DIAG=if ('DIAG' %in% names(CONFIG)) gsub('\'|\"','',CONFIG$DIAG) else 'TBD'
  
  DPinput=read.table(FILES['DPinput'],sep='\t',header=T,stringsAsFactors=F)
  DPinput=DPinput[!is.na(DPinput$subclonal.fraction) & DPinput$subclonal.fraction>0,]
  DPclusters=read.table(FILES['DPclusters'],sep='\t',header=T,stringsAsFactors=F)
  stopifnot(nrow(DPinput)==nrow(DPclusters))
  DPinput$cluster=DPclusters$most.likely.cluster
  rm(DPclusters)
  DPinput=DPinput[which(DPinput$cluster %in% names(which(table(DPinput$cluster)>10))),]
  
  CovT=as.numeric(readLines(FILES['CoverageT']))
  CovN=as.numeric(readLines(FILES['CoverageN']))
  MAPD=as.numeric(readLines(FILES['MAPD']))
  BBmetrics=read.table(FILES['CelPlo'],stringsAsFactors=F,sep='\t',header=T)
  NRPCC=(BBmetrics$purity*CovT)/(BBmetrics$purity*BBmetrics$ploidy+(1-BBmetrics$purity)*2)
  BB=read.table(FILES['BB'],header=T,stringsAsFactors=F,sep='\t')[,c(1:3,7,8:13)]
  BB=BB[which(BB$chr %in% 1:22),]
  BB=BB[which(!is.na(BB$nMaj1_A)),]
  INDEX=which(is.na(BB$frac2_A))
  if (length(INDEX)>0) {
    BB$nMaj2_A[INDEX]=BB$nMin2_A[INDEX]=BB$frac2_A[INDEX]=0
  }; rm(INDEX)
  BB$unrounded_major_cn=BB$nMaj1_A*BB$frac1_A+BB$nMaj2_A*BB$frac2_A
  BB$unrounded_minor_cn=BB$nMin1_A*BB$frac1_A+BB$nMin2_A*BB$frac2_A
  BB$major_cn=round(BB$unrounded_major_cn)
  BB$minor_cn=round(BB$unrounded_minor_cn)
  BB$width=BB$endpos-BB$startpos+1
  medianCNDiff=median(BB$major_cn-BB$minor_cn)
  Size50pCCF=sum(BB$width[which(BB$frac1_A>0.45 & BB$frac1_A<0.55)])
  fracSubclonalCNS=sum(BB$width[which(BB$frac1_A!=1)])/sum(BB$width)
  modeMajor=getModeAllele(BB,'major_cn')
  modeMinor=getModeAllele(BB,'minor_cn')
  HD_size=sum(BB$width[which(BB$nMaj1_A==0 & BB$nMin1_A==0 & BB$frac1_A>0.2)])
  HD_fraction=HD_size/sum(BB$width)*100
  ploidy=ploidy_from_BB(BB)
  fracHomo=sum(BB$width[which(BB$minor_cn==0)])/sum(BB$width)
  WGD=assessWGD(ploidy,fracHomo)
  if (WGD) {
    GI=1-sum(BB$width[which(BB$major_cn==2 & BB$minor_cn==2)])/sum(BB$width)
  } else {
    GI=1-sum(BB$width[which(BB$major_cn==1 & BB$minor_cn==1)])/sum(BB$width)
  }
  
  COMBS=c()
  for (i in 0:4) {
    for (j in 0:4) {
      if (j<i) next
      COMBS=c(COMBS,sum(BB$width[which(BB$major_cn==j & BB$minor_cn==i)])/sum(BB$width)*100)
      names(COMBS)[length(COMBS)]=paste0(j,'+',i)
    }
  }
  COMBS=c(COMBS,100-sum(COMBS))
  names(COMBS)[length(COMBS)]='Others'
  
  BBprofile=readPNG(FILES['BBprofile'])
  Sunrise=readPNG(FILES['Sunrise'])
  DPplot=readPNG(FILES['DPplot'])
  
  pdf(paste0(CONFIG$OUTPUT_DIR,'/Refit_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_IDcard.pdf'),pointsize=8,width=7+7/3)
  layout(matrix(c(1,1,1,2,2,3,4,5,6,7,7,8),ncol=3,byrow=T),heights=c(0.1,0.25,0.4,0.25),widths=c(rep(1/3,3)))
  par(mar=c(0,0,0,0))
  plot(NULL,ylim=c(0,1),xlim=c(0,1),axes=F,xlab='',ylab='')
  text(0.5,2/3,paste0(CONFIG$TUMOURNAME,' (',CONFIG$PARTICIPANTID,') - ',
                      'Run ',CONFIG$RUN,' - ',
                      ifelse(CONFIG$IS_MALE=="TRUE",'Male','Female')),cex=3)
  text(0.5,1/3,CONFIG$DIAG,cex=3)
  myPlotImage(BBprofile)
  par(mar=c(4.2,2.1,2.1,1.05))
  barplot(rev(COMBS),horiz=T,las=2,xlim=c(0,100),xlab='Genome fraction (%)',cex.lab=1.25)
  PassQC=ifelse(file.exists(paste0(CONFIG$OUTPUT_DIR,'/Refit_',CONFIG$RUN,'/PASS')),'PASS','FAIL')
  myPlotImage(DPplot,YLIM=c(0.1,0.9),TEXT=PassQC,COL=ifelse(PassQC=='PASS','green','red'))
  myPlotImage(Sunrise,XLIM=c(0.1,0.9))
  
  par(mar=c(4.2,4.2,1.05,1.05))
  if (nrow(DPinput)>0) {
  hist(DPinput$mutation.copy.number,breaks=seq(0,ceiling(max(DPinput$mutation.copy.number)),0.1),col='grey',main='',xlab='Mutation copy-number',
       cex.lab=1.25,xlim=c(0,6))
  } else {
    plot.new()
  }
  
  if (nrow(DPinput)>0) {
    plot_CNseg_mutCN(BB,DPinput)
  } else {
    plot.new()
  }

  nSVs=if (file.exists(FILES['SV'])) length(readLines(FILES['SV']))-1 else NA
  
  par(mar=rep(0.25,4))
  plot(NULL,ylim=c(0,1),xlim=c(0,1),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
  TO_PRINT=c(paste0('HD fraction (%) = ',round(HD_fraction*100,2)),
             paste0('HD size (bp) = ',formatC(HD_size,format="d",big.mark=",")),
             paste0('50% CCF size (bp) = ',formatC(Size50pCCF,format="d",big.mark=",")),
             paste0('Subclonal CNS (%) = ',round(fracSubclonalCNS*100,2)),
             paste0('Mode major/minor A = ',modeMajor,' / ',modeMinor),
             paste0('Median CN diff = ',round(medianCNDiff,4)),
             paste0('LOH fraction (%) = ',round(fracHomo*100,2)),
             paste0('WGD = ',WGD),
             paste0('GI = ',round(GI,4)),
             paste0('Coverage T/N = ',round(CovT,2),' / ',round(CovN,2)),
             paste0('MAPD = ',round(MAPD,4)),
             paste0('NRPCC = ',round(NRPCC,2)),
             paste0('SNVs = ',nrow(DPinput)),
             paste0('SVs = ',nSVs))
  text(rep(0.5,length(TO_PRINT)),seq(1,0,length.out=length(TO_PRINT)+2)[2:(length(TO_PRINT)+1)],TO_PRINT,cex=1.5)
  box()
  
  dev.off()
  
  OUT=list(Tumour_name=CONFIG$TUMOURNAME,
           Patient_ID=CONFIG$PARTICIPANTID,
           Sex=ifelse(CONFIG$IS_MALE=="TRUE",'Male','Female'),
           Diagnosis=CONFIG$DIAG,
           Run=CONFIG$RUN,
           Purity=BBmetrics$purity,
           Ploidy=BBmetrics$ploidy,
           HD_fraction=HD_fraction,
           HD_size=HD_size,
           CCF50_size=Size50pCCF,
           fracSubclonalCNS=fracSubclonalCNS,
           Mode_major=modeMajor,
           Mode_minor=modeMinor,
           Median_CN_diff=medianCNDiff,
           LOH_fraction=fracHomo,
           WGD=WGD,
           GI=GI,
           Coverage_T=CovT,
           Coverage_N=CovN,
           MAPD=MAPD,
           NRPCC=NRPCC,
           SNV=nrow(DPinput),
           SV=nSVs)
  
  return(OUT)
}

ARGS=commandArgs(trailingOnly=TRUE)
CONFIG=ARGS[1]
stopifnot(file.exists(CONFIG))
CONFIG=strsplit(readLines(CONFIG),'=')
names(CONFIG)=sapply(CONFIG,function(x) x[1])
CONFIG=lapply(CONFIG,function(x) x[2])
RUN=as.numeric(CONFIG$RUN)
stopifnot(RUN %in% 1:4)

setwd(paste0(CONFIG$OUTPUT_DIR,'/Refit_',RUN))

print("GenerateDPinputt2BamAC")
GenerateDPinputt2BamAC(
  OUTPUT_DIR=CONFIG$OUTPUT_DIR,
  TUMOURNAME=CONFIG$TUMOURNAME,
  IS_MALE=as.logical(CONFIG$IS_MALE),
  RUN=RUN)
       
print("RunDPClust")
RunDPClust(
  OUTPUT_DIR=CONFIG$OUTPUT_DIR,
  PARTICIPANTID=CONFIG$PARTICIPANTID,
  TUMOURNAME=CONFIG$TUMOURNAME,
  NORMALNAME=CONFIG$NORMALNAME,
  RUN=RUN)
       
print("AssessBBDPCRun")
AssessBBDPCRun(
  sampledir=CONFIG$OUTPUT_DIR,
  participantid=CONFIG$PARTICIPANTID,
  tumourplatekey=CONFIG$TUMOURNAME,
  normalplatekey=CONFIG$NORMALNAME,
  RUN=RUN,
  lowerclonal=0.95,
  upperclonal=1.05)
       
print("RunCNAqc")
RunCNAqc(
  sampledir=CONFIG$OUTPUT_DIR,
  TUMOURNAME=CONFIG$TUMOURNAME,
  RUN=RUN)

print('GenerateIDcard')
IDcard=GenerateIDcard(CONFIG)
save(IDcard,file=paste0(CONFIG$OUTPUT_DIR,'/Refit_',CONFIG$RUN,'/',CONFIG$TUMOURNAME,'_IDcard.Rdata'))

setwd(CONFIG$OUTPUT_DIR)
if (file.exists(paste0('Refit_',RUN,'/PASS'))) {
  system(paste0('ln -s L-FitCopyNumber_',RUN,' L-FitCopyNumber'),wait=T)
  system(paste0('ln -s M-CallSubclones_',RUN,' M-CallSubclones'),wait=T)
  system(paste0('ln -s Refit_',RUN,' Refit'),wait=T)
} else if (RUN %in% 1:2 && file.exists(paste0('Refit_',RUN,'/FAIL')) && file.exists(paste0('Refit_',RUN,'/',CONFIG$TUMOURNAME,'_metrics_run',RUN))) {
  WGD=IDcard$WGD
  DP_input_file=paste0('Refit_1/',CONFIG$TUMOURNAME,'_DPinput.txt')
  metrics_run=read.table(paste0('Refit_1/',CONFIG$TUMOURNAME,'_metrics_run1'),sep=',',stringsAsFactors=F,header=T)
  new_rho=metrics_run$purity*metrics_run$position.topccfcluster
  if (new_rho<as.numeric(CONFIG$MIN_RHO)) {
    new_rho=as.numeric(CONFIG$MIN_RHO)
  } else if (new_rho>as.numeric(CONFIG$MAX_RHO)) {
    new_rho=as.numeric(CONFIG$MAX_RHO)
  }
  new_psi=reestimate_ploidy(metrics_run$purity,
                            new_rho,
                            metrics_run$ploidy,
                            WGD)
  if (new_psi<as.numeric(CONFIG$MIN_PLOIDY)) {
    new_psi=as.numeric(CONFIG$MIN_PLOIDY)
  } else if (new_psi>as.numeric(CONFIG$MAX_PLOIDY)) {
    new_psi=as.numeric(CONFIG$MAX_PLOIDY)
  }
  dir.create(paste0("L-FitCopyNumber_",RUN+1))
  dir.create(paste0("M-CallSubclones_",RUN+1))
  dir.create(paste0("Refit_",RUN+1))
  CONFIG$RUN=RUN+1
  CONFIG$USE_PRESET_RHO_PSI='T'
  CONFIG$PRESET_RHO=new_rho
  CONFIG$PRESET_PSI=new_psi
  writeLines(paste0(names(CONFIG),'=',CONFIG),con=paste0(CONFIG$TUMOURNAME,'_configfile.txt'))
}