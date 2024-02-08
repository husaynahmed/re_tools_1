#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Battenberg))

args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
CONFIG=args[1]
stopifnot(file.exists(CONFIG))
CONFIG=strsplit(readLines(CONFIG),'=')
names(CONFIG)=sapply(CONFIG,function(x) x[1])
CONFIG=lapply(CONFIG,function(x) x[2])

setwd(CONFIG$OUTPUT_DIR)

CHR_NAME=as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (CHR_NAME==23) CHR_NAME='X'

print('Running: generate.impute.input.wgs')
generate.impute.input.wgs(chrom=CHR_NAME,
  tumour.allele.counts.file=paste0(CONFIG$OUTPUT_DIR,"/A-GetAlleleCounts/",CONFIG$TUMOURNAME,"_alleleFrequencies_chr",CHR_NAME, ".txt"),
  normal.allele.counts.file=paste0(CONFIG$OUTPUT_DIR,"/A-GetAlleleCounts/",CONFIG$NORMALNAME,"_alleleFrequencies_chr",CHR_NAME, ".txt"),
  output.file=paste0(CONFIG$OUTPUT_DIR,"/D-GenerateImputeInputFromAlleleFrequencies/",CONFIG$TUMOURNAME, "_impute_input_chr",CHR_NAME, ".txt"),
  imputeinfofile=CONFIG$IMPUTEINFOFILE,
  is.male=as.logical(CONFIG$IS_MALE),
  problemLociFile=CONFIG$PROBLEMLOCI,
  useLociFile=NA)

print('Convert Beagle5 -> Impute')
vcfbeagle=Battenberg:::convert.impute.input.to.beagle.input(
  imputeinput=paste0(CONFIG$OUTPUT_DIR,"/D-GenerateImputeInputFromAlleleFrequencies/",CONFIG$TUMOURNAME,"_impute_input_chr",CHR_NAME,".txt"),
  chrom=CHR_NAME)
vcfbeagle_path=paste0(CONFIG$OUTPUT_DIR,"/E-RunBeagle5/",CONFIG$TUMOURNAME,"_beagle5_input_chr",CHR_NAME,".vcf")
outbeagle_path=paste0(CONFIG$OUTPUT_DIR,"/E-RunBeagle5/",CONFIG$TUMOURNAME,"_beagle5_output_chr",CHR_NAME)
Battenberg:::writevcf.beagle(vcfbeagle, filepath=vcfbeagle_path,genomereference="GRCh37")
Battenberg:::run.beagle5(beaglejar=CONFIG$BEAGLE_JAR,
  vcfpath=vcfbeagle_path,
  reffile=gsub("CHROMNAME",CHR_NAME,CONFIG$BEAGLE_REF),
  outpath=outbeagle_path,
  plinkfile=gsub("CHROMNAME",CHR_NAME,CONFIG$BEAGLE_PLINK),
  maxheap.gb=15,
  nthreads=1,
  window=40,
  overlap=4,
  javajre='java')
Battenberg:::writebeagle.as.impute(
  vcf=paste0(outbeagle_path,".vcf.gz"),
  outfile=paste0(CONFIG$OUTPUT_DIR,"/E-RunBeagle5/",CONFIG$TUMOURNAME,"_impute_output_chr",CHR_NAME,"_allHaplotypeInfo.txt"))

print('Running: GetChromosomeBAFs')
GetChromosomeBAFs(chrom=CHR_NAME,
  SNP_file=paste0(CONFIG$OUTPUT_DIR,"/A-GetAlleleCounts/",CONFIG$TUMOURNAME, "_alleleFrequencies_chr", CHR_NAME, ".txt"),
  haplotypeFile=paste0(CONFIG$OUTPUT_DIR,"/E-RunBeagle5/",CONFIG$TUMOURNAME, "_impute_output_chr", CHR_NAME, "_allHaplotypeInfo.txt"),
  samplename=CONFIG$TUMOURNAME,
  outfile=paste0(CONFIG$OUTPUT_DIR,"/G-GetHaplotypedBAFs/",CONFIG$TUMOURNAME, "_chr",CHR_NAME,"_heterozygousMutBAFs_haplotyped.txt"),
  chr_names=get.chrom.names(CONFIG$IMPUTEINFOFILE, as.logical(CONFIG$IS_MALE)),
  minCounts=as.numeric(CONFIG$MIN_NORMAL_DEPTH))

print('Running: plot.haplotype.data')
plot.haplotype.data(haplotyped.baf.file=paste0(CONFIG$OUTPUT_DIR,"/G-GetHaplotypedBAFs/",CONFIG$TUMOURNAME, "_chr",CHR_NAME,"_heterozygousMutBAFs_haplotyped.txt"),
  imageFileName=paste0(CONFIG$OUTPUT_DIR,"/I-PlotHaplotypedData/",CONFIG$TUMOURNAME,"_chr",CHR_NAME,"_heterozygousData.png"),
  samplename=CONFIG$TUMOURNAME,
  chrom=CHR_NAME,
  chr_names=get.chrom.names(CONFIG$IMPUTEINFOFILE, as.logical(CONFIG$IS_MALE)))

print('All done.')
