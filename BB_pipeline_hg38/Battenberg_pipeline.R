#!/usr/bin/env Rscript

CreateScriptsDirs=function(resultsdir,participantid,tumourplatekey,tumourBAM,normalplatekey,normalBAM,
                           VCF_snv,VCF_sv,gender,Diag,
                           PIPELINE_DIR,REFERENCE_DIR,REFERENCE_FASTA) {
  submit=function(JOBNAME,TUMOURPLATEKEY,WORKDIR,WALLTIME,PIPELINE_DIR,SCRIPT,CPU=1,MEM=NULL,ARRAY=NULL,DEP=NULL) {
    LOG_PREFIX=paste0(WORKDIR,'/logs/',JOBNAME,'.%A',ifelse(!is.null(ARRAY),'.%a',''))
    query=paste0(JOBNAME,'=$(sbatch',
                 ' --job-name=',JOBNAME,'_',TUMOURPLATEKEY,
                 ifelse(!is.null(DEP),paste0(' --depend=',DEP[1],':',paste0("${",DEP[2:length(DEP)],'}',collapse=':')),''),
                 ' --time=',WALLTIME,
                 ' --ntasks=1',
                 ' --cpus-per-task=',CPU,
                 ifelse(!is.null(MEM),paste0(' --mem=',MEM),''),
                 ifelse(!is.null(ARRAY),paste0(' --array=',ARRAY),''),
                 ' -o ',LOG_PREFIX,'.out',
                 ' -e ',LOG_PREFIX,'.err',
                 ' ',PIPELINE_DIR,'/',SCRIPT,
                 ' ',WORKDIR,'/',TUMOURPLATEKEY,'_configfile.txt)')
    query=c(query,paste0('if [[ "$',JOBNAME,'" =~ Submitted\\ batch\\ job\\ ([0-9]+) ]]; then echo "${',JOBNAME,'}" && ',JOBNAME,'="${BASH_REMATCH[1]}"; elif [[ "$',JOBNAME,'" =~ ([0-9]+) ]]; then echo "${',JOBNAME,'}"; else echo "sbatch failed:" && echo "${',JOBNAME,'}" && exit 1; fi'))
    return(query)
  }
  
  stopifnot(file.exists(resultsdir))
  stopifnot(participantid!='')
  stopifnot(tumourplatekey!='')
  stopifnot(normalplatekey!='')
  stopifnot(file.exists(VCF_snv))
  if (!file.exists(VCF_sv)) warning(paste0(VCF_sv,' does not exist! Pursuing...'))
  stopifnot(file.exists(tumourBAM))
  stopifnot(file.exists(normalBAM))
  stopifnot(gender %in% c('male','female'))
  stopifnot(file.exists(PIPELINE_DIR))

  isMale=ifelse(gender=="male","TRUE","FALSE")

  # create overall sample dir
  sampledir = paste0(resultsdir,"/P_",participantid,"_T_",tumourplatekey)
  dir.create(sampledir)

  # create subdirs & subsubdirs
  dir.create(paste0(sampledir,"/A-GetAlleleCounts"))
  dir.create(paste0(sampledir,"/Additional_data"))
  dir.create(paste0(sampledir,"/B-RunBAFLogR"))
  dir.create(paste0(sampledir,"/C-RunGCcorrect"))
  dir.create(paste0(sampledir,"/D-GenerateImputeInputFromAlleleFrequencies"))
  dir.create(paste0(sampledir,"/E-RunBeagle5"))
  dir.create(paste0(sampledir,"/G-GetHaplotypedBAFs"))
  dir.create(paste0(sampledir,"/I-PlotHaplotypedData"))
  dir.create(paste0(sampledir,"/J-CombineBAFfiles"))
  dir.create(paste0(sampledir,"/K-SegmentBAFphased"))
  dir.create(paste0(sampledir,"/L-FitCopyNumber_1"))
  dir.create(paste0(sampledir,"/M-CallSubclones_1"))
  dir.create(paste0(sampledir,"/Refit_1"))
  dir.create(paste0(sampledir,"/logs"))
  
  # create basic config file (parameters will be added later if reruns are required)
  config=c()
  config=c(config,paste0("OUTPUT_DIR=",sampledir))
  config=c(config,paste0("PIPELINE_DIR=",PIPELINE_DIR))
  config=c(config,paste0("PARTICIPANTID=",participantid))
  config=c(config,paste0("TUMOURNAME=",tumourplatekey))
  config=c(config,paste0("NORMALNAME=",normalplatekey))
  config=c(config,paste0("IS_MALE=",isMale))
  config=c(config,paste0("DIAG='",Diag,"'"))
  config=c(config,paste0("TUMOURBAM=",tumourBAM))
  config=c(config,paste0("NORMALBAM=",normalBAM))
  config=c(config,"MIN_NORMAL_DEPTH=10")
  config=c(config,"PLATFORM_GAMMA=1")
  config=c(config,"SEGMENTATION_GAMMA=10")
  config=c(config,"PHASING_GAMMA=3")
  config=c(config,"CLONALITY_DIST_METRIC=0")
  config=c(config,"ASCAT_DIST_METRIC=1")
  config=c(config,"MIN_PLOIDY=1.6")
  config=c(config,"MAX_PLOIDY=4.8")
  config=c(config,"MIN_RHO=0.10") # was 0.13
  config=c(config,"MAX_RHO=1.02")
  config=c(config,"MIN_GOODNESS_OF_FIT=0.63")
  config=c(config,"BALANCED_THRESHOLD=0.51")
  config=c(config,paste0("REFERENCE_FASTA=",REFERENCE_FASTA))
  config=c(config,paste0("IMPUTEINFOFILE=",REFERENCE_DIR,"/impute_info_CAMP.txt"))
  config=c(config,paste0("SEED=",as.integer(Sys.time())))
  config=c(config,"MAX_CN_STATE=250")
  config=c(config,paste0("SV_BREAKPOINTS_FILE=",VCF_sv))
  config=c(config,paste0("VCFFILEPATH=",VCF_snv))
  config=c(config,paste0("PROBLEMLOCI=",REFERENCE_DIR,"/probloci.hg38_22072022.txt.gz"))
  config=c(config,paste0("G1000_PREFIX_AC=",REFERENCE_DIR,"/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"))
  config=c(config,paste0("G1000_PREFIX=",REFERENCE_DIR,"/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr"))
  config=c(config,paste0("GCCORRECTPREFIX=",REFERENCE_DIR,"/GC_correction_hg38/1000G_GC_chr"))
  config=c(config,paste0("RTCORRECTPREFIX=",REFERENCE_DIR,"/RT_correction_hg38/1000G_RT_chr"))
  config=c(config,"ALLELECOUNTER=/camp/apps/eb/software/alleleCount/4.0.0-foss-2016b/bin/alleleCounter")
  config=c(config,"BEAGLE_JAR=",REFERENCE_DIR,"/beagle.18May20.d20.jar")
  config=c(config,paste0("BEAGLE_REF=",REFERENCE_DIR,"/beagle5/chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz"))
  config=c(config,paste0("BEAGLE_PLINK=",REFERENCE_DIR,"/beagle5/plink.chrCHROMNAME.GRCh38.map"))
  config=c(config,"RUN=1")
  config=c(config,"USE_PRESET_RHO_PSI=F")
  config=c(config,"PRESET_RHO=NA")
  config=c(config,"PRESET_PSI=NA")
  
  writeLines(text=config,con=paste0(sampledir,"/",tumourplatekey,"_configfile.txt"))
  
  # Make config file for peaks run (make this more manipulable later, just using standard for now)
  system(paste0("cp ",PIPELINE_DIR,"/peakconfigfile.R ",sampledir,"/",tumourplatekey,"_peakconfigfile.R"))
  
  header=c()
  header=c(header,"#!/usr/bin/env bash")
  header=c(header,'')
  
  header=c(header,"# Get allele frequencies")
  header=c(header,submit('GetAlleleFrequencies',tumourplatekey,sampledir,'04:00:00',PIPELINE_DIR,
                         '1_GetAlleleCounts.sh',MEM='6GB',ARRAY='1-46'))

  header=c(header,"# Get BAF and logR + correction")
  header=c(header,submit('RunBAFLogR',tumourplatekey,sampledir,'04:00:00',PIPELINE_DIR,
                         '2_RunBAFLogR.R',DEP=c('afterok','GetAlleleFrequencies'),MEM='42GB'))
  header=c(header,'')
  
  header=c(header,"# Get additional data")
  header=c(header,submit('AdditionalData',tumourplatekey,sampledir,'04:00:00',PIPELINE_DIR,
                         '3_AdditionalData.R',DEP=c('afterok','RunBAFLogR'),MEM='12GB'))
  header=c(header,'')

  header=c(header,'# Perform phasing')
  header=c(header,submit('Phasing',tumourplatekey,sampledir,'04:00:00',PIPELINE_DIR,
                         '4_Phasing.R',DEP=c('afterok','GetAlleleFrequencies'),MEM='24GB',ARRAY='1-23'))
  header=c(header,'')
  
  header=c(header,'# Combine BAF + segmentation')
  header=c(header,submit('CombineBAFandSegment',tumourplatekey,sampledir,'04:00:00',PIPELINE_DIR,
                         '5_CombineBAFandSegment.R',DEP=c('afterok','Phasing','AdditionalData')))
  header=c(header,'')
  
  header=c(header,'# Clonal and subclonal CNA calling')
  header=c(header,submit('CNAcalling',tumourplatekey,sampledir,'24:00:00',PIPELINE_DIR,
                         '6_CNAcalling.R',MEM='12GB',DEP=c('afterok','CombineBAFandSegment')))
  header=c(header,'')
  
  header=c(header,'# QC profile')
  header=c(header,submit('QCprofile',tumourplatekey,sampledir,'12:00:00',PIPELINE_DIR,
                         '7_QCprofile.R',MEM='8GB',DEP=c('afterok','CNAcalling')))
  header=c(header,'')
  
  header=c(header,'# Clean folder')
  header=c(header,submit('Clean',tumourplatekey,sampledir,'04:00:00',PIPELINE_DIR,
                         '8_Clean.sh',DEP=c('afterok','CNAcalling')))

  writeLines(text=header,con=paste0(sampledir,'/',tumourplatekey,"_submission.sh"))
  system(paste0("chmod 744 ",sampledir,'/',tumourplatekey,"_submission.sh"))
}

########
# CORE #
########
setwd('/camp/project/proj-vanloo-secure/Hartwig/Battenberg/Battenberg_pipeline/')
METADATA=read.table('metadata_cleaned.tsv',sep='\t',header=T,stringsAsFactors=F)
OUT_dir="/camp/project/proj-vanloo-secure/Hartwig/Battenberg/Battenberg_profiles"
PIPELINE_DIR="/camp/project/proj-vanloo-secure/Hartwig/Battenberg/Battenberg_pipeline/steps"
REFERENCE_DIR="/camp/project/proj-vanloo/reference_files/human/references/Battenberg_hg38"
REFERENCE_FASTA="/camp/project/proj-vanloo/reference_files/human/references/alignment/hs38DH/hs38DH.fa" # Must be the Fasta file used when running BWA
rownames(METADATA)=METADATA$sampleId

METADATA$tumourBAM='PATH/TO/TUMOUR/CRAMorBAM'
METADATA$normalBAM='PATH/TO/NORMAL/CRAMorBAM'
METADATA$VCF_snv='PATH/TO/SOMATIC/SNVs' # Used to be paste0('/TBD/',METADATA$sampleId[i],'.purple.somatic.vcf.gz')
METADATA$VCF_sv='PATH/TO/SOMATIC/SVs'# Used to be paste0('/TBD/',METADATA$sampleId[i],'.purple.sv.vcf.gz')

for (i in 1:nrow(METADATA)) {
  print(paste0('Processing: ',METADATA$patientId[i],' (',METADATA$sampleId[i],'; ',i,'/',nrow(METADATA),')'))
  CreateScriptsDirs(resultsdir=OUT_dir,
                    participantid=METADATA$patientId[i],
                    tumourplatekey=METADATA$sampleId[i],
                    tumourBAM=METADATA$tumourBAM[i],
                    normalplatekey=paste0(METADATA$patientId[i],'_N'),
                    normalBAM=METADATA$normalBAM[i],
                    VCF_snv=METADATA$VCF_snv[i],
                    VCF_sv=METADATA$VCF_sv[i],
                    gender=METADATA$sex[i],
                    Diag=paste0(METADATA$primaryTumorLocation[i],' (',METADATA$biopsySite[i],')'),
                    PIPELINE_DIR=PIPELINE_DIR,
                    REFERENCE_DIR=REFERENCE_DIR,
                    REFERENCE_FASTA=REFERENCE_FASTA)
}; rm(i)

print('All done.')
