# -------------------------------------------------
# peak calling config file
# written by Alex J. Cornish
# -------------------------------------------------
chrs = paste0("chr", 1:22) # only include autosomes to be safe
states = c("1:0", "1:1", "2:0", "2:1", "2:2") # copy number states to consider        
thres.peaks.y = 0.3 # peaks in copy-state-specific raw VAF distribution with density less than this not considered, orignally chose by looking at bi-modal log distribution, but think that was wrong
thres.peak.prop.n = 0.05 # proportion of mutations in karyotype for it to be considered


