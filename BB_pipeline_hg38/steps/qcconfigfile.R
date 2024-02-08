# -------------------------------------------------
# QC config file
# written by Alex J. Cornish & Anna Frangou
# -------------------------------------------------

# general
chrs = paste0("chr", c(1:22, "X")) # chromosomes to consider
autosomes = paste0("chr", c(1:22)) # autosomes to consider
pval.subclonal = 0.05 # P-value at which segments are called as subclonal

# thresholds
thres.purity.diff = 0.05 # samples rerun if purity difference is greater than this
thres.chrsizeincorrect.tol = 0.5 
thres.homodel.homodellargest = 1E7 
thres.propmuts = 0.01 
thres.propmuts.superclonal.or.tetra = 0.05 
thres.clonalpeak.lower = 0.95 
thres.clonalpeak.upper = 1.05 
thres.50pcpeak.lower = 0.45 
thres.50pcpeak.upper = 0.55 
thres.clonalpeak.lower.wide = 0.9 
thres.clonalpeak.upper.wide = 1.1 
thres.50pcpeak.lower.wide = 0.4 
thres.50pcpeak.upper.wide = 0.6 
thres.incorrecttetraploid.copynumber22or33 = 0.2 
thres.incorrecttetraploid.peak.diff = 0.05 
thres.incorrecttetraploid.cnodd = 0.15 
tol.around.point5 = 0.05 
thres.incorrectdiploid.aroundpoint5 = 0.15 
