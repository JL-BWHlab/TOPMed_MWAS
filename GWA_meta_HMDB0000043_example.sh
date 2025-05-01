
# --------------------------------------------------------
# 
# GWAS meta-analysis combining ARIC, SOL, NHS, CHS. FHS, WHI*3
# Analysis were run on each metabolite
# Code for HMDB0000043, shown as example
# Paths are masked
# 
# --------------------------------------------------------

# cd /scratchfolder/GWAS_Meta_Metabolites
# /homedir/PackagesTools/METAL/build/bin/metal

# IV weighted 
CLEAR
SCHEME STDERR
EFFECT_PRINT_PRECISION 12
STDERR_PRINT_PRECISION 12

# Custom filters can be used to select SNPs for inclusion in the meta-analysis.
# Exclude low MAF many of which showed extreme SE
ADDFILTER EAF > 0.01
ADDFILTER EAF < 0.99

# Eolnames
MARKER   SNP
ALLELE   A1 A2
EFFECT   Beta
STDERR   SE
PVAL     P

# ARIC
PROCESS  /pathtofolder/ARIC.GWAS.Metabolites/ARIC_EA_HMDB0000043_summary.txt.gz
PROCESS  /pathtofolder/ARIC.GWAS.Metabolites/ARIC_AA_HMDB0000043_summary.txt.gz

# SOL
PROCESS  /pathtofolder/SOL.GWAS.Metabolites/SOL_HMDB0000043_summary.txt.gz

# NHS/HPFS
PROCESS  /pathtofolder/NHS.HPFS.GWAS.Metabolites/NHS.HPFS_HMDB0000043_summary.txt.gz

# CHS 
PROCESS  /pathtofolder/CHS.GWAS.Metabolites/CHS_HMDB0000043_summary.txt.gz

# FHS
PROCESS  /pathtofolder/FHS.GWAS.Metabolites/FHS_HMDB0000043_summary.txt.gz

# WHI
PROCESS  /pathtofolder/WHI.GWAS.Metabolites/GARNET/WHI_GARNET_EA_HMDB0000043_summary.txt.gz
PROCESS  /pathtofolder/WHI.GWAS.Metabolites/SHARe/WHI_SHARe_AA_HMDB0000043_summary.txt.gz	
PROCESS  /pathtofolder/WHI.GWAS.Metabolites/WHIMS/WHI_WHIMS_EA_HMDB0000043_summary.txt.gz

# Meta-analysis
OUTFILE Meta_HMDB0000043_ .tbl
ANALYZE HETEROGENEITY



