
# ---------------------------------------------------------------------
# 
#   Example code - genetic r2
#   1. each code run for 10 metabs
#   2. similar code were batch run 46 times to get results for all metabolites
#   2. summary.res from 46 batch run were then combined
#   Path were masked
# 
# ---------------------------------------------------------------------

# module load R

seq = 1
seq.use = ((seq-1)*10+1):(seq*10)


library(data.table)
library(readr)
library(dplyr)
library(ieugwasr)
library(genetics.binaRies)

source("/home/TOPMed_T2D/Basefile/Functions.R")

# ---------------------------------------------------------------------------------
# load and set up 
# ---------------------------------------------------------------------------------

# load index file for all GWAS of metabolites in each cohort/meta-analysis

metab.list = read.csv("/home/TOPMed_T2D/Basefile/TOPMed_GWAS.of.Metabolites_recorded.csv")
head(metab.list)
dim(metab.list)

metab.list$uID = as.character(metab.list$uID) # uID is the harmonized metab's ID across cohorts
metab.list.use = metab.list$uID[seq.use]

metab.list$path.file = as.character(metab.list$path.file) # path.file is a var in the index file showing path and name to the GWAS summ stats
filepath.use = metab.list$path.file[seq.use]

length(metab.list.use)
length(filepath.use)
data.frame(metab.list.use=metab.list.use,filepath.use=filepath.use) # check


# ---------------------------------------------------------------------------------
# Lead SNP of metabolites
# ---------------------------------------------------------------------------------

# MAF from 1000 G

maf = load(file="/pathto_reffiles/1KG_p3v5a_hg19/chr.all.1kg.phase3.v5a.afmeu.freq.frq.RData")

# for each metabolite
# run clumping -> get top loci at p =5e-8/258 =1.091703e-10 & P<5e-8 
# if no genome-wide sig, get at 1e-7, 1e-6, then run a test for only the top 1 snp

summary.res = data.frame(uID=metab.list.use, No.1e.10=NA,r2.1e.10=NA,
	No.5e.8=NA,r2.5e.8=NA,No.1e.7=NA,r2.1e.7=NA,No.1e.6=NA,r2.1e.6=NA,P.top1=NA,r2.top1=NA)
p.list = c(1.091703e-10,5e-8,1e-7,1e-6)


for(i in 1:length(metab.list.use)) {

	# ---------------------------------
	# read in the file
	# ---------------------------------

	file = fread(filepath.use[i])
	file$SNP = as.character(file$SNP)
	file$P = as.numeric(file$P)

	file = merge(file,maf[,c("SNP","MAF")],by="SNP")

	# file for clumping

	file.p = as.data.frame(file[,c("SNP","P")])
	setnames(file.p,c("SNP","P"),c("rsid","pval"))

	PRUNED.LIST = NULL

	# ---------------------------------
	# first run the pre-set p values
	# ---------------------------------

	for(j in 1:length(p.list)) {

		p.use = p.list[j]
		
		# run clumping

		if( min(file.p$pval)<p.use ) {

			# clumping
			clump.list = ld_clump(dat = file.p,
			  		      	   	  clump_kb = 1000,
			  		      	   	  clump_r2 = 0.01,
			  		      	   	  clump_p  = p.use,
			  		      	   	  plink_bin = genetics.binaRies::get_plink_binary(),
			  		      	   	  bfile = "/pathto_reffiles/1KG_p3v5a_hg19/chr.all.1kg.phase3.v5a.afmeu")
	
			# get data to run r2
			clump.file = merge(file,data.frame(SNP=clump.list$rsid),by="SNP")
			r2 = 0

			for(k in 1:dim(clump.file)[1]) {

				beta = clump.file$Beta[k]
				p = clump.file$MAF[k]

				r2 = r2 + beta^2*2*p*(1-p)

			}

			# write results
			summary.res[i,j*2]=dim(clump.file)[1]
			summary.res[i,j*2+1]=r2


		} else {

			clump.file = NA
			r2 = 0

			# write results
			summary.res[i,j*2]=0
			summary.res[i,j*2+1]=0

		}


		# add to clump list

		pruned.list = list(clump.file)
		names(pruned.list) = paste("clump.list",j,sep='')
		PRUNED.LIST = append(PRUNED.LIST,pruned.list)

		print( paste("Metabolite ",i, ": ", metab.list.use[i], " | Completed for P < ",p.use,sep='') )

	}

	# ---------------------------------
	# column 10 and 11 - run the top SNP for 
	# those without genome-wide significant SNPs
	# ---------------------------------

	if( min(file[which(file$MAF>0.05),"P"])>5e-8 ) {

		# get the line
		topsnp = file[which(file$MAF>0.05),]
		topsnp = topsnp[which(topsnp$P==min(topsnp$P)),][1,]

		beta = topsnp$Beta
		p = topsnp$MAF
		r2 = beta^2*2*p*(1-p)

		summary.res[i,10]=topsnp$P
		summary.res[i,11]=r2

		pruned.list = list(topsnp)
		names(pruned.list) = paste("clump.list",4,sep='')


	} else {

		summary.res[i,10]=NA
		summary.res[i,11]=NA

		pruned.list = NA
		names(pruned.list) = paste("clump.list",4,sep='')


	}

	PRUNED.LIST = append(PRUNED.LIST,pruned.list)

	# rename all list

	names(PRUNED.LIST) = c("P.1.1e10","P.5e.8","P.1e.7","P.1e.6","topSNP")


	# ---------------------------------
	# write the pruned list out
	# ---------------------------------

	# clump list
	save(PRUNED.LIST,file=paste("/path_to_intermediatefile/r2.Clump.list_",metab.list.use[i],".RData",sep=''))
	
	# write out results each in case interfered
	write.csv(summary.res, paste("/path_to_results/Results_r2_",seq,".csv",sep=''), row.names=FALSE)

	# print the progress
	print( paste("Metabolite ",i, ": ", metab.list.use[i], " | Completed All",sep='') )

}


