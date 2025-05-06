

#---------------------------------------------------------------------------
#
#   Example code - MBE, IVW, weighted median, egger, for all T2D-associated metabs
#   1. each code run for 20 metabs
#   2. similar code were batch run 23 times to get results for all metabolites
#   2. summary.res from 23 batch run were then combined
#
#---------------------------------------------------------------------------

# module load R

seq = 1
seq.use = ((seq-1)*20+1):(seq*20)

library(data.table)
library(readr)
library(dplyr)
library(ieugwasr)
library(genetics.binaRies)
library(MendelianRandomization)
source("/home/TOPMed_T2D/Basefile/Functions.R")

# ---------------------------------------------------------------------------------
# load and set up 
# ---------------------------------------------------------------------------------

# load index file for all GWAS of metabolites in each cohort/meta-analysis

metab.list = read.csv("/n/home10/jli0525/TOPMed_T2D/Basefile/TOPMed_GWAS.of.Metabolites_recorded.csv")
head(metab.list)
dim(metab.list) # 458 31

metab.list$uID = as.character(metab.list$uID)
metab.list.use = metab.list$uID[seq.use]


# ---------------------------------------------------------------------------------
# Read in the outcomes
# ---------------------------------------------------------------------------------

X2 = fread("/pathto_traits_summstats/T2D.NatGen2022.TA_Use.txt.gz")
X2$N = X2$Ncase + X2$Ncntl

# exclude the HLA regions

X2$HLA = ifelse( (X2$Chr==6 & X2$BP>29000000 & X2$BP<34000000),1,0)
table(X2$HLA,useNA="ifany")

X2 = X2[which(X2$HLA==0),]

head(X2)
dim(X2)

# ---------------------------------------------------------------------------------
# For each metabolite - 
# Read in and merge with outcomes
# ---------------------------------------------------------------------------------

for (k in 1:length(metab.list.use)) {

	# define exposure

	exposure = metab.list.use[k]

	X1 = fread(metab.list[which(metab.list$uID==exposure),"path.file"])

# ---------------------------------------------------------------------------------
# merging with outcome
# ---------------------------------------------------------------------------------


	X1.k = X1[,c("SNP","A1","A2","Beta","SE","P")]
	names(X1.k) = c("SNP","t1_A1","t1_A2","t1_Beta","t1_Se","t1_P")
	X1.k$t1_N = metab.list[which(metab.list$uID==exposure),"all.N"]

	X2.k = X2[,c("SNP","A1","A2","Beta","Se","P","N")]
	names(X2.k) = c("SNP","t2_A1","t2_A2","t2_Beta","t2_Se","t2_P","t2_N")

	X = merge(X1.k, X2.k, by="SNP")
	X = as.data.frame(X)

	# match alleles between two files
	
	X$dir = effect_strand ("t1_A1","t1_A2","t2_A1","t2_A2", X)
	X = X[!is.na(X$dir),]
	X$t2_Beta = X$t2_Beta * X$dir

	head(X)
	dim(X)

	print( paste("Metabolite ",k, ": ", exposure, " | Completed merging",sep='') )


# ---------------------------------------
#  LD Pruning
# ---------------------------------------

	pval_thresh =c( 5e-8, 1e-7, 1e-6, 1e-5 )


	# Pruning on t1_P

	dat1 = data.frame(rsid=X[,"SNP"],pval=as.numeric(X[,"t1_P"]))

	snp.list = NA

	for(i in 1:4) {

		if( min(dat1$pval)<pval_thresh[i] ) { 

		    clump.list = ld_clump(dat = dat1,
			  		      	clump_kb = 1000,
			  		      	clump_r2 = 0.01,
			  		      	clump_p  = pval_thresh[i],
			  		      	plink_bin = genetics.binaRies::get_plink_binary(),
			  		      	bfile = "/pathto_reffiles/1KG_p3v5a_hg19/chr.all.1kg.phase3.v5a.afmeu")

		    snp.list = as.character(clump.list$rsid)

		    if(length(snp.list)>=3) {

		    	pruned1 = data.frame(dir = "metab->T2D", SNP=snp.list, P.thr=pval_thresh[i])

		    	break

		    	print( paste(exposure,i))

		   	}
		}

		
		if( i==4 & ( min(dat1$pval)>=pval_thresh[i] | length(snp.list)<3 ) ) {

			pruned1 = data.frame(dir = "metab->T2D", SNP=NA, P.thr=pval_thresh[i])

		}

	}


	# Pruning on t2_P

	dat2 = data.frame(rsid=X[,"SNP"],pval=as.numeric(X[,"t2_P"]))

	snp.list = NA

	for(i in 1:4) {

		if( min(dat2$pval)<pval_thresh[i] ) { 

		    clump.list = ld_clump(dat = dat2,
			  		      	clump_kb = 1000,
			  		      	clump_r2 = 0.01,
			  		      	clump_p  = pval_thresh[i],
			  		      	plink_bin = genetics.binaRies::get_plink_binary(),
			  		      	bfile = "/pathto_reffiles/1KG_p3v5a_hg19/chr.all.1kg.phase3.v5a.afmeu")

		    snp.list = as.character(clump.list$rsid)

		    if(length(snp.list)>=3) {

		    	pruned2 = data.frame(dir = "T2D->metab", SNP=snp.list, P.thr=pval_thresh[i])

		    	break

		    	print( paste(exposure,i))

		   	}
		}

		
		if( i==4 & ( min(dat2$pval)>=pval_thresh[i] | length(snp.list)<3 ) ) {

			pruned2 = data.frame(dir = "T2D->metab", SNP=NA, P.thr=pval_thresh[i])

		}

	}



	# Saving Pruning results

	save(pruned1,pruned2 ,file=paste("/path_to_intermediatefile/PRUNED_e_",exposure,"_o_T2D.RData",sep=''))

	print( paste("Metabolite ",k, ": ", exposure, " | Completed LD Pruning",sep='') )



# ---------------------------------------
#  Fit MBE and IVW
# ---------------------------------------

	get = data.frame(exp=c(rep(exposure,4),rep("T2D",4)),out=c(rep("T2D",4),rep(exposure,4)),
		p.thr=NA,snpno=NA,method=rep(c("IVW","MBE","Median","Egger"),2), 
		beta=NA,se=NA,p=NA,Int=NA,Int.se=NA,Int.p=NA)


	# for exposure -> outcome

		if( dim(pruned1[!is.na(pruned1$SNP),])[1]>=3 ) {

			# get parameters

			get[1:4,"p.thr"]  = min(pruned1$P.thr)
			get[1:4,"snpno"]  = dim(pruned1)[1]

		
			# get data

			MRuse <- merge(X,data.frame(SNP=as.character(pruned1$SNP),by="SNP"))

			input1 = mr_input(
			 		bx   = MRuse[,"t1_Beta"], 
		     		bxse = MRuse[,"t1_Se"], 
		     		by   = MRuse[,"t2_Beta"], 
		     		byse = MRuse[,"t2_Se"],
		     		snps = MRuse[,"SNP"], 
		     		effect_allele = MRuse[,"t1_A1"], 
		     		other_allele = MRuse[,"t1_A2"] )

			
			MRfit_ivw   = mr_ivw(input1)
			MRfit_mbe   = mr_mbe(input1)
			MRfit_wmed  = mr_median(input1, weighting="weighted")
			MRfit_egger = mr_egger(input1)


			# get the res 

			get[1,c("beta","se","p")] = c(MRfit_ivw$Estimate, MRfit_ivw$StdError, MRfit_ivw$Pvalue)
			get[2,c("beta","se","p")] = c(MRfit_mbe$Estimate, MRfit_mbe$StdError, MRfit_mbe$Pvalue)
			get[3,c("beta","se","p")] = c(MRfit_wmed$Estimate, MRfit_wmed$StdError, MRfit_wmed$Pvalue)
			get[4,c("beta","se","p")] = c(MRfit_egger$Estimate, MRfit_egger$StdError.Est, MRfit_egger$Pvalue.Est)

			get[1:4,"Int"] = MRfit_egger$Intercept
			get[1:4,"Int.se"] = MRfit_egger$StdError.Int
			get[1:4,"Int.p"] = MRfit_egger$Pvalue.Int


		} else {

			print( paste("Metabolite ",k, ": ", exposure, " | Skipped as less then 3 SNPs" ,sep='') )

		}

	
	# for exposure -> outcome

		if( dim(pruned2[!is.na(pruned2$SNP),])[1]>=3 ) {

			# get parameters

			get[5:8,"p.thr"]  = min(pruned2$P.thr)
			get[5:8,"snpno"]  = dim(pruned2)[1]

		
			# get data

			MRuse <- merge(X,data.frame(SNP=as.character(pruned2$SNP),by="SNP"))

			input2 = mr_input(
			 		bx   = MRuse[,"t2_Beta"], 
		     		bxse = MRuse[,"t2_Se"], 
		     		by   = MRuse[,"t1_Beta"], 
		     		byse = MRuse[,"t1_Se"],
		     		snps = MRuse[,"SNP"], 
		     		effect_allele = MRuse[,"t1_A1"], 
		     		other_allele = MRuse[,"t1_A2"] )

			
			MRfit_ivw   = mr_ivw(input2)
			MRfit_mbe   = mr_mbe(input2)
			MRfit_wmed  = mr_median(input2, weighting="weighted")
			MRfit_egger = mr_egger(input2)


			# get the res 

			get[5,c("beta","se","p")] = c(MRfit_ivw$Estimate, MRfit_ivw$StdError, MRfit_ivw$Pvalue)
			get[6,c("beta","se","p")] = c(MRfit_mbe$Estimate, MRfit_mbe$StdError, MRfit_mbe$Pvalue)
			get[7,c("beta","se","p")] = c(MRfit_wmed$Estimate, MRfit_wmed$StdError, MRfit_wmed$Pvalue)
			get[8,c("beta","se","p")] = c(MRfit_egger$Estimate, MRfit_egger$StdError.Est, MRfit_egger$Pvalue.Est)

			get[5:8,"Int"] = MRfit_egger$Intercept
			get[5:8,"Int.se"] = MRfit_egger$StdError.Int
			get[5:8,"Int.p"] = MRfit_egger$Pvalue.Int


		} else {

			print( paste("Metabolite ",k, ": ", exposure, " | Skipped as less then 3 SNPs" ,sep='') )

		}


	# write results
	write.csv(get, paste("/path_to_results/Results_e_",exposure,"_o_T2D.csv",sep=''), row.names=FALSE)

}






