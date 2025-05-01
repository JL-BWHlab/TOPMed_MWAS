
#------------------------------------------------------------------------------------------------------
#
# Metabolome wide meta-analysis with T2D across cohort
# Example code - Data processing/checking of each cohort were removed/masked
#
#------------------------------------------------------------------------------------------------------

library(metafor)
library(data.table)

# read in metabolite harmonization file

list_use = read.csv("~/Combined_Broad_Metabolome_Metabs_harmonization_date.csv",header=T)

head(list_use)
names(list_use)
dim(list_use)
 

#----------------------------------------------------------------------------------------------------------
#	Step 1 - Read in results from each cohorts (removed from example code) 
#----------------------------------------------------------------------------------------------------------

# ......
# Data processing/checking of each cohort were removed/masked
# General info in each study cohort after reading/checking

#------------
#     nhs
#------------

# results to use: nhs
# identifyer in results: uID
# identifyer in list: MetabID_NHSI
# models: 1,2,3,4,5
nhs$uID = as.character(nhs$uID)
list_use$MetabID_NHSI = as.character(list_use$MetabID_NHSI)

#------------
# nhs2
#------------

# results to use: nh2
# identifyer in results: uID
# identifyer in list: MetabID_NHSII
# models: 1,2,3,4,5
nh2$uID = as.character(nh2$uID)
list_use$MetabID_NHSII = as.character(list_use$MetabID_NHSII)

#------------
# hpfs
#------------

# results to use: hpfs
# identifyer in results: uID
# identifyer in list: MetabID_HPFS
# models: 1,2,3,4,5
hpfs$uID = as.character(hpfs$uID)
list_use$MetabID_HPFS = as.character(list_use$MetabID_HPFS)

#------------
# predimed
#------------

# results to use: predimed
# identifyer in results: uID
# identifyer in list: MetabID_PREDIMED
# models: 1,2,3,4,5
predimed$uID = as.character(predimed$uID)
list_use$MetabID_PREDIMED = as.character(list_use$MetabID_PREDIMED)

#------------
#  WHI-white
#------------

# results to use: whi_wh
# identifyer in results: uID
# identifyer in list: MetabID_WHI
# models: 1,2,3,4,5
whi_wh$uID = as.character(whi_wh$uID)
list_use$MetabID_WHI = as.character(list_use$MetabID_WHI)

#------------
# WHI-non-white
#------------

# results to use: whi_nw
# identifyer in results: uID
# identifyer in list: MetabID_WHI
# models: 1,2,3,4,5
whi_nw$uID = as.character(whi_nw$uID)
list_use$MetabID_WHI = as.character(list_use$MetabID_WHI)

#------------
#  WHI-black
#------------

# results to use: whi_aa
# identifyer in results: uID
# identifyer in list: MetabID_WHI
# models: 1,2,3,4,5
whi_aa$uID = as.character(whi_aa$uID)
list_use$MetabID_WHI = as.character(list_use$MetabID_WHI)

#------------
# SOL 
#------------

# results to use: solcox
# identifyer in results: uID
# identifyer in list: MetabID_SOL
# models: 1,2,3,4,5
solcox$uID = as.character(solcox$uID)
list_use$MetabID_SOL = as.character(list_use$MetabID_SOL)

#------------
# ARIC-black
#------------

# results to use: aric_aa
# identifyer in results: uID
# identifyer in list: MetabID_ARIC
# models: 1,2,3,4,5
aric_aa$uID = as.character(aric_aa$uID)
list_use$MetabID_ARIC = as.character(list_use$MetabID_ARIC)

#------------
# ARIC-white
#------------

# results to use: aric_wh
# identifyer in results: uID
# identifyer in list: MetabID_ARIC
# models: 1,2,3,4,5
aric_wh$uID = as.character(aric_wh$uID)
list_use$MetabID_ARIC = as.character(list_use$MetabID_ARIC)

#------------
#    FHS
#------------

# file to use: fhs
# identifyer in file: uID
# identifyer in list: MetabID_FHS
# models: 1,2,3,4,5
fhs$uID = as.character(fhs$uID)
list_use$MetabID_FHS = as.character(list_use$MetabID_FHS)

#------------
# MESA-white
#------------

# file to use: mesa_wh
# identifyer in file: uID
# identifyer in list: MetabID_MESA
# models: 1,2,3,4,5
mesa_wh$uID = as.character(mesa_wh$uID)
list_use$MetabID_MESA = as.character(list_use$MetabID_MESA)

#------------
# MESA-non-white
#------------

# file to use: mesa_nw
# identifyer in file: uID
# identifyer in list: MetabID_MESA
# models: 1,2,3,4,5
mesa_nw$uID = as.character(mesa_nw$uID)
list_use$MetabID_MESA = as.character(list_use$MetabID_MESA)

#------------
# BPRHS
#------------

# file to use: bprhs
# identifyer in file: uID
# identifyer in list: MetabID_BPRHS
# models: 1,2,3,4,5
bprhs$uID = as.character(bprhs$uID)
list_use$MetabID_BPRHS = as.character(list_use$MetabID_BPRHS)


#----------------------------------------------------------------------------------------------------------
# Step 2 -  Meta-analysis of all cohort datasets
#----------------------------------------------------------------------------------------------------------

metause <- function (estlist,semlist,plist,nlist,ncaselist) {

	dati = data.frame(est=estlist,sem=semlist,p=plist,n=nlist,ncase=ncaselist)
	dati$study=ifelse(!is.na(dati$est),1,0)
	dati$dir=ifelse(is.na(dati$est),"o",ifelse(dati$est>0,"+","-"))

	# if more than one dataset have beta - run meta 

	if( length(dati$est[!is.na(dati$est)])>1 ) {

		fit = rma(yi = dati$est[!is.na(dati$est)], sei = dati$sem[!is.na(dati$est)], method="FE",weighted=TRUE)
		metai = c(as.numeric(c(fit$beta, fit$se, fit$pval)), sum(dati$study,na.rm=TRUE),
				  sum(dati$n,na.rm=TRUE), sum(dati$ncase,na.rm=TRUE), paste(dati$dir,collapse=""), as.numeric(fit$QEp))

	# if one study have beta - use the one with data

	} else if ( length(dati$est[!is.na(dati$est)])==1 ) {

		index = which(!is.na(dati$est))

		metai = c(as.numeric(dati[index,c("est","sem","p")]), sum(dati$study,na.rm=TRUE),
				  sum(dati$n,na.rm=TRUE), sum(dati$ncase,na.rm=TRUE), paste(dati$dir,collapse=""), NA)

	# if one study have beta - use the one with data

	} else if ( length(dati$est[!is.na(dati$est)])==0 ) {

		metai = c(NA,NA,NA, sum(dati$study,na.rm=TRUE),
				  sum(dati$n,na.rm=TRUE), sum(dati$ncase,na.rm=TRUE), paste(dati$dir,collapse=""), NA)

	}

	return(metai)

}


# for each model, first extract orgainal data then meta

for(modeluse in 1:5 ) {

	coefsingle = c("Est","sem","P","N","Ncase")
	coefcombined = c("Est","sem","P","Nset","N","Ncase","Direction","HetP")
	
	Result = as.data.frame(matrix(NA,dim(list_use)[1],104))

	colnames(Result) = c("uID",

		paste("nhs_",coefsingle,sep=''), paste("nh2_",coefsingle,sep=''), paste("hpfs_",coefsingle,sep=''),
		paste("predi_",coefsingle,sep=''), paste("fhs_",coefsingle,sep=''), 
		paste("whi_wh_",coefsingle,sep=''), paste("whi_nw_",coefsingle,sep=''), paste("whi_aa_",coefsingle,sep=''),
		paste("aric_wh_",coefsingle,sep=''), paste("aric_aa_",coefsingle,sep=''),
		paste("mesa_wh_",coefsingle,sep=''), paste("mesa_nw_",coefsingle,sep=''),
		paste("sol_",coefsingle,sep=''), paste("bprhs_",coefsingle,sep=''),

		paste("All_",coefcombined,sep=''), 
		paste("wh_",coefcombined,sep=''),paste("aa_",coefcombined,sep=''),paste("ha_",coefcombined,sep=''),"HetP_ethnic"

		) # 473 104

	Result$uID = list_use$uID

	IDlist = list_use[c("uID","MetabID_NHSI","MetabID_NHSII","MetabID_HPFS","MetabID_PREDIMED","MetabID_WHI","MetabID_SOL","MetabID_ARIC",
		"MetabID_FHS","MetabID_MESA","MetabID_BPRHS","MetaboliteAnnot_B","MetaboliteAnnot_M","SuperP","SubP","Broad","Metabolome","AllN")]
	
	Result = merge(IDlist,Result,by="uID") # 473 121


	for(i in 1:dim(list_use)[1]) {

	#------------
	#     nhs
	#------------

	# results to use: nhs
	# identifyer in results: uID
	# identifyer in list: MetabID_NHSI
	# models: 1,2,3,4,5

	nhsmetabi = Result$MetabID_NHSI[i]

	if( is.na(nhsmetabi) ) {
		
		Result[i,paste("nhs_",coefsingle,sep='')] =NA

		} else if( !is.na(nhsmetabi) ) {

		Result[i,paste("nhs_",coefsingle,sep='')] = 
			nhs[which(nhs$uID==nhsmetabi), c(paste("nhs_",c("Est","sem","P"),modeluse,sep=''),"nhs_N","nhs_Ncase")]
		
	}


	#------------
	# nhs2
	#------------

	# results to use: nh2
	# identifyer in results: uID
	# identifyer in list: MetabID_NHSII
	# models: 1,2,3,4,5

	nh2metabi = Result$MetabID_NHSII[i]

	if( is.na(nh2metabi) ) {
		
		Result[i,paste("nh2_",coefsingle,sep='')] =NA

		} else if( !is.na(nh2metabi) ) {

		Result[i,paste("nh2_",coefsingle,sep='')] = 
			nh2[which(nh2$uID==nh2metabi), c(paste("nh2_",c("Est","sem","P"),modeluse,sep=''),"nh2_N","nh2_Ncase")]
		
	}


	#------------
	# hpfs
	#------------

	# results to use: hpfs
	# identifyer in results: uID
	# identifyer in list: MetabID_HPFS
	# models: 1,2,3,4,5

	hpfsmetabi = Result$MetabID_HPFS[i]

	if( is.na(hpfsmetabi) ) {
		
		Result[i,paste("hpfs_",coefsingle,sep='')] =NA

		} else if( !is.na(hpfsmetabi) ) {

		Result[i,paste("hpfs_",coefsingle,sep='')] = 
			hpfs[which(hpfs$uID==hpfsmetabi), c(paste("hpfs_",c("Est","sem","P"),modeluse,sep=''),"hpfs_N","hpfs_Ncase")]
		
	}


	#------------
	# predimed
	#------------

	# results to use: predimed
	# identifyer in results: uID
	# identifyer in list: MetabID_PREDIMED
	# models: 1,2,3,4,5

	predimetabi = Result$MetabID_PREDIMED[i]

	if( is.na(predimetabi) ) {
		
		Result[i,paste("predi_",coefsingle,sep='')] =NA

		} else if( !is.na(predimetabi) ) {

		Result[i,paste("predi_",coefsingle,sep='')] = 
			predimed[which(predimed$uID==predimetabi), c(paste("predi_",c("Est","sem","P"),modeluse,sep=''),"predi_N","predi_Ncase")]
		
	}


	#------------
	#  WHI-white
	#------------

	# results to use: whi_wh
	# identifyer in results: uID
	# identifyer in list: MetabID_WHI
	# models: 1,2,3,4,5

	whimetabi = Result$MetabID_WHI[i]

	if( is.na(whimetabi) ) {
		
		Result[i,paste("whi_wh_",coefsingle,sep='')] =NA

		} else if( !is.na(whimetabi) ) {

		Result[i,paste("whi_wh_",coefsingle,sep='')] = 
			whi_wh[which(whi_wh$uID==whimetabi), c(paste("whi_wh_",c("Est","sem","P"),modeluse,sep=''),"whi_wh_N","whi_wh_Ncase")]
		
	}


	#------------
	# WHI-non-white
	#------------

	# results to use: whi_nw
	# identifyer in results: uID
	# identifyer in list: MetabID_WHI
	# models: 1,2,3,4,5

	if( is.na(whimetabi) ) {
		
		Result[i,paste("whi_nw_",coefsingle,sep='')] =NA

		} else if( !is.na(whimetabi) ) {

		Result[i,paste("whi_nw_",coefsingle,sep='')] = 
			whi_nw[which(whi_nw$uID==whimetabi), c(paste("whi_nw_",c("Est","sem","P"),modeluse,sep=''),"whi_nw_N","whi_nw_Ncase")]
		
	}


	#------------
	# WHI-aa
	#------------

	# results to use: whi_aa
	# identifyer in results: uID
	# identifyer in list: MetabID_WHI
	# models: 1,2,3,4,5

	if( is.na(whimetabi) ) {
		
		Result[i,paste("whi_aa_",coefsingle,sep='')] =NA

		} else if( !is.na(whimetabi) ) {

		Result[i,paste("whi_aa_",coefsingle,sep='')] = 
			whi_aa[which(whi_aa$uID==whimetabi), c(paste("whi_aa_",c("Est","sem","P"),modeluse,sep=''),"whi_aa_N","whi_aa_Ncase")]
		
	}



	#------------
	# SOL 
	#------------

	# results to use: solcox
	# identifyer in results: uID
	# identifyer in list: MetabID_SOL
	# models: 1,2,3,4,5

	solmetabi = Result$MetabID_SOL[i]

	if( is.na(solmetabi) ) {
		
		Result[i,paste("sol_",coefsingle,sep='')] =NA

		} else if( !is.na(solmetabi) ) {

		Result[i,paste("sol_",coefsingle,sep='')] = 
			solcox[which(solcox$uID==solmetabi), c(paste("sol_",c("Est","sem","P"),modeluse,sep=''),"sol_N","sol_Ncase")]
		
	}


	#------------
	# ARIC-white
	#------------

	# results to use: aric_wh
	# identifyer in results: uID
	# identifyer in list: MetabID_ARIC
	# models: 1,2,3,4,5

	aricmetabi = Result$MetabID_ARIC[i]

	if( is.na(aricmetabi) ) {
		
		Result[i,paste("aric_wh_",coefsingle,sep='')] =NA

		} else if( !is.na(aricmetabi) ) {

		Result[i,paste("aric_wh_",coefsingle,sep='')] = 
			aric_wh[which(aric_wh$uID==aricmetabi), c(paste("aric_wh_",c("Est","sem","P"),modeluse,sep=''),"aric_wh_N","aric_wh_Ncase")]
		
	}

	#------------
	# ARIC-black
	#------------

	# results to use: aric_aa
	# identifyer in results: uID
	# identifyer in list: MetabID_ARIC
	# models: 1,2,3,4,5

	if( is.na(aricmetabi) ) {
		
		Result[i,paste("aric_aa_",coefsingle,sep='')] =NA

		} else if( !is.na(aricmetabi) ) {

		Result[i,paste("aric_aa_",coefsingle,sep='')] = 
			aric_aa[which(aric_aa$uID==aricmetabi), c(paste("aric_aa_",c("Est","sem","P"),modeluse,sep=''),"aric_aa_N","aric_aa_Ncase")]
		
	}


	#------------
	#    FHS
	#------------

	# file to use: fhs
	# identifyer in file: uID
	# identifyer in list: MetabID_FHS
	# models: 1,2,3,4,5

	fhsmetabi = Result$MetabID_FHS[i]

	if( is.na(fhsmetabi) ) {
		
		Result[i,paste("fhs_",coefsingle,sep='')] =NA

		} else if( !is.na(fhsmetabi) ) {

		Result[i,paste("fhs_",coefsingle,sep='')] = 
			fhs[which(fhs$uID==fhsmetabi), c(paste("fhs_",c("Est","sem","P"),modeluse,sep=''),"fhs_N","fhs_Ncase")]
		
	}


	#------------
	# MESA-white
	#------------

	# file to use: mesa_wh
	# identifyer in file: uID
	# identifyer in list: MetabID_MESA
	# models: 1,2,3,4,5

	mesametabi = Result$MetabID_MESA[i]

	if( is.na(mesametabi) ) {
		
		Result[i,paste("mesa_wh_",coefsingle,sep='')] =NA

		} else if( !is.na(mesametabi) ) {

		Result[i,paste("mesa_wh_",coefsingle,sep='')] = 
			mesa_wh[which(mesa_wh$uID==mesametabi), c(paste("mesa_wh_",c("Est","sem","P"),modeluse,sep=''),"mesa_wh_N","mesa_wh_Ncase")]
		
	}


	#------------
	# MESA-non-white
	#------------

	# file to use: mesa_nw
	# identifyer in file: uID
	# identifyer in list: MetabID_MESA
	# models: 1,2,3,4,5

	if( is.na(mesametabi) ) {
		
		Result[i,paste("mesa_nw_",coefsingle,sep='')] =NA

		} else if( !is.na(mesametabi) ) {

		Result[i,paste("mesa_nw_",coefsingle,sep='')] = 
			mesa_nw[which(mesa_nw$uID==mesametabi), c(paste("mesa_nw_",c("Est","sem","P"),modeluse,sep=''),"mesa_nw_N","mesa_nw_Ncase")]
	
	}


	#------------
	# BPRHS
	#------------

	# file to use: bprhs
	# identifyer in file: uID
	# identifyer in list: MetabID_BPRHS
	# models: 1,2,3,4,5

	bprhsmetabi = Result$MetabID_BPRHS[i]

	if( is.na(bprhsmetabi) ) {
		
		Result[i,paste("bprhs_",coefsingle,sep='')] =NA

		} else if( !is.na(bprhsmetabi) ) {

		Result[i,paste("bprhs_",coefsingle,sep='')] = 
			bprhs[which(bprhs$uID==bprhsmetabi), c(paste("bprhs_",c("Est","sem","P"),modeluse,sep=''),"bprhs_N","bprhs_Ncase")]
		
	}


	#------------
	# combine all
	#------------

	allprex = c("nhs_","nh2_","hpfs_","predi_","fhs_","whi_wh_","whi_nw_","aric_wh_","aric_aa_","mesa_wh_","mesa_nw_","sol_","bprhs_") # 13
	whprex  = c("nhs_","nh2_","hpfs_","predi_","fhs_","whi_wh_","aric_wh_","mesa_wh_") # 8
	aaprex  = c("whi_aa_","aric_aa_") # 2
	haprex  = c("mesa_nw_","sol_","bprhs_") # 3


	# apply the meta-function

	allfit = metause( estlist  = as.numeric(Result[i,paste(allprex,"Est"  ,sep='')]),
	   				  semlist  = as.numeric(Result[i,paste(allprex,"sem"  ,sep='')]),
	   				  plist    = as.numeric(Result[i,paste(allprex,"P"    ,sep='')]),
	   				  nlist    = as.numeric(Result[i,paste(allprex,"N"    ,sep='')]),
	   				  ncaselist= as.numeric(Result[i,paste(allprex,"Ncase",sep='')]))

	whfit  = metause( estlist  = as.numeric(Result[i,paste(whprex, "Est"  ,sep='')]),
	   				  semlist  = as.numeric(Result[i,paste(whprex, "sem"  ,sep='')]),
	   				  plist    = as.numeric(Result[i,paste(whprex, "P"    ,sep='')]),
	   				  nlist    = as.numeric(Result[i,paste(whprex, "N"    ,sep='')]),
	   				  ncaselist= as.numeric(Result[i,paste(whprex, "Ncase",sep='')]))

	aafit  = metause( estlist  = as.numeric(Result[i,paste(aaprex, "Est"  ,sep='')]),
	   				  semlist  = as.numeric(Result[i,paste(aaprex, "sem"  ,sep='')]),
	   				  plist    = as.numeric(Result[i,paste(aaprex, "P"    ,sep='')]),
	   				  nlist    = as.numeric(Result[i,paste(aaprex, "N"    ,sep='')]),
	   				  ncaselist= as.numeric(Result[i,paste(aaprex, "Ncase",sep='')]))

	hafit  = metause( estlist  = as.numeric(Result[i,paste(haprex, "Est"  ,sep='')]),
	   				  semlist  = as.numeric(Result[i,paste(haprex, "sem"  ,sep='')]),
	   				  plist    = as.numeric(Result[i,paste(haprex, "P"    ,sep='')]),
	   				  nlist    = as.numeric(Result[i,paste(haprex, "N"    ,sep='')]),
	   				  ncaselist= as.numeric(Result[i,paste(haprex, "Ncase",sep='')]))

	Result[i,paste("All_",c("Est","sem","P","HetP"),sep='')] = as.numeric(allfit[c(1:3,8)])
	Result[i,paste("All_",c("Nset","N","Ncase"),sep='')] = as.numeric(allfit[4:6])
	Result[i,paste("All_",c("Direction"),sep='')] = paste("d:",as.character(allfit[7])) # make sure in excel it still look good if start with "-"

	Result[i,paste("wh_", c("Est","sem","P","HetP"),sep='')] = as.numeric(whfit[c(1:3,8)])
	Result[i,paste("wh_", c("Nset","N","Ncase"),sep='')] = as.numeric(whfit[4:6])
	Result[i,paste("wh_", c("Direction"),sep='')] = paste("d:",as.character(whfit[7]))

	Result[i,paste("aa_", c("Est","sem","P","HetP"),sep='')] = as.numeric(aafit[c(1:3,8)])
	Result[i,paste("aa_", c("Nset","N","Ncase"),sep='')] = as.numeric(aafit[4:6])
	Result[i,paste("aa_", c("Direction"),sep='')] = paste("d:",as.character(aafit[7]))

	Result[i,paste("ha_", c("Est","sem","P","HetP"),sep='')] = as.numeric(hafit[c(1:3,8)])
	Result[i,paste("ha_", c("Nset","N","Ncase"),sep='')] = as.numeric(hafit[4:6])
	Result[i,paste("ha_", c("Direction"),sep='')] = paste("d:",as.character(hafit[7]))

	# test for the ethnic hetp

	racei = data.frame( est=as.numeric(Result[i,c("wh_Est","aa_Est","ha_Est")]),
						sem=as.numeric(Result[i,c("wh_sem","aa_sem","ha_sem")]))

	index = !is.na(racei$est)

	if( length(racei$est[index])>1 ) {
		Result[i,"HetP_ethnic"] = as.numeric(rma(yi = as.numeric(racei$est[index]), 
			sei = as.numeric(racei$sem[index]), method="FE",weighted=TRUE)$QEp)
	} else if( length(racei$est[index])==1 ) {
		Result[i,"HetP_ethnic"] = NA
	}


	#------------------------------------
	# print the process
	#------------------------------------

	print( paste("Model:",modeluse, " -- Metabolite:",i,sep='') )

	}

	write.csv(Result,paste("Meta.Results_Model",modeluse,"_date.csv",sep=''),row.names=FALSE)
	print( paste("Complete Model:",modeluse, " -- Total metabolite:",dim(Result)[1],sep='') )

}



#------------------------------------------------------------------------
# Step 3 -  Check and filter on the results
#------------------------------------------------------------------------

Model1 = read.csv(paste("Meta.Results_Model",1,"_date.csv",sep=''),header=T)
Model2 = read.csv(paste("Meta.Results_Model",2,"_date.csv",sep=''),header=T)
Model3 = read.csv(paste("Meta.Results_Model",3,"_date.csv",sep=''),header=T)
Model4 = read.csv(paste("Meta.Results_Model",4,"_date.csv",sep=''),header=T)
Model5 = read.csv(paste("Meta.Results_Model",5,"_date.csv",sep=''),header=T)
Model6 = read.csv(paste("Meta.Results_Model",6,"_date.csv",sep=''),header=T)
Model7 = read.csv(paste("Meta.Results_Model",7,"_date.csv",sep=''),header=T)

# calculate FDR

Model1$All_FDR = p.adjust(Model1$All_P,method='fdr',n=length(Model1$All_P))
Model2$All_FDR = p.adjust(Model2$All_P,method='fdr',n=length(Model2$All_P))
Model3$All_FDR = p.adjust(Model3$All_P,method='fdr',n=length(Model3$All_P))
Model4$All_FDR = p.adjust(Model4$All_P,method='fdr',n=length(Model4$All_P))
Model5$All_FDR = p.adjust(Model5$All_P,method='fdr',n=length(Model5$All_P))
Model6$All_FDR = p.adjust(Model6$All_P,method='fdr',n=length(Model6$All_P))
Model7$All_FDR = p.adjust(Model7$All_P,method='fdr',n=length(Model7$All_P))


# check all significant metabs

c( dim(Model1[which(Model1$All_FDR<0.05),])[1],
   dim(Model2[which(Model2$All_FDR<0.05),])[1],
   dim(Model3[which(Model3$All_FDR<0.05),])[1],
   dim(Model4[which(Model4$All_FDR<0.05),])[1],
   dim(Model5[which(Model5$All_FDR<0.05),])[1])
# 307 235 300 280 296

c( dim(Model1[which(Model1$All_FDR<0.01),])[1],
   dim(Model2[which(Model2$All_FDR<0.01),])[1],
   dim(Model3[which(Model3$All_FDR<0.01),])[1],
   dim(Model4[which(Model4$All_FDR<0.01),])[1],
   dim(Model5[which(Model5$All_FDR<0.01),])[1])
# 257 206 259 237 251


# re-write all results 

write.csv(Model1,paste("Meta.Results_Model",1,"_date.csv",sep=''),row.names=FALSE)
write.csv(Model2,paste("Meta.Results_Model",2,"_date.csv",sep=''),row.names=FALSE)
write.csv(Model3,paste("Meta.Results_Model",3,"_date.csv",sep=''),row.names=FALSE)
write.csv(Model4,paste("Meta.Results_Model",4,"_date.csv",sep=''),row.names=FALSE)
write.csv(Model5,paste("Meta.Results_Model",5,"_date.csv",sep=''),row.names=FALSE)










