#---------------------------------------------------------------------------
#
#   Example code - coloc with GTEx tissuse-specific gene expression
#   From Buu Truong
#
#---------------------------------------------------------------------------

coloc_probs_all = NULL

# Loop over all genes
for (gene_name_i in 1:length(genelist)) {
 
    print(gene_name_i)
    message2 = ""
    
    gene_name = genelist[gene_name_i]
    print(gene_name)
    
    if (gene_name == "")
        next()
    writeLines("READING GENE FILE")

    # Get ENSG ID
    ensg_genename = gene_annot_df %>% filter(NAME == gene_name)
    ensg_genename = ensg_genename$ENSGID
    if (length(ensg_genename)>1) ensg_genename = ensg_genename[1]
    #   ensg_genename
    if (is.na(ensg_genename))
        next()
    
    # Get gene boundaries
    genes.fin = gene_annot_df %>% 
        filter(ENSGID == ensg_genename) %>%
        select(START, END, NAME, CHR)
    colnames(genes.fin) = c("Gene_start_(bp)", "Gene_end_(bp)", "gene", "CHR")
    chr = as.numeric(genes.fin$CHR)

    # Get GWAS data and lead SNPs
    gwas = sumstat %>% filter(CHR == chr)
    gwas.fin = gwas[which(gwas[,"P"]<gwas.p.threshold),]
    
    snps_gwas = gwas.fin[,"varID"]
    pos_gwas = as.numeric(gwas.fin[match(snps_gwas,gwas.fin[,"varID"]),"BP"])
    
    snps_gwas2 = NULL
    tmp = unique(snps_gwas[which(pos_gwas>=(genes.fin[1,"Gene_start_(bp)"]-gene.boundary.window.size)&pos_gwas<=(genes.fin[1,"Gene_end_(bp)"]+gene.boundary.window.size))])
    
    snps_gwas2 = c(snps_gwas2,tmp)
    lead_snps_gwas = sort(unique(snps_gwas2))
    length(lead_snps_gwas)
    
    # Get SNPs for eQTL dataset (p-value < chosen threshold)
    ensg_genename_new = unlist(strsplit(ensg_genename, split="[.]"))[1]
    eqtl.fin = eqtl[which(eqtl$eqtl_gene_id %in% c(ensg_genename, ensg_genename_new)),]
    
    if (nrow(eqtl.fin) == 0) next()

    tmp = eqtl.fin[which(eqtl.fin[,"pval_nominal"]<eqtl.p.threshold),]
    lead_snps_eqtl = eqtl.fin[which(eqtl.fin[,"pval_nominal"]<eqtl.p.threshold),"varID"]
    lead_snps_eqtl = eqtl.fin[,"varID"]

    # Get common lead SNPs between GWAS and eQTL datasets
    common_lead_snps = intersect(lead_snps_gwas,lead_snps_eqtl)
    if (length(common_lead_snps) == 0) next()
    
    lead_snps = lead_snps_gwas
    gwas.fin2 = gwas.fin[which(gwas.fin[,"varID"]%in%lead_snps),]
    gwas.fin2 = gwas.fin2[order(gwas.fin2$BP),]

    if(length(lead_snps)>1){
        rows = i = 1
        while(i<=length(lead_snps)-1) {
        i = getleadsnp.fn(gwas.fin2,baseRow=i)
        rows = c(rows,i)
        }
        snp = NA
        for(i in 1:(length(rows)-1)){
        if(rows[i+1]-rows[i]>1) {
            tmp.dat = gwas.fin2[(rows[i]:(rows[i+1]-1)),]
            tmp = tmp.dat[which(tmp.dat[,"P"]==min(tmp.dat[,"P"])),"varID"]
        }
        else {
            if (i<(length(rows)-1))
            tmp = gwas.fin2[rows[i],"varID"]; 
            if(i==(length(rows)-1)) 
            tmp = c(gwas.fin2[rows[i],"varID"],gwas.fin2[rows[i+1],"varID"])
        }
        snp=c(snp,tmp)
        }
        lead_snps = unique(snp[-1])
    }
    lead_snps
    
    if(length(lead_snps)>1){
        snp2 = NA
        rows2 = which(gwas.fin2[,"varID"]%in%lead_snps)
        gwas.fin2 = cbind(gwas.fin2,rows=c(1:nrow(gwas.fin2)))
        i = 1
        while(i<=(length(rows2)-1)){
            tmp.dat2 = gwas.fin2[c(rows2[i],(rows2[i+1])),]
            if(tmp.dat2[2,"BP"]-tmp.dat2[1,"BP"]<lead.snp.window.size) {
                index = which(tmp.dat2[,"P"]!=min(tmp.dat2[,"P"]));

                if (length(index) == 0) index = 2

                tmp2 = tmp.dat2[-index,"varID"]; 
                rows2 = rows2[-which(rows2==tmp.dat2[index,"rows"])]
            } else {
                tmp2 = c(tmp.dat2[1,"varID"],tmp.dat2[2,"varID"]); 
                i = i+1
            }
            snp2 = unique(c(snp2,tmp2))
        }
        lead_snps = unique(snp2[-1])
    }
    message("Lead SNPs -> ",paste0(lead_snps,collapse=","))

    if(gwas.response.type=="quant") coloc_probs_header = c("Region","Gene","Gene.Start","Gene.End","Trait","Tissue","GWAS.Data","CHR","GWAS.Lead.SNP","GWAS.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","ENSG_gene","Region.Num.SNPs","Region.Num.GWAS.Sig.SNPs","Best.Causal.SNP","TopSNP_eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4","Coloc2.Ratio","Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4", "Message")
    if(gwas.response.type=="cc") coloc_probs_header = c("Region","Gene","Gene.Start","Gene.End","Trait","Tissue","GWAS.Data","CHR","GWAS.Lead.SNP","GWAS.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.NCASES","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","ENSG_gene","Region.Num.SNPs","Region.Num.GWAS.Sig.SNPs","Best.Causal.SNP","TopSNP_eQTL","eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4","Coloc2.Ratio","Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4", "Message")
    
    # Loop over all lead SNPs and prepare datasets for coloc (and GCTA-COJO to identify secondary signals, if necessary)

    for(j in 1:length(lead_snps)){ 
        
        varid = lead_snps[j]
        pos = gwas.fin2[which(gwas.fin2$varID == varid),"BP"]
        
        coloc_probs = matrix(,1,length(coloc_probs_header))
        colnames(coloc_probs) = coloc_probs_header
        coloc_probs = data.frame(coloc_probs)
        for(i in 1:ncol(coloc_probs)) coloc_probs[,i] = as.character(unlist(coloc_probs[,i]))
        
        # Get boundary rows for GWAS dataset
        base = which(gwas$varID==varid)
        boundaryRows = data.frame(getWindow.fn(baseRow=base,data=gwas,size=lead.snp.window.size,jump=100))  ## Get boundary rows for GWAS dataset; these will yield all the GWAS SNPs within chosen window size on either side of lead SNP
        if(boundaryRows[[1]]<0) boundaryRows[[1]]=1
        x_gwas = gwas[boundaryRows[[1]]:boundaryRows[[2]],]
        if(length(which(x_gwas[,"FREQ"]%in%c(0,1,NA)))>0) 
        x_gwas = x_gwas[-which(x_gwas[,"FREQ"]%in%c(0,1,NA)),]
        snp = x_gwas[which(x_gwas$varID==varid),"varID"]
        
        # Count number of SNPs (significant vs total) that lie in the region surrounding lead SNP
        num.snps.sig = nrow(x_gwas[which(x_gwas[,"P"]<gwas.p.threshold),])
        num.snps = nrow(x_gwas)
        
        eqtl.fin2 = eqtl.fin[which(
        eqtl.fin[,"BP"] > gwas[boundaryRows[[1]],"BP"] &
            eqtl.fin[,"BP"]<=gwas[boundaryRows[[2]],"BP"]),]

        idx = which(eqtl.fin[,"varID"]%in%x_gwas[,"varID"])
        if (length(idx) == 0) next()
        x_eqtl = cbind(
            eqtl.fin[idx,],
            "N"=N_eqtl_tissue
        )
        x_eqtl = x_eqtl %>%
            rowwise() %>%
            mutate(maf = min(FREQ, 1-FREQ))
        x_eqtl = as.data.frame(x_eqtl)

        # Prepare GWAS and eQTL datasets for colocalization analysis
        if(gwas.response.type=="quant")  
            gwas.coloc = x_gwas[,c("CHR","varID","BP","A1","FREQ","BETA","SE","P","N")];
        if(gwas.response.type=="cc")  {
            x_gwas[,"PROPCASES"] = x_gwas[,"NCASES"]/x_gwas[,"N"]
            gwas.coloc = x_gwas[,c("CHR","varID","BP","A1","FREQ","BETA","SE","P","N","PROPCASES")]
        }
        eqtl.coloc = x_eqtl[,c("CHR","varID","BP","A1","maf","slope","slope_se","pval_nominal","N")]
        colnames(eqtl.coloc) = c("Chr","SNP","bp","refA","freq","b","se","p","n")
        if(gwas.response.type=="quant")  colnames(gwas.coloc) = colnames(eqtl.coloc)
        if(gwas.response.type=="cc") colnames(gwas.coloc) = c("Chr","SNP","bp","refA","freq","b","se","p","n","s")
        
        eqtl.coloc = eqtl.coloc[which(eqtl.coloc[,"SNP"]%in%gwas.coloc[,"SNP"]),]
        gwas.coloc = gwas.coloc[which(gwas.coloc[,"SNP"]%in%eqtl.coloc[,"SNP"]),]
        gwas.coloc = gwas.coloc[match(eqtl.coloc[,"SNP"],gwas.coloc[,"SNP"]),]
        if(length(which(is.na(gwas.coloc[,"SNP"])))>0) 
        gwas.coloc = gwas.coloc[-which(is.na(gwas.coloc[,"SNP"])),]
        eqtl.coloc = eqtl.coloc[match(gwas.coloc[,"SNP"],eqtl.coloc[,"SNP"]),]
        if(length(which(is.na(eqtl.coloc[,"SNP"])))>0) 
        eqtl.coloc = eqtl.coloc[-which(is.na(eqtl.coloc[,"SNP"])),]
        (all.equal(as.character(unlist(eqtl.coloc[,"SNP"])),as.character(unlist(gwas.coloc[,"SNP"]))))
        
        # Convert allele frequencies to minor allele frequency
        eqtl.coloc$maf = ifelse(eqtl.coloc[,"freq"]<0.5,eqtl.coloc[,"freq"],1-eqtl.coloc[,"freq"])
        gwas.coloc$maf = ifelse(gwas.coloc[,"freq"]<0.5,gwas.coloc[,"freq"],1-gwas.coloc[,"freq"])

        ## Run colocalization analysis between GWAS and eQTL datasets (identify primary signal)
        
        eqtl_ordered = eqtl.coloc[order(eqtl.coloc[,"p"],decreasing=FALSE),]
        eqtl_ordered2 = eqtl_ordered[which(eqtl_ordered[,"p"]<eqtl.p.threshold),]
        if(nrow(eqtl_ordered2)>0) 
        eqtl_topsnps = eqtl_ordered2[,"SNP"]
        region = paste0(gwas[boundaryRows[[1]],"BP"],":",gwas[boundaryRows[[2]],"BP"])
        eqtl_topsnp = eqtl.coloc[which(eqtl.coloc[,"p"]==min(eqtl.coloc[,"p"],na.rm=T)),"SNP"][1]
        
        # Prepare eQTL and GWAS datasets for colocalization analysis 
        # Run colocalization analysis between GWAS and eQTL datasets
        
        if(nrow(eqtl.coloc)>0){
        if(gwas.response.type=="quant") {

            coloc_data = data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
            colnames(coloc_data) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
            if(length(which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])))>0) coloc_data = coloc_data[-which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])),]
            
            coloc_data = coloc_data[which(!duplicated(coloc_data$gwas.SNP)),]
            
            coloc = tryCatch.W.E(coloc.abf(
            dataset1=list(
                snp=coloc_data[,"gwas.SNP"], 
                pvalues=coloc_data[,"gwas.p"], 
                N=round(median(coloc_data[,"gwas.n"]),digits=0), 
                MAF=coloc_data[,"gwas.maf"],
                type=gwas.response.type), 
            dataset2=list(
                snp=coloc_data[,"eqtl.SNP"],
                pvalues=coloc_data[,"eqtl.p"], 
                N=median(coloc_data[,"eqtl.n"]), 
                type="quant",
                MAF=coloc_data[,"eqtl.maf"]),
            p1=coloc.p1, p2=coloc.p2, p12=coloc.p12))
        }
        if(gwas.response.type=="cc") {

            coloc_data = data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se","s")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
            colnames(coloc_data) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","gwas.s","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
            if(length(which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])))>0) coloc_data = coloc_data[-which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])),]
            coloc = tryCatch.W.E(coloc.abf(dataset1=list(snp=coloc_data[,"gwas.SNP"],pvalues=coloc_data[,"gwas.p"], N=round(median(coloc_data[,"gwas.n"]),digits=0), MAF=coloc_data[,"gwas.maf"],s=coloc_data[,"gwas.s"],type=gwas.response.type), dataset2=list(snp=coloc_data[,"eqtl.SNP"],pvalues=coloc_data[,"eqtl.p"], N=median(coloc_data[,"eqtl.n"]), type="quant",MAF=coloc_data[,"eqtl.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12))

        }
            warn = inherits(coloc$warning,"warning")
            error = inherits(coloc$value,"error")
        }
        
        if (error)
            message("Error coloc 1")
        
        message = "No support (primary signal) for colocalization in this region."
        best.causal.snp = NULL
        
        if(!error) 
        if(coloc$value[[1]][[5]]<0.5 & coloc$value[[1]][[6]]>0.70)
            message = "Strong support (primary signal) for colocalization in this region!"
        
        if(!error)
            best.causal.snp = as.character(unlist(coloc$value[[2]][which(coloc$value[[2]][,"SNP.PP.H4"]==max(coloc$value[[2]][,"SNP.PP.H4"])),"snp"]))
        
        if(!error) 
            if(coloc$value[[1]][[5]]<0.5&coloc$value[[1]][[6]]>0.5&coloc$value[[1]][[6]]<=0.75)
                message = "Good support (primary signal) for colocalization in this region!"
        
        if(!error)
            if(coloc$value[[1]][[6]]<0.2) 
                message = "No support (primary signal) for colocalization in this region."
        
        if(!error)
            if(coloc$value[[1]][[5]]<0.5&coloc$value[[1]][[6]]>0.2&coloc$value[[1]][[6]]<=0.5) 
                message = "Moderate support (primary signal) for colocalization in this region."
        
        if(!error)
            if(coloc$value[[1]][[5]]>0.5) 
                message = "Evidence of LD contamination found in this region."
        
        ## Prepare output dataset
        coloc_probs[,"Message"] = message
        coloc_probs[,"GWAS.Lead.SNP"] = varid
        coloc_probs[,"Trait"] = trait
        coloc_probs[,"Tissue"] = tissue
        coloc_probs[,"GWAS.BP"] = pos
        coloc_probs[,"GWAS.P"] = x_gwas[which(x_gwas[,"varID"]==varid),"P"][1]
        coloc_probs[,"GWAS.Beta"] = x_gwas[which(x_gwas[,"varID"]==varid),"BETA"][1]
        coloc_probs[,"GWAS.SE"] = x_gwas[which(x_gwas[,"varID"]==varid),"SE"][1]
        coloc_probs[,"GWAS.A1"] = x_gwas[which(x_gwas[,"varID"]==varid),"A1"][1]
        coloc_probs[,"GWAS.A2"] = x_gwas[which(x_gwas[,"varID"]==varid),"A2"][1]
        coloc_probs[,"GWAS.Freq"] = x_gwas[which(x_gwas[,"varID"]==varid),"FREQ"][1]
        coloc_probs[,"GWAS.MAF"] = coloc_data[which(coloc_data[,"gwas.SNP"]==varid),"gwas.maf"]
        coloc_probs[,"ENSG_gene"] = ensg_genename
        coloc_probs[,"Gene"] = gene_name
        coloc_probs[,"TopSNP_eQTL"] = eqtl_topsnp
        coloc_probs[,"Best.Causal.SNP"] = best.causal.snp[1]
        coloc_probs[,"eQTL.Beta"] = x_eqtl[which(x_eqtl[,"varID"]==varid),"slope"][1]
        coloc_probs[,"CHR"] = x_eqtl[which(x_eqtl[,"varID"]==varid),"CHR"][1]
        coloc_probs[,"eQTL.SE"] = x_eqtl[which(x_eqtl[,"varID"]==varid),"slope_se"][1]
        coloc_probs[,"eQTL.P"] = x_eqtl[which(x_eqtl[,"varID"]==varid),"pval_nominal"][1]
        coloc_probs[,"eQTL.A1"] = x_eqtl[which(x_eqtl[,"varID"]==varid),"A1"][1]
        coloc_probs[,"eQTL.A2"] = x_eqtl[which(x_eqtl[,"varID"]==varid),"A2"][1]
        coloc_probs[,"eQTL.MAF"] = coloc_data[which(coloc_data[,"eqtl.SNP"]==varid),"eqtl.maf"]
        coloc_probs[,"Region"] = region
        coloc_probs[,"Region.Num.SNPs"] = num.snps
        coloc_probs[,"Region.Num.GWAS.Sig.SNPs"] = num.snps.sig
        coloc_probs[,"GWAS.N"] =  x_gwas[which(x_gwas[,"varID"]==varid),"N"][1]
        coloc_probs[,"eQTL.N"] =  x_eqtl[which(x_eqtl[,"varID"]==varid),"N"][1]
        if(gwas.response.type=="cc")   coloc_probs[,"GWAS.NCASES"] =  x_gwas[which(x_gwas[,"varID"]==varid),"NCASES"][1]
        
        if(error==FALSE){
        coloc_probs[,"Coloc.Ratio"] = signif(coloc$value[[1]]["PP.H4.abf"]/(coloc$value[[1]]["PP.H3.abf"]+coloc$value[[1]]["PP.H4.abf"]),digits=4)
        coloc_probs[,c("Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4")] = signif(coloc$value[[1]][-1],digits=4)
        coloc_probs_all = rbind(coloc_probs_all, coloc_probs)
        }
        
    }
    
    
}

if (!is.null(coloc_probs_all)) {
    coloc_probs_all$tissue = tissue_list
    write.table(coloc_probs_all, file=paste0(opt$out, "_", tissue_name, "_coloc.txt") ,row.names=F, quote=F,sep="\t")
}

writeLines("Finished all")