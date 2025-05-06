#----------------------------------------------------------------------------------------------------------
#
#  Example code - Figure 6, metabolites' assoc with risk factor vs. with T2D risk, distribution of mediation proportions
#  Data-specific / path details are removed/masked 
#
#----------------------------------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(plyr)

# load results data table - "all.res": include each metabolite's
# 1) association with T2D; 2) association and r2 with lifesytle/diet, from meta-analyses

head(all.res)
dim(all.res)


# --------------------------------------------------------------------
# Figure 6A-C and SFigures 
# Risk factor - metabolites - T2D, 2-way association comparison
# --------------------------------------------------------------------

all.res[1:2,]
dim(all.res) # 469 792

# risk factors

modifiable.factor = c("current.smk","act","BMI","redm","processedm","totalfish","poultry","totaldairy","egg","vege","nutlegume","totalfruits","potato","wholegrain","totalrefined","sugardrinks","coffeetea","alcohol")
modifiable.lables = c("Current smoking","Physical activity","BMI","Red meat","Processed meat","Total fish","Poultry","Total dairy","Egg","Vegetables","Nuts and legumes","Total fruits","Potato","Whole grains","Total refined grains","Sugary drinks","Coffee and tea","Alcohol")


# check correlations

check = data.frame(var=modifiable.factor,cor.all=NA,cor.not.sig=NA,cor.t2d.sig=NA)

for(i in 1:length(modifiable.factor)) {

  check[i,2] = cor(all.res[,"M2_All_Est"],all.res[,paste(modifiable.factor[i],"_beta",sep='')])
  check[i,3] = cor(all.res[which(all.res$M2.sig=="Not associated"),"M2_All_Est"],
      all.res[which(all.res$M2.sig=="Not associated"),paste(modifiable.factor[i],"_beta",sep='')])
  check[i,4] = cor(all.res[which(all.res$M2.sig=="T2D associated"),"M2_All_Est"],
      all.res[which(all.res$M2.sig=="T2D associated"),paste(modifiable.factor[i],"_beta",sep='')])

}

print(check)


# Amino Acid           -> #E74B3C // red
# Carbohydrate/Energy  -> #E67E22 // orange
# Fatty Acids/CEs      -> #F1C40D // yellow
# Glycerolipids        -> #239B56 // green
# Phospholipids        -> #0F6556 // dark green
# Other Lipids         -> #0070C0 // blue
# Nucleotide           -> #7130A0 // purple
# Other Pathways       -> #212F3D // dark dark blue

# T2D associated -> #595959 // dark grey
# Not associated -> #B8BBBC // light grey

all.res$general.cat = as.character(all.res$general.cat)
all.res$M2.sig = as.character(all.res$M2.sig)


pdf("Assoc_beta_risk.factor_T2D_adj.color.pdf",width=7.5,height=5)

for(i in 1:length(modifiable.factor)) {

  factori = modifiable.factor[i]
  modifiable.lablei = modifiable.lables[i]

  # limit to y-axis range but don't want to exclude dots due to it
  
  toplot = all.res[,c(paste(factori,c("_beta","_p"),sep=''),"M2_All_Est","M2.sig","general.cat")]
  toplot[which(toplot$M2_All_Est>0.5),"M2_All_Est"]=0.5
  toplot[which(toplot$M2_All_Est<(-0.5)),"M2_All_Est"]= -0.5
  
  names(toplot)[1:2] = c("beta","p")
  toplot$fdr = p.adjust(toplot$p,method='fdr',n=length(toplot$p))


  # general.cat = detailed biochemical categories - only plot it for metabs related to T2D and risk factors

  dat1 = toplot[which(toplot$M2.sig=="T2D associated" & toplot$fdr<0.05),]  # associated with T2D and food
  dat2 = toplot[which(toplot$M2.sig=="T2D associated" & toplot$fdr>=0.05),] # associated with T2D but not food
  dat3 = toplot[which(toplot$M2.sig=="Not associated"),]  # not associated with T2D 

  dat2$general.cat = "Associted with T2D but not RF"
  dat3$general.cat = "Not associted with T2D"
  
  dati = rbind(rbind(dat1,dat2),dat3)
  

  # keep all categories consistant across RF to plot consistant color, if missing, add a line with null data

  x = setdiff(c(unique(all.res$general.cat),"Associted with T2D but not RF","Not associted with T2D"),unique(dati$general.cat))
  if(length(x)>0) {
    add=data.frame(beta=NA,p=1,M2_All_Est=NA,M2.sig=NA,general.cat=x,fdr=1)
    dati = rbind(dati,add)
  }

  # Amino Acid           -> #E74B3C // red
  
  # Associted with T2D but not RF  -> #595959 // dark grey
  
  # Carbohydrate/Energy  -> #E67E22 // orange
  # Fatty Acids/CEs      -> #F1C40D // yellow
  # Glycerolipids        -> #239B56 // green

  ########## Not associted       -> #B8BBBC // light grey
  # Not associted with T2D       -> #B8BBBC // light grey
   
  # Nucleotide           -> #7130A0 // purple
  # Other Lipids         -> #0070C0 // blue
  # Other Pathways       -> #212F3D // dark dark blue
  # Phospholipids        -> #0F6556 // dark green

  ########## T2D associted       -> #595959 // dark grey
  
  # set limit for x and y axis; for coffee, set to c(-0.10,0.18) to limit blanks on the left side
  
  lim.x.max = max(abs(dat1$beta),na.rm=TRUE)
  lim.x.min = -1*lim.x.max

  lim.y.max = max(abs(dati$M2_All_Est),na.rm=TRUE)
  lim.y.min = -1*lim.y.max

  print(
    ggplot(dati, aes(x=beta, y=M2_All_Est),aes(color=as.factor(general.cat))) + 
    scale_color_manual(values=c("#E74B3C","#595959","#E67E22","#F1C40D","#239B56","#B8BBBC","#B8BBBC","#7130A0","#0070C0","#212F3D","#0F6556","#595959")) + 
    
    geom_point(shape=16,size = 2, alpha=0.8, aes(color=as.factor(general.cat))) + 
    geom_vline(xintercept = 0, color="#424242", size=0.5) + 
    geom_hline(yintercept = 0, color="#424242", size=0.5) +
    geom_smooth(aes(color=as.factor(M2.sig)), method=lm, se=FALSE, size=0.8) + 

    xlab( paste("Beta of ",modifiable.lables[i],sep='') )+
    ylab( "log(RR) of incident T2D, Model 2" )+
    xlim( c(lim.x.min,lim.x.max) )+
    ylim( c(lim.y.min,lim.y.max) )+

    scale_fill_continuous(guide = guide_legend(title = NULL)) + 

    theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill=NA, colour = "black", size=0.5),
        axis.text = element_text(size = 14))
    )

}

dev.off()
  

# --------------------------------------------------------------------
# Figure 6D-F
# BMI, PA, coffee/tea, mediation proportion by biochemical categories
# --------------------------------------------------------------------

med = read.csv("path_to_result_tables/mediation_results_meta_analysis.csv")
head(med)
dim(med)

# food group to be plotted

rfs = c("BMI","act","coffeetea")
namesinm = c("std.log_BMIcont","std.log_actcont","coffeetea")
correct.dir = c(1,-1,-1)

all.res$BMI_fdr = p.adjust(all.res$BMI_p, method="fdr", dim(all.res)[1])
all.res$act_fdr = p.adjust(all.res$act_p, method="fdr", dim(all.res)[1])
all.res$coffeetea_fdr = p.adjust(all.res$coffeetea_p, method="fdr", dim(all.res)[1])


# Amino Acid           -> #E74B3C // red
# Carbohydrate/Energy  -> #E67E22 // orange
# Fatty Acids/CEs      -> #F1C40D // yellow
# Glycerolipids        -> #239B56 // green
# Nucleotide           -> #7130A0 // purple
# Other Lipids         -> #0070C0 // blue
# Other Pathways       -> #212F3D // dark dark blue
# Phospholipids        -> #0F6556 // dark green

########## Not associted       -> #B8BBBC // light grey

colorsprf=list(c("#E74B3C","#E67E22","#F1C40D","#239B56","#7130A0","#0070C0","#212F3D","#0F6556","#B8BBBC"),
               c("#E74B3C","#F1C40D","#239B56","#0F6556","#B8BBBC"),
               c("#E74B3C","#F1C40D","#239B56","#0070C0","#212F3D","#0F6556","#B8BBBC"))
lims = c(32,32,32)

check.metab.n = data.frame(risk.factor=rfs, suppose=NA,merged=NA,step2.passDIR=NA,step2.passTEp=NA,sig=NA)


pdf("RiskFactor.metab.T2D_mediation_8.5.pdf",width=8,height=5.8)

for(i in 1:length(rfs)) {
  
  # get metabolites suppose to be in mediation
  metab.list1 = all.res[which(all.res$M2_All_FDR<0.05 & all.res[,paste(rfs[i],"_fdr",sep='')]<0.05 & 
                                sign(all.res$M2_All_Est)*sign(all.res[,paste(rfs[i],"_beta",sep='')])==correct.dir[i]),
                        c("uID","MetaboliteAnnot","Figure2.SuperP","Figure2.SubCat","general.cat")]
  
  # merge with mediation results
  dati = merge(med[which(med$Exposure == namesinm[i]),],metab.list1,by.x="Name",by.y="uID")
  
  check.metab.n[,"suppose"]=length(metab.list1)
  check.metab.n[,"merged"]=dim(dati)[1]
  
  # further filtering - align between total effect and indirect effect direction 
  dati = dati[which(sign(dati$Est_Meta_Total)*sign(dati$Est_Meta_Indirect)>0),]
  check.metab.n[i,"step2.passDIR"]=dim(dati)[1]
  
  # further filtering - total effect significant p<0.05
  dati = dati[which(sign(dati$Est_Meta_Total)==correct.dir[i] & dati$P_Meta_Total<0.05),]
  check.metab.n[i,"step2.passTEp"]=dim(dati)[1]
  
  
  # plot on qualified metabolites
  
  if (dim(dati)[1]>0) {
    
    dati$fdr_indirect = p.adjust(dati$P_Meta_Indirect,method="fdr",dim(dati)[1])
    dati$sig = ifelse(dati$fdr_indirect<0.05,"FDR<0.05","None")
    meamvalue = ddply(dati, "sig", summarise, grp.mean=mean(Rrtnde_RRtnie.1_Rrte.1))
    
    dati$mark = dati$general.cat
    dati$mark[which(dati$sig=="None")] = "z-not significant"
    
    check.metab.n[i,"sig"]=dim(dati[which(dati$sig=="FDR<0.05"),])[1]
  
    # check on which cat is not available
    x = setdiff(c(unique(all.res$general.cat),"z-not significant"),unique(dati$mark))
  
    p=ggplot(dati, aes(x=Rrtnde_RRtnie.1_Rrte.1, fill=mark, color=mark)) + 
      geom_histogram(alpha=0.9,bins=50,binwidth=0.8) + # position="dodge",
      scale_color_manual(values=colorsprf[[i]])+
      scale_fill_manual(values=colorsprf[[i]])+
      geom_vline(data=meamvalue, aes(xintercept=grp.mean),linetype="dashed",size=1)+
      xlim(0,lims[i])+
      xlab("Proportions mediated")+
      ggtitle(rfs[i])+
      theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
            panel.border = element_rect(fill=NA, colour = "black", size=0.5),
            axis.text = element_text(size = 16))
    
    print(p)
    
  }
}

dev.off()
