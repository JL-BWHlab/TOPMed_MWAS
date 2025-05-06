#----------------------------------------------------------------------------------------------------------
#
#  Example code - Figure 5, r2 of risk factors comparing T2D-associated vs. non-associated metabolites
#  Data-specific / path details are removed/masked 
#
#----------------------------------------------------------------------------------------------------------

library(ggplot2)

# load results data table - "all.res": include each metabolite's
# 1) association with T2D; 2) association and r2 with lifesytle/diet, from meta-analyses

head(all.res)
dim(all.res)

#------------------------------------
# consistant categories presentation as genetic r2

all.res$general.cat = all.res$SuperP
all.res$general.cat[which(all.res$SuperP=="Lipid")] = all.res[which(all.res$SuperP=="Lipid"),"SubP.use"]

all.res$general.cat[which(all.res$general.cat=="Diacylglycerol" | 
                          all.res$general.cat=="Triacylglycerol" | 
                          all.res$general.cat=="Other Glycerolipids")] = "Glycerolipids"

all.res$general.cat[which(all.res$general.cat=="Phosphatidylcholine" | 
                          all.res$general.cat=="Phosphatidylethanolamine" | 
                          all.res$general.cat=="Other Phospholipids" | 
                          all.res$general.cat=="Lysophospholipid" |
                          all.res$general.cat=="Plasmalogen")] = "Phospholipids"

all.res$general.cat[which(all.res$general.cat=="Primary Bile Acid Metabolism" | 
                          all.res$general.cat=="Secondary Bile Acid Metabolism" |
                          all.res$general.cat=="Carnitines" |
                          all.res$general.cat=="Sphingolipid Metabolism" | 
                          all.res$general.cat=="Other Lipid Signaling")] = "Other Lipids"

all.res$general.cat[which(all.res$general.cat=="Cholesterol ester" | 
                          all.res$general.cat=="Fatty Acids" )] = "Fatty Acids/CEs"

all.res$general.cat[which(all.res$general.cat=="Energy" | 
                          all.res$general.cat=="Carbohydrate" )] = "Carbohydrate/Energy"

all.res$general.cat[which(all.res$general.cat=="Steroids" |
                          all.res$general.cat=="Cofactors and Vitamins" |
                          all.res$general.cat=="Xenobiotics")] = "Other Pathways"

all.res$general.cat = factor(all.res$general.cat, 
                             levels = c("Amino Acid", "Carbohydrate/Energy", "Fatty Acids/CEs", "Glycerolipids", 
                                        "Phospholipids", "Other Lipids","Nucleotide","Other Pathways"))


all.res$M2.sig = ifelse(all.res$M2_All_FDR<0.05,"T2D associated","Not associated")
all.res$general.cat.sig = paste(all.res$general.cat,"-",all.res$M2.sig)

all.res$general.cat.sig = factor(all.res$general.cat.sig, 
                                 levels = c("Amino Acid - T2D associated",         "Amino Acid - Not associated",
                                            "Carbohydrate/Energy - T2D associated","Carbohydrate/Energy - Not associated",
                                            "Fatty Acids/CEs - T2D associated",    "Fatty Acids/CEs - Not associated",
                                            "Glycerolipids - T2D associated",      "Glycerolipids - Not associated",
                                            "Phospholipids - T2D associated",      "Phospholipids - Not associated",
                                            "Other Lipids - T2D associated",       "Other Lipids - Not associated",
                                            "Nucleotide - T2D associated",         "Nucleotide - Not associated",
                                            "Other Pathways - T2D associated",     "Other Pathways - Not associated"))

table(all.res$general.cat.sig)

head(all.res)
dim(all.res)


#---------------------------------------------------------------------------
# Compare averaged diet/lifesytle r2 between T2D-associated vs. non-associated metabolties
#---------------------------------------------------------------------------

foods = c("redm","processedm","totalfish","poultry","totaldairy","egg","vege","nutlegume","totalfruits",
    "potato","wholegrain","totalrefined","sugardrinks","coffeetea","alcohol","15foods")

covs = c("age","current.smk","act","BMI")

allvar = c(covs,foods)

check.cat = c("all","Amino Acid", "Carbohydrate/Energy", "Fatty Acids/CEs","Glycerolipids",
              "Phospholipids","Other Lipids","Nucleotide","Other Pathways")

behav.char = c("dir.avg_current.smk_r2","dir.avg_act_r2",
  "dir.avg_redm_r2","dir.avg_processedm_r2","dir.avg_totalfish_r2",
  "dir.avg_poultry_r2","dir.avg_totaldairy_r2","dir.avg_egg_r2",
  "dir.avg_vege_r2","dir.avg_nutlegume_r2","dir.avg_totalfruits_r2",
  "dir.avg_potato_r2","dir.avg_wholegrain_r2","dir.avg_totalrefined_r2",
  "dir.avg_sugardrinks_r2","dir.avg_coffeetea_r2","dir.avg_alcohol_r2")

behav.name = c("current.smk","act","redm","processedm","totalfish",
  "poultry","totaldairy","egg","vege","nutlegume","totalfruits",
  "potato","wholegrain","totalrefined","sugardrinks","coffeetea","alcohol")


all.res$dir.avg_total.behav_r2 = rowSums(all.res[,behav.char])


check.avg = data.frame(cat=NA,var=NA,t2d.n=NA,t2d.q25=NA,t2d.median=NA,t2d.q75=NA,t2d.iqr=NA,t2d.l=NA,t2d.h=NA,
                                     non.n=NA,non.q25=NA,non.median=NA,non.q75=NA,non.iqr=NA,non.l=NA,non.h=NA,p.nonpara=NA)
check.avg.var = c("sol_sex_r2",paste("dir.avg_",allvar,"_r2",sep=''),"dir.avg_total.behav_r2")

x=0


# for each category

for(i in 1:length(check.cat)) {

  if(check.cat[i]=="all") { 
    dati = all.res
  } else {
    dati = all.res[which(all.res$general.cat==check.cat[i]),]
  }


  # for each var

  for(j in 1:length(check.avg.var)) {

    datij = dati[,c("M2.sig",check.avg.var[j])]
    names(datij)[2] = "avg.r2"
    datij = datij[!is.na(datij$avg.r2),]

    x=x+1
    
    check.avg[x,"cat"] = check.cat[i]
    check.avg[x,"var"] = check.avg.var[j]

    t2d.i = datij[which(datij$M2.sig=="T2D associated"),]
    non.i = datij[which(datij$M2.sig=="Not associated"),]

    # ----------------
    # t2d.i

    check.avg[x,"t2d.n"] = dim(t2d.i)[1]
    check.avg[x,c("t2d.q25","t2d.median","t2d.q75")] = quantile(t2d.i$avg.r2,c(0.25,0.5,0.75))
    check.avg[x,c("t2d.iqr")] = check.avg[x,c("t2d.q75")] - check.avg[x,c("t2d.q25")]
    check.avg[x,c("t2d.l")] = max(0, check.avg[x,c("t2d.q25")] - 1.5* check.avg[x,c("t2d.iqr")])
    check.avg[x,c("t2d.h")] = check.avg[x,c("t2d.q75")] + 1.5* check.avg[x,c("t2d.iqr")]


    # ----------------
    # non.i

    check.avg[x,"non.n"] = dim(non.i)[1]
    check.avg[x,c("non.q25","non.median","non.q75")] = quantile(non.i$avg.r2,c(0.25,0.5,0.75))
    check.avg[x,c("non.iqr")] = check.avg[x,c("non.q75")] - check.avg[x,c("non.q25")]
    check.avg[x,c("non.l")] = max(0, check.avg[x,c("non.q25")] - 1.5* check.avg[x,c("non.iqr")])
    check.avg[x,c("non.h")] = check.avg[x,c("non.q75")] + 1.5* check.avg[x,c("non.iqr")]

    # ----------------
    # compare p

    if( dim(t2d.i)[1]>1 & dim(non.i)[1]>1 ) {

      check.avg[x,"p.nonpara"] = wilcox.test(datij$avg.r2~datij$M2.sig)$p.value

    } else {

      check.avg[x,"p.nonpara"] = NA

    }


  }

}

head(check.avg)
dim(check.avg)

write.csv(check.avg,"Dir.avg_lifestyles.foodgroups.comparision_dated.csv",row.names=FALSE)


#----------------------------------------------------------------------------------------------------------
# Figure 5A - for all metabolites
#----------------------------------------------------------------------------------------------------------

# name major factors, age, sex, BMI, and combined behavior

big.char = c("dir.avg_age_r2","sol_sex_r2","dir.avg_BMI_r2","dir.avg_total.behav_r2")
big.name = c("age","sex","BMI","behaviors")

# figure

pdf("All_lifestyle_foods_r2.pdf", height=11,width=8)


for(k in 1:length(check.cat)) {


	if(check.cat[k]=="all") { 
		use.dat = all.res
	} else {
		use.dat = all.res[which(all.res$general.cat==check.cat[k]),]
  	}


	# -----------------
	# major factors

	maxlim=0

	for(i in 1:length(big.char)) {

  		dati = data.frame(var = rep(big.name[i],dim(use.dat)[1]),
                    r2  = use.dat[,big.char[i]],
                    T2D.assoc = use.dat[,"M2.sig"])
		if (i==1) {
		  DAT.big = dati
		} else if(i>1) {
		  DAT.big = rbind(DAT.big,dati)
		}

  		maxlim = max(maxlim,
    	quantile(dati$r2,0.75,na.rm=TRUE)+1.5*(quantile(dati$r2,0.75,na.rm=TRUE)-quantile(dati$r2,0.25,na.rm=TRUE)))

	}

	DAT.big$var = factor(DAT.big$var, 
  		levels = c("age","sex","BMI","behaviors"),
  		labels = c("Age","Sex","BMI","Behavioral factors"))

	DAT.big$T2D.assoc = factor(DAT.big$T2D.assoc, 
  		levels = c("T2D associated","Not associated"))


	p1 = DAT.big %>%
  		 drop_na()%>%
  		 ggplot(aes(x=r2,y=var, fill=T2D.assoc)) +
  		 geom_boxplot(outlier.shape=NA, width=.6, lwd=0.6) + # 
  		 scale_fill_viridis(discrete = TRUE) + 
  		 scale_fill_manual(name = "T2D Association", values = c("#FB7C6F","#3D9DE1")) + 
  		 # geom_point(aes(colour = T2D.assoc),alpha=0.5,shape=16,size=2,
  		 #             position = position_jitterdodge()) + 
  		 scale_x_continuous(limits = c(0, 0.3)) +
  		 scale_y_discrete(limits=rev) + 

  		 ylab( "" ) +
  		 xlab( "" ) +
  		 xlim( c(0,maxlim+0.005) ) +
  		 # ggtitle( "Determination factor of r2" ) +
  		 # scale_fill_continuous(guide = guide_legend(title = NULL)) + 
    		 theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
        		 panel.border = element_rect(fill=NA, colour = "black", size=0.5),
        		 axis.text = element_text(size = 16),
        		 axis.title = element_text(size = 16))


	# -----------------
	# break down behavior factors

	maxlim=0


	for(i in 1:length(behav.char)) {

  		dati = data.frame(var = rep(behav.name[i],dim(use.dat)[1]),
                    r2  = use.dat[,behav.char[i]],
                    T2D.assoc = use.dat[,"M2.sig"])
  		if (i==1) {
  			DAT.minior = dati
  		} else if(i>1) {
  			DAT.minior = rbind(DAT.minior,dati)
  		}

  		maxlim = max(maxlim,
  		quantile(dati$r2,0.75,na.rm=TRUE)+1.5*(quantile(dati$r2,0.75,na.rm=TRUE)-quantile(dati$r2,0.25,na.rm=TRUE)))

	}

	DAT.minior$var = factor(DAT.minior$var, 
  		levels = c("current.smk","act","redm","processedm","totalfish",
  			"poultry","totaldairy","egg","vege","nutlegume","totalfruits",
  			"potato","wholegrain","totalrefined","sugardrinks","coffeetea","alcohol"),
  		labels = c("Current smoking","Physical activity","Red meat","Processed meat","Total fish",
  			"Poultry","Total dairy","Egg","vegetables","Nuts and legumes","Total fruits",
  			"Potato","Whole grain","Refined grain","Sugary drinks","Coffee and tea","Alcohol"))

	DAT.minior$T2D.assoc = factor(DAT.minior$T2D.assoc, 
  		levels = c("T2D associated","Not associated"))


	p2 = DAT.minior %>%
  		 drop_na()%>%
  		 ggplot(aes(x=r2,y=var, fill=T2D.assoc)) +
  		 geom_boxplot(outlier.shape=NA, width=.8, lwd=0.6) + # 
  		 scale_fill_viridis(discrete = TRUE) + 
  		 scale_fill_manual(name = "T2D Association", values = c("#FB7C6F","#3D9DE1")) + 
		 # geom_point(aes(colour = T2D.assoc),alpha=0.3,shape=16,size=2,
		 #              position = position_jitterdodge()) + 
  		 scale_x_continuous(limits = c(0, 0.3)) +
  		 scale_y_discrete(limits=rev) + 

  		 ylab( "" ) +
  		 xlab( "r-square") +
  		 xlim( c(0,maxlim+0.002) ) +
  		 # ggtitle( "Determination factor of r2" ) +
  		 # scale_fill_continuous(guide = guide_legend(title = NULL)) + 
  		   theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
  		       panel.border = element_rect(fill=NA, colour = "black", size=0.5),
  		       axis.text = element_text(size = 16),
  		       axis.title = element_text(size = 16))

	combine= ggarrange(p1, p2, heights = c(2.5, 6),
          ncol = 1, nrow = 2)

	print(annotate_figure(combine,top = text_grob(check.cat[k], color = "Black", face = "bold", size = 16)))

}
  
dev.off()


#----------------------------------------------------------------------------------------------------------
# Figure 5B - for significantly differnt factors - plot biochemical categories
#----------------------------------------------------------------------------------------------------------

# identify nominally significant differences / with enrichment in T2D-associated metabolites

head(check.avg)

forplot = check.avg[which(check.avg$p.nonpara<0.05 & check.avg$var!="dir.avg_15foods_r2" & 
  check.avg$var!="dir.avg_total.behav_r2" & check.avg$var!="sol_sex_r2" & check.avg$cat!="all"),]
forplot = forplot[which(forplot$t2d.median>=forplot$non.median),]

head(forplot)
dim(forplot)


forplot$cat = as.character(forplot$cat)
forplot$var = as.character(forplot$var)

all.res$general.cat = as.character(all.res$general.cat)
all.res$M2.sig = as.character(all.res$M2.sig)

table(forplot$var)

# -----------------
# figure 5B

check.max= data.frame(varname=unique(forplot$var),max=0)

for(i in 1:dim(forplot)[1]) {

 	dati = all.res[which(all.res$general.cat==forplot$cat[i]), c("general.cat",forplot$var[i],"M2.sig")]
 	dati$var = forplot$var[i]
 	names(dati) = c("Cat","r2","T2D.assoc","Var")

 	if (i==1) {
 		DAT = dati
 	} else if(i>1) {
 		DAT = rbind(DAT,dati)
 	}

 	rangemax = quantile(dati$r2,0.75,na.rm=TRUE)+1.5*(quantile(dati$r2,0.75,na.rm=TRUE)-quantile(dati$r2,0.25,na.rm=TRUE))
 	check.max[which(check.max$varname==forplot$var[i]),"max"] = max(check.max[which(check.max$varname==forplot$var[i]),"max"],rangemax)

}

head(DAT)
dim(DAT)

head(check.max)
dim(check.max)


table(DAT[,c("Var","Cat")])

DAT$Cat = fct_rev(factor(DAT$Cat,
	levels = c("Amino Acid", "Carbohydrate/Energy", "Fatty Acids/CEs", "Glycerolipids", "Phospholipids", "Other Lipids","Other Pathways")))

DAT$Var = factor(DAT$Var,
	levels = c("dir.avg_BMI_r2","dir.avg_act_r2","dir.avg_redm_r2","dir.avg_processedm_r2",
    	"dir.avg_vege_r2","dir.avg_totalfruits_r2","dir.avg_potato_r2","dir.avg_sugardrinks_r2","dir.avg_coffeetea_r2","dir.avg_alcohol_r2"),
	labels = c("BMI","Physical activity","Red meat","Processed meat","Vegetables","Total fruits","Potato",
    	"Sugary drinks","Coffee and tea","Alcohol"))

DAT$T2D.assoc = factor(DAT$T2D.assoc,
	levels = c("T2D associated","Not associated"))

# use cutoff 0.03, 0.04, and 0.08

p=  ggplot(DAT, aes(x=r2,y=Cat, fill=T2D.assoc)) +
	geom_boxplot(outlier.shape=NA, width=.8, lwd=0.6) + # 
	facet_wrap(~Var , scale="free", ncol = 1,strip.position = "left") + 

	scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
	scale_fill_manual(name = "T2D Association", values = c("#FB7C6F","#3D9DE1")) + 
	scale_x_continuous(n.breaks=4) + 
	xlim( c(0,0.08) ) +

	ylab( "" ) +
	xlab( "r-square") +

	# scale_fill_continuous(guide = guide_legend(title = NULL)) + 
	theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill=NA, colour = "black", size=0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text.y = element_blank()
          )


# convert ggplot object to grob object
gp <- ggplotGrob(p)

# optional: take a look at the grob object's layout
gtable::gtable_show_layout(gp)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.rows <- gp$layout$t[grepl("panel", gp$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
y.var <- sapply(ggplot_build(p)$layout$panel_scales_y,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp$heights[facet.rows] <- gp$heights[facet.rows] * y.var

# plot result

pdf("All_lifestyle_foods_r2_significant.cats_0.08.pdf", height=11,width=8)
grid::grid.draw(gp)
dev.off()


