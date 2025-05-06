#----------------------------------------------------------------------------------------------------------
#
#	Example code - Figure 4A: genetic r2
# Use Model2 for defining T2D related and T2D non-related metabolites – organized by categories
# Data-specific / path details are removed/masked 
#
#----------------------------------------------------------------------------------------------------------

library(ggplot2)

# load meta-analysis result table - "res.use" - association between each metabolite and T2D
# include results from all models, for meta-analysis of all cohorts and by race/ethnicity

head(res.use)
dim(res.use)

# a var for catagories combining super and sub pathways

res.use$general.cat = res.use$SuperP
res.use$general.cat[which(res.use$SuperP=="Lipid")] = res.use[which(res.use$SuperP=="Lipid"),"SubP.use"]

res.use$general.cat[which(res.use$general.cat=="Diacylglycerol" | 
						  res.use$general.cat=="Triacylglycerol" | 
						  res.use$general.cat=="Other Glycerolipids")] = "Glycerolipids"

res.use$general.cat[which(res.use$general.cat=="Phosphatidylcholine" | 
						  res.use$general.cat=="Phosphatidylethanolamine" | 
						  res.use$general.cat=="Other Phospholipids" | 
						  res.use$general.cat=="Lysophospholipid" |
						  res.use$general.cat=="Plasmalogen")] = "Phospholipids"

res.use$general.cat[which(res.use$general.cat=="Primary Bile Acid Metabolism" | 
						  res.use$general.cat=="Secondary Bile Acid Metabolism" |
						  res.use$general.cat=="Carnitines" |
						  res.use$general.cat=="Sphingolipid Metabolism" | 
						  res.use$general.cat=="Other Lipid Signaling")] = "Other Lipids"

res.use$general.cat[which(res.use$general.cat=="Cholesterol ester" | 
                            res.use$general.cat=="Fatty Acids" )] = "Fatty Acids/CEs"

res.use$general.cat[which(res.use$general.cat=="Energy" | 
                            res.use$general.cat=="Carbohydrate" )] = "Carbohydrate/Energy"

res.use$general.cat[which(res.use$general.cat=="Steroids" |
						  res.use$general.cat=="Cofactors and Vitamins" |
						  res.use$general.cat=="Xenobiotics")] = "Other Pathways"


# define levels 

table(res.use$general.cat,useNA="ifany")
res.use$general.cat = factor(res.use$general.cat, 
		levels = c("Amino Acid", "Carbohydrate/Energy", "Fatty Acids/CEs", "Glycerolipids", "Phospholipids", "Other Lipids","Nucleotide","Other Pathways"))

# draw with 458 metabolites with genetic data

res.use.r2 = res.use[!is.na(res.use$genetics_r2),]
res.use.r2$M2.sig = ifelse(res.use.r2$M2_All_FDR<0.05,"T2D associated","Not associated")
res.use.r2$general.cat.sig = paste(res.use.r2$general.cat,"-",res.use.r2$M2.sig)

res.use.r2$general.cat.sig = factor(res.use.r2$general.cat.sig, 
	levels = c("Amino Acid - T2D associated",         "Amino Acid - Not associated",
			       "Carbohydrate/Energy - T2D associated","Carbohydrate/Energy - Not associated",
			       "Fatty Acids/CEs - T2D associated",    "Fatty Acids/CEs - Not associated",
			       "Glycerolipids - T2D associated",      "Glycerolipids - Not associated",
			       "Phospholipids - T2D associated",      "Phospholipids - Not associated",
			       "Other Lipids - T2D associated",       "Other Lipids - Not associated",
			       "Nucleotide - T2D associated",         "Nucleotide - Not associated",
			       "Other Pathways - T2D associated",     "Other Pathways - Not associated"))


head(res.use.r2)
dim(res.use.r2)

summary(res.use.r2$genetics_r2)

# Figure 4A 

pdf("Genetic_r2_by.cat_by.assoc_updated.pdf",width=11,height=6)

ggplot(res.use.r2, aes(y = general.cat.sig, x = genetics_r2)) +
  geom_jitter(height = 0.3, aes(colour = general.cat.sig), size=2) + 
  scale_color_manual(name = "Side", values = c("#E74C3C","#F59A90",  # Amino Acid
  											   "#E67E22","#F5B47A",  # Carbohydrate/Energy
  											   "#F1C40F","#F6DC76",  # Fatty Acids/CEs
  											   "#239B56","#81EAAE",  # Glycerolipids
  											   "#0E6655","#78EBD5",  # Phospholipids
  											   "#0070C0","#7BC3F5",  # Other Lipids
  											   "#7030A0","#C894EF",  # Nucleotide
  											   "#212F3D","#5D6D7E")) + # Other Pathways
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median, 
               geom = "crossbar", 
               width = 0.8) +
  scale_y_discrete(name = "",limits=rev) + 
  scale_x_continuous(name = "") + 
  coord_cartesian(xlim = c(0, 0.68)) + 
  theme_light() + 
  theme_minimal() +

  xlab( "Genetic explaned r-square" )+
  scale_fill_continuous(guide = guide_legend(title = NULL)) + 
	 	theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
	 			panel.border = element_rect(fill=NA, colour = "black", size=0.5),
	 			axis.text = element_text(size = 16))

dev.off()


# calculate p for difference in r2 distribution, comp T2D-associated vs. non-associated metabolites

get.diff.p = data.frame(cats = as.character(unique(res.use.r2$general.cat)),p=NA)

for(i in 1:length(unique(unique(res.use.r2$general.cat)))) {

	dati = res.use.r2[which(res.use.r2$general.cat==get.diff.p[i,"cats"]),]
	get.diff.p[i,"p"] = wilcox.test(dati$genetics_r2 ~ dati$M2.sig)$p.val
	
}

print(get.diff.p)
wilcox.test(res.use.r2$genetics_r2 ~ res.use.r2$M2.sig)$p.val





