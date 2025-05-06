#---------------------------------------------------------------------------
#
# Example code - Figure 2
# 2A: association by biochemical categories; 2B: comparing results between race/ethnicity
# Data-specific / path details are removed/masked 
#
#---------------------------------------------------------------------------

library(ggplot2)

# codes related to loading data omitted
# load meta-analysis result table - "res.use" - association between each metabolite and T2D
# include results from all models, for meta-analysis of all cohorts and by race/ethnicity

#--------------------------------------------------------------------
# Figure 2A: Plot categories by significance of associations with T2D
#--------------------------------------------------------------------

# use model 2
res.use$color.directions = ifelse(res.use$M2_All_FDR<0.05,ifelse(res.use$M2_All_Est>=0,1,-1),0)
table(res.use$color.directions)


table(res.use$SuperP)

# get all counts with SuperP and SubP.use

data_names = res.use[,c("SuperP","SubP.use")]
data_names = data_names[!duplicated(data_names$SubP.use),]
data_names = data_names[order(data_names$SuperP,data_names$SubP.use),]

data_names$non_sig = NA
data_names$neg = NA
data_names$pos = NA


for(i in 1:dim(data_names)[1]) {

  SuperP.i = as.character(data_names[i,"SuperP"])
  SubP.use.i = as.character(data_names[i,"SubP.use"])

  nonsig = res.use[which(res.use$SuperP==SuperP.i & res.use$SubP.use==SubP.use.i & res.use$color.directions==0),]
  neg = res.use[which(res.use$SuperP==SuperP.i & res.use$SubP.use==SubP.use.i & res.use$color.directions==-1),]
  pos = res.use[which(res.use$SuperP==SuperP.i & res.use$SubP.use==SubP.use.i & res.use$color.directions==1),]

  data_names[i,"non_sig"] = dim(nonsig)[1]
  data_names[i,"neg"] = dim(neg)[1]
  data_names[i,"pos"] = dim(pos)[1]

}

print(data_names)
write.csv(data_names,"Figure2A_dated.csv")


#--------------------------------------------------------------------
# Figure 2B: Comparing results between race/ethnicity 
#--------------------------------------------------------------------

# Functional for comparison figures

figure.compare <- function(Est.names, N.names, Ncase.names, color.cat, color.pale, panel.name1, panel.name2, data, figurename) {
  
  dati = data[,c(Est.names,N.names, Ncase.names, color.cat)]
  names(dati) = c("est1", "est2", "n1", "n2", "Ncase1", "Ncase2", "color.cat")
  
  highest = max(abs(c(dati$est1, dati$est2)),na.rm=TRUE)
  lowest = -1 * highest
  
  pdf(figurename,width=6,height=4.5)
  print(ggplot(dati, aes(x=est1, y=est2, color=as.factor(color.cat))) + 
          scale_color_manual(values=color.pale) + 
          geom_point(shape=16,size = 2, aes(color=as.factor(color.cat) )) + 
          geom_vline(xintercept = 0, color="#424242", size=0.5) + 
          geom_hline(yintercept = 0, color="#424242", size=0.5) +
          xlab( paste("beta of", panel.name1) )+
          ylab( paste("beta of", panel.name2) )+
          xlim( c(lowest,highest) )+
          ylim( c(lowest,highest) )+
          ggtitle( paste("Compare",panel.name1,"vs.", panel.name2 ) ) +
          scale_fill_continuous(guide = guide_legend(title = NULL)) + 
          theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
                panel.border = element_rect(fill=NA, colour = "black", size=0.5),
                axis.text = element_text(size = 16)))
  dev.off()
  
  checklist = as.matrix(data.frame(set1=c( min(dati$n1), round(mean(dati$n1),0), median(dati$n1), max(dati$n1)),
                         set1.case = c( min(dati$Ncase1), round(mean(dati$Ncase1),0), median(dati$Ncase1), max(dati$Ncase1)),
                         set2=c( min(dati$n2), round(mean(dati$n2),0), median(dati$n2), max(dati$n2)),
                         set2.case = c( min(dati$Ncase2), round(mean(dati$Ncase2),0), median(dati$Ncase2), max(dati$Ncase2))))
  checklist = t(checklist)
  colnames(checklist) = c("min","mean","median","max")
  
  print( paste("correlation =", round(cor(dati$est1, dati$est2, use="complete.obs"), digit=4) ) )
  print( checklist )
  
}

# -----------------------------------------------------
# Comparison between non-Hispanic White vs. others
# -----------------------------------------------------

dat = res.use
dat$sig.cat = ifelse(  dat$M2_All_FDR< 0.05,"1. sig all", 
              ifelse( (dat$M2_All_FDR>=0.05 & dat$M2_wh_FDR<0.05), "2. in White", 
              ifelse( (dat$M2_All_FDR>=0.05 & dat$M2_nw_FDR<0.05), "3. in non-White", "4. not significant")))

Est.names = c("M2_wh_Est","M2_nw_Est")
N.names = c("M2_wh_N","M2_nw_N")
Ncase.names = c("M2_wh_Ncase","M2_nw_Ncase")
color.cat = "sig.cat"

color.pale = c("#EC6466",
               "#0070C0",
               "#00B050",
               "#808080")

panel.name1 = "non-Hispanic White" 
panel.name2 = "Others"

figure.compare(Est.names,N.names,Ncase.names,color.cat,color.pale, panel.name1, panel.name2, dat, "White_vs_nonWhite.pdf")


# -----------------------------------------------
# Comparison between non-Hispanic White vs. HA
# -----------------------------------------------

dat$sig.cat = ifelse(  dat$M2_All_FDR< 0.05,"1. sig all", 
              ifelse( (dat$M2_All_FDR>=0.05 & dat$M2_wh_FDR<0.05), "2. in White", 
              ifelse( (dat$M2_All_FDR>=0.05 & dat$M2_ha_FDR<0.05), "3. in HA", "4. not significant")))

Est.names = c("M2_wh_Est","M2_ha_Est")
N.names = c("M2_wh_N","M2_ha_N")
Ncase.names = c("M2_wh_Ncase","M2_ha_Ncase")
color.cat = "sig.cat"

color.pale = c("#EC6466",
               "#0070C0",
               "#FFC000",
               "#808080")

panel.name1 = "non-Hispanic White" 
panel.name2 = "US Hispanic/Latinos"

figure.compare(Est.names,N.names,Ncase.names,color.cat,color.pale, panel.name1, panel.name2, dat, "White_vs_HispanicLatino.pdf")


# -----------------------------------------------
# Comparison between non-Hispanic White vs. AA
# -----------------------------------------------

dat$sig.cat = ifelse(  dat$M2_All_FDR< 0.05,"1. sig all", 
              ifelse( (dat$M2_All_FDR>=0.05 & dat$M2_wh_FDR<0.05), "2. in White", 
              ifelse( (dat$M2_All_FDR>=0.05 & dat$M2_aa_FDR<0.05), "3. in AA", "4. not significant")))

Est.names = c("M2_wh_Est","M2_aa_Est")
N.names = c("M2_wh_N","M2_aa_N")
Ncase.names = c("M2_wh_Ncase","M2_aa_Ncase")
color.cat = "sig.cat"

color.pale = c("#EC6466",
               "#0070C0",
               "#7030A0",
               "#808080")

panel.name1 = "non-Hispanic White" 
panel.name2 = "Black"

figure.compare(Est.names,N.names,Ncase.names,color.cat,color.pale, panel.name1, panel.name2, dat, "White_vs_AA.pdf")


