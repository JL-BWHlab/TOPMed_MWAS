#----------------------------------------------------------------------------------------------------------
#
#  Example code - Figure 7, metabolomic signature's association with T2D, sankey plot
#  Data-specific / path details are removed/masked 
#
#----------------------------------------------------------------------------------------------------------

library(metafor)
library(data.table)
library(ggplot2)
library(psych)
library(graph)
library(network)
library(htmltools)
library(networkD3)

setwd("path_to_signature_res")

# load results data tables
# Individual cohort results harmonization & meta-analysis emoved/masked 

#------------------------------
# Figure 7D: Plot incident rate
#------------------------------

# T2D incident rate by deciles of metabolomic signature / a rbind from each individual cohort

all.t2d.rate = read.csv("Signature_T2D_rate_all.csv",header=T)
head(all.t2d.rate)
dim(all.t2d.rate)

t2d.rate.plot$decile = all.t2d.rate$decile/10-0.05 # decile 1-10 to 0.5-0.95 for x-axis
print(t2d.rate.plot)
dim(t2d.rate.plot)

pdf("Incident_T2Drate_allcohort.pdf",width=6,height=4.5)

print(ggplot(t2d.rate.plot, aes(x=decile, y=rate, color=as.factor(cohort))) + 
        scale_color_manual(values=c("#A02A93","#166082","#4EA72E","#0D9ED5","#E97133","#186B24","#a14e00","#001070")) +
        geom_point(shape=16,size = 3, aes(color=as.factor(cohort) )) + 
        geom_smooth(method="loess" , color="#424242", fill="#d3d5d5", se=TRUE, method.args = list(degree = 1), size=1) +
        # geom_vline(xintercept = 0, color="#424242", size=0.5) + 
        # geom_hline(yintercept = 0, color="#424242", size=0.5) +
        # xlab( "Percentile of metabolomic signature" )+
        ylab( "Incident T2D rate" )+
        xlim( c(0,1) )+
        ylim( c(0,0.6) )+
        ggtitle( "Incident rate my metabolimic signature" ) +
        scale_fill_continuous(guide = guide_legend(title = NULL)) + 
        theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), 
              panel.border = element_rect(fill=NA, colour = "black", size=1),
              axis.text = element_text(size = 16)))

dev.off()


#----------------------------------------------------------------------------------------------------------
# Figure 7F: Sankey plot
#----------------------------------------------------------------------------------------------------------

# Load the signature model

model = data.frame(read_excel("Summary_MetabSignature_AUC_20250303.xlsx", sheet = "Models"))[,c(1,3)]
names(model)[2] = "FinalModel"
model$FinalModel = as.numeric(model$FinalModel)
head(model)
dim(model)

# risk factor-metabolite association results from meta-analysis

res = read.csv("path_to_T2D_risk_factor_analysis_res/All.Assoc.T2D.Foodgroups.lifestyle.csv")
res = res[,c("uID","MetaboliteAnnot","Lipid.Type","Lipid.CL","Lipid.DB","SuperP","SubP.use",
             "M2_All_Est","M2_All_sem","M2_All_P","M2_All_FDR","M2_All_Nset","M2_All_N","M2_All_Ncase","M2_All_Direction",
             names(res)[712:812])]
head(res)
dim(res)

# merge data

signat_1 = merge(res,model,by="uID",all.x=T)
head(signat_1)
dim(signat_1)

signat_1use = signat_1[!is.na(signat_1$FinalModel),]
head(signat_1use)
dim(signat_1use)


# load the metabolomic signature' association with risk factors and T2D

sig1 = read.csv("Risk.factor_Signature_MetaAssoc.csv",header=T)
head(sig1)
dim(sig1)

sig2 = read.csv("Signature_T2D_MetaAssoc.csv",header=T)
head(sig2)
dim(sig2)


#-------------------------------------------
# Sankey plot for the following risk factors
#-------------------------------------------

sig1[which(sig1$meta_FDR<0.05),"RiskFactor"]

# act
# alcohol
# BMI
# coffeetea
# redm
# sugardrinks

# source=risk factors | target=MetS+metabolite | link=association
# source=MetS+metabolite | target=increase/decrease T2D | link=association strengths

use.modifiable.factor = c("BMI","act","coffeetea","redm","sugardrinks","alcohol")
use.modifiable.lables = c("BMI","Physical activity","Coffee and tea","Red meat","Sugary drinks","Wine")

# association related to metabolomic signature

set0 = sig1[which(sig1$meta_FDR<0.05),c("RiskFactor","meta_Beta")]
set0 = data.frame(source=set0$RiskFactor,target="Metabolomic Signature",value=abs(set0$meta_Beta),assoc_type=sign(set0$meta_Beta))
rownames(set0) = set0$source
set0 = set0[use.modifiable.factor,]
set0$source = use.modifiable.lables
print(set0)

# association related to metabolites

set1 = NULL

for(i in 1:length(use.modifiable.factor)) {
  
  pvar = paste(use.modifiable.factor[i],"_fdr",sep='')
  bvar = paste(use.modifiable.factor[i],"_beta",sep='')
  
  dati = signat_1use[which(signat_1use[,pvar]<0.05),]
  dati = dati[order(dati[,bvar]),]
  
  if ( dim(dati)[1]>0 ) {
    
    seti = data.frame(source=use.modifiable.lables[i],target=dati$MetaboliteAnnot,value=abs(dati[,bvar]),
                      assoc_type=sign(dati[,bvar]))
    set1 = rbind(set1,seti)
    
  }
  
}

head(set1)
dim(set1)

# merge the left part 

set1 = rbind(set0,set1)
head(set1,n=12)
dim(set1)


# right part - by M2 association

seti.0 = data.frame(source="Metabolomic Signature",target="Higher T2D risk",value=abs(sig2[which(sig2$Var=="INT_persd"),"meta_Est2"]),assoc_type=1)

seti.a = signat_1use[which(signat_1use$M2_All_Est>0),]
seti.b = signat_1use[which(signat_1use$M2_All_Est<0),]

set2.a = data.frame(source=seti.a$MetaboliteAnnot,target="Higher T2D risk",value=abs(seti.a$M2_All_Est),assoc_type= 1)
set2.b = data.frame(source=seti.b$MetaboliteAnnot,target="Lower T2D risk",value=abs(seti.b$M2_All_Est),assoc_type= -1)

set2 = rbind(seti.0,rbind(set2.a,set2.b))

# set up the Link files

Link = rbind(set1,set2)

# check Nodes names

Nodes = NULL
for(i in 1:dim(Link)[1]) {
  Nodes = c(Nodes,as.character(Link[i,c("source","target")]))
}
Nodes = unique(Nodes)
Nodes = data.frame(order=1:length(Nodes)-1,name=Nodes)


for(i in 1:dim(Link)[1]) {
  
  si = Link[i,"source"]
  Link[i,"source"] = Nodes[which(Nodes$name==si),"order"]
  
  ti = Link[i,"target"]
  Link[i,"target"] = Nodes[which(Nodes$name==ti),"order"]
  
}

# link

Link$source = as.numeric(Link$source)
Link$target = as.numeric(Link$target)
Link$value = round(Link$value*100,digit=2)
Link$assoc_type = ifelse(Link$assoc_type==1,"postive","negtive")

# node group

Nodes$group = "Metabolites"

for(i in 1:dim(Nodes)[1]) {
  
  if ( length(intersect(Nodes$name[i],c("BMI","Red meat","Sugary drinks")))>0) {
      Nodes$group[i] = "bad lifestyle"
  
  } else if ( length(intersect(Nodes$name[i],c("Physical activity","Coffee and tea","Wine")))>0) {
    Nodes$group[i] = "good lifestyle"
    
  } else if (Nodes$name[i] =="Higher T2D risk" | Nodes$name[i] =="Lower T2D risk" | Nodes$name[i] =="Metabolomic Signature") {
      Nodes$group[i] = Nodes$name[i]
  
  }
}

Nodes=data.frame(Nodes[,c("name","group")])

# check

length(unique(c(Link$source,Link$target)))
dim(Nodes)

# plot

SigSankey = list(Nodes=Nodes,Links=Link)


# "bad lifestyle","Metabolomic Signature","good lifestyle","Higher T2D risk","Lower T2D risk","Metabolites","postive","negtive"
# "#8c53ff","#FF9300","#00B050","#FF5552","#4D8ECD","#B5B5B5","#FFC9C8","#B7D0EA"

my_color <- 'd3.scaleOrdinal() .domain(["bad lifestyle","Metabolomic Signature","good lifestyle","Higher T2D risk","Lower T2D risk","Metabolites","postive","negtive"]) .range(["#8c53ff","#FF9300","#00B050","#FF5552","#4D8ECD","#B5B5B5","#FFC9C8","#B7D0EA"])'

sankeyNetwork(Links = SigSankey$Links, Nodes = SigSankey$Nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", LinkGroup = "assoc_type",NodeGroup="group",colourScale=my_color,
              fontSize = 16, fontFamily="Arial",nodeWidth = 30, height = 1000, width = 650) %>% saveNetwork(file = "Sankey_selected_riskfactors.html")




