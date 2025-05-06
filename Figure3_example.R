#----------------------------------------------------------------------------------------------------------
#
#	 Example code - Figure 3, display significant metabolites using model 2 results
#  Data-specific / path details are removed/masked 
#
#----------------------------------------------------------------------------------------------------------

library(ggplot2)

# load meta-analysis result table - "res.use" - association between each metabolite and T2D
# include results from all models, for meta-analysis of all cohorts and by race/ethnicity

head(res.use)
dim(res.use)

	
	#--------------------------------------------------------------
	# for sub-figure A
	#--------------------------------------------------------------
	# note that figure 3 was initially figure 2 thus variable names were figure2

	dat = res.use[which(res.use$Figure2==2),]
	head(dat)
	dim(dat) # 118 272
	
	
	# 1. p values as height, color as association magnitudes and directions // winsorize the results at a limit to make sure meanningful results demonstration
	
	dat$height = -log10(dat$M2_All_FDR)
	dat$height[which(dat$height>20)]=20 # if more sig than 1e-20, draw at 1e-20
	
	dat$colors = dat$M2_All_Est
	dat$colors[which(abs(dat$colors)>0.3)] = sign(dat$colors[which(abs(dat$colors)>0.3)])*0.3 # if stronger than 0.3, draw at 0.3
	
	
	# 2. seperate super paths with a missing colum --> add 1 line after each Figure2.SuperP
	#    clear way for y-axis --> add 3 lines after the final Figure2.SuperP
	
	to_add = data.frame(matrix(NA,1*length(unique(dat$Figure2.SuperP)),ncol(dat)))
	colnames(to_add) = colnames(dat)
	to_add$Figure2.SuperP = unique(dat$Figure2.SuperP)
	to_add$Figure2.SubCat = "zzz_grid"
	
	to_add = to_add[order(to_add$Figure2.SuperP),]
	to_add = rbind(to_add,to_add[dim(to_add)[1],])
	to_add = rbind(to_add,to_add[dim(to_add)[1],])
	
	dat = rbind(dat,to_add)
	
	
	# 3. order data by effect size (M2_All_Est), defined seq (Figure2.seq), and super/sub pathway (Figure2.SubCat & Figure2.SuperP)
	
	dat = dat[order(dat$M2_All_Est),]
	dat = dat[order(dat$Figure2.seq),]
	dat = dat[order(dat$Figure2.SubCat),]
	dat = dat[order(dat$Figure2.SuperP),]
	
	dat$seq = 1:dim(dat)[1]
	

	# 4. angle of the direction of metab annot 
	
	angle = 90 - 360 * (dat$seq-0.5) /dim(dat)[1]
	dat$hjust = ifelse( angle < -90, 1, 0)
	dat$angle = ifelse(angle < -90, angle+180, angle)
	
		
	# 5. prepare a data frame for base lines
	
	dat$Figure2.SuperP = as.character(dat$Figure2.SuperP)
	cats = unique(dat$Figure2.SuperP)
	cats = cats[order(cats)]
	base_data = data.frame(group=cats,start=NA,end=NA,title=NA)
	
	for(i in 1:length(cats)) {
	  
	  substracti = dat[which(dat[,"Figure2.SuperP"]==cats[i]),]
	  base_data[i,"start"]=min(substracti$seq)
	  
	  if(i<length(cats)) {
	    base_data[i,"end"]=max(substracti$seq)-1 
	    
	  } else if (i==length(cats)) {
	    base_data[i,"end"]=max(substracti$seq)-3
	    
	  }
	  
	  base_data[i,"title"] = mean( c(base_data[i,"start"],base_data[i,"end"]) )
	  
	}
	
	base_data = base_data[order(base_data$start),]
	
	
	# 6. prepare a data frame for grid (scales)
	
	grid_data = data.frame(group=c("first","second"),start=NA,end=NA)
	grid_data[1,c("start","end")] = c(base_data[dim(base_data)[1],"end"]+0.5,base_data[dim(base_data)[1],"end"]+1.5)
	grid_data[2,c("start","end")] = c(max(dat$seq)-0.5,max(dat$seq)+0.5)
	
	
# ---------
# 7. plot figure

  p2= ggplot(dat) +

	# add the bars

    geom_bar( aes(x=seq, y=height, fill=colors), stat="identity") + 
    coord_polar() +

    # add color pattele

    scale_fill_gradient2(low = "#006FC0", mid = "white", high = "#FF5552", midpoint = 0, na.value = "grey90", n.breaks=10) + 
    ylim(c(-8,23)) + 
    theme_light() + 
    theme_minimal() +

    theme(
		legend.position = "none",
		axis.text = element_blank(),
		axis.title = element_blank(),
		panel.grid = element_blank(),
		plot.margin = unit(rep(-1,4), "cm") 
  	) +


    # add value = 5,10,15 lines 

    geom_segment(data=grid_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = start, y = 5, xend = end, yend = 5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = start, y = 10, xend = end, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = start, y = 15, xend = end, yend = 15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = start, y = 20, xend = end, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    

    # Add text showing the value of each 100/75/50/25 lines

    ggplot2::annotate("text", x = rep(max(dat$seq)-0.5,5), y = c(0, 5, 10, 15,20), label = c("0", "5", "10", "15", "20") , color="black", size=3 , angle=0, fontface="bold", hjust=1) +


  	# Add labels on top of each bar

  	geom_text(data=dat, aes(x=seq, y=height+0.5, label=MetaboliteAnnot, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= dat$angle, inherit.aes = FALSE ) +


  	# Add base line information

  	geom_segment(data=base_data, aes(x = start, y = -0.5*23/18, xend = end, yend = -0.5*23/18), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )


ggsave(p2, file="Figure2A.Assoc.Lipidomics_updated.png", width=10, height=10)



#--------------------------------------------------------------
# for sub-figure B
#--------------------------------------------------------------

	dat = res.use[which(res.use$Figure2==1),]
	head(dat)
	dim(dat) # 117 272

	
	# 1. p values as height, color as association magnitudes and directions // winsorize the results at a limit to make sure meanningful results demonstration
	
	dat$height = -log10(dat$M2_All_FDR)
	dat$height[which(dat$height>15)]=15 # if more sig than 1e-15, draw at 1e-15
	
	dat$colors = dat$M2_All_Est
	dat$colors[which(abs(dat$colors)>0.3)] = sign(dat$colors[which(abs(dat$colors)>0.3)])*0.3 # if stronger than 0.3, draw at 0.3
	
	
	# 2. seperate super paths with a missing colum --> add 1 line after each Figure2.SuperP
	#    clear way for y-axis --> add 3 lines after the final Figure2.SuperP
	
	to_add = data.frame(matrix(NA,1*length(unique(dat$Figure2.SuperP)),ncol(dat)))
	colnames(to_add) = colnames(dat)
	to_add$Figure2.SuperP = unique(dat$Figure2.SuperP)
	to_add$Figure2.SubCat = "zzz_grid"
	
	to_add = to_add[order(to_add$Figure2.SuperP),]
	to_add = rbind(to_add,to_add[dim(to_add)[1],])
	to_add = rbind(to_add,to_add[dim(to_add)[1],])
	
	dat = rbind(dat,to_add)
	
	
	# 3. order data by effect size (M2_All_Est), defined seq (Figure2.seq), and super/sub pathway (Figure2.SubCat & Figure2.SuperP)

	dat = dat[order(dat$M2_All_Est),]
	dat = dat[order(dat$Figure2.seq),]
	dat = dat[order(dat$Figure2.SubCat),]
	dat = dat[order(dat$Figure2.SuperP),]
	
	dat$seq = 1:dim(dat)[1]
	

	# 4. angle of the direction of metab annot 
	
	angle = 90 - 360 * (dat$seq-0.5) /dim(dat)[1]
	dat$hjust = ifelse( angle < -90, 1, 0)
	dat$angle = ifelse(angle < -90, angle+180, angle)
	
	
	# 5. prepare a data frame for base lines
	
	dat$Figure2.SuperP = as.character(dat$Figure2.SuperP)
	cats = unique(dat$Figure2.SuperP)
	cats = cats[order(cats)]
	base_data = data.frame(group=cats,start=NA,end=NA,title=NA)
	
	for(i in 1:length(cats)) {
	  
	  substracti = dat[which(dat[,"Figure2.SuperP"]==cats[i]),]
	  base_data[i,"start"]=min(substracti$seq)
	  
	  if(i<length(cats)) {
	    base_data[i,"end"]=max(substracti$seq)-1 
	    
	  } else if (i==length(cats)) {
	    base_data[i,"end"]=max(substracti$seq)-3
	    
	  }
	  
	  base_data[i,"title"] = mean( c(base_data[i,"start"],base_data[i,"end"]) )
	  
	}
	
	base_data = base_data[order(base_data$start),]
	
	
	# 6. prepare a data frame for grid (scales)
	
	grid_data = data.frame(group=c("first","second"),start=NA,end=NA)
	grid_data[1,c("start","end")] = c(base_data[dim(base_data)[1],"end"]+0.5,base_data[dim(base_data)[1],"end"]+1.5)
	grid_data[2,c("start","end")] = c(max(dat$seq)-0.5,max(dat$seq)+0.5)
	

# ---------
# 7. plot figure
	
	p1= ggplot(dat) +
	  
	  # add the bars
	  
	  geom_bar( aes(x=seq, y=height, fill=colors), stat="identity") + 
	  coord_polar() +
	  
	  # add color patellae
	  
	  scale_fill_gradient2(low = "#006FC0", mid = "white", high = "#FF5552", midpoint = 0, na.value = "grey90", n.breaks=10) + 
	  ylim(c(-6,18)) + 
	  theme_light() + 
	  theme_minimal() +
	  
	  theme(
	    legend.position = "none",
	    axis.text = element_blank(),
	    axis.title = element_blank(),
	    panel.grid = element_blank(),
	    plot.margin = unit(rep(-1,4), "cm") 
	  ) +
	  
	  
	  # add value = 5,10,15 lines 
	  
	  geom_segment(data=grid_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
	  geom_segment(data=grid_data, aes(x = start, y = 5, xend = end, yend = 5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
	  geom_segment(data=grid_data, aes(x = start, y = 10, xend = end, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
	  geom_segment(data=grid_data, aes(x = start, y = 15, xend = end, yend = 15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
	  
	  
	  # Add text showing the value of each 100/75/50/25 lines
	  
	  ggplot2::annotate("text", x = rep(max(dat$seq)-0.5,5), y = c(0, 5, 10, 15, 20), label = c("0", "5", "10", "15", "20") , color="black", size=3 , angle=0, fontface="bold", hjust=1) +
	  
	  
	  # Add labels on top of each bar
	  
	  geom_text(data=dat, aes(x=seq, y=height+0.5, label=MetaboliteAnnot, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= dat$angle, inherit.aes = FALSE ) +
	  	  
	  
	  # Add base line information
	  
	  geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )
	
	
	ggsave(p1, file="Figure2B.Assoc.OtherMetabs_updated.png", width=10, height=10)
	
	
	
