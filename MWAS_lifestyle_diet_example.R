
#------------------------------------------------------------------------------------------------------
#
# Metabolome wide association analysis with T2D within a cohort
# Example code - name of cohort and cohort-specific data processing were removed/masked
#
#------------------------------------------------------------------------------------------------------

library(survival)
library(sas7bdat)
library(data.table)
library(mma)
library(readxl)


# load QCed, imputated, and batch normalized (inverse-normal-transformed) data

load(file="~/DataUse_Calculate.Food.r2.RData")

dim(all.int.use) # data table, row= participant, columns =metabolite and phenotypes
length(avai.metabs) # list of names of metabolites in this dataset


# inverse normal transformation by batch

int = function(var,batch) {
  
  dat = data.frame(y=var,b=batch,inty=NA)
  dat$y=as.numeric(var)
  dat$b=as.character(batch)

  batchn = unique(dat$b)
  for(i in 1:length(batchn)) {

       yi = as.numeric(dat$y[which(dat$b==batchn[i])])
       intyi = dat$inty[which(dat$b==batchn[i])]

       if(length(yi[!is.na(yi)])>0) {
              indexi = !is.na(yi)
              intyi[indexi] = qnorm( rank(yi[indexi]) / (length(yi[indexi])+1),mean=0,sd=1)
       }
       
       dat$inty[which(dat$b==batchn[i])] = intyi
  }

  return(dat$inty)

}

# var list

var.needed = c("id","cohort","ageyr","sex","current.smoking","actcont","BMIcont","fast","endpoint","caco",
                    "redm","processedm","totalfish","poultry","totaldairy","egg","vege","nutlegume","totalfruits",
                    "potato","wholegrain","totalrefined","sugardrinks","coffeetea","alcohol","ahei","calor")

outcome.needed = c("hrtbase","strbase","canbase","dbbase","deadbase","t2dbasediag","cancerbasediag","cvdbasediag",
                   "fhxdb","phxchol","phxhbp","WHRcont","mnpmh","multivtuse","aspirinuse","antiinfuse","antihtuse","antihluse",
                   "diabetes","ptime_diabetes","CVD","ptime_CVD","CHD","stroke","CVDdeath","ptime_CVDdeath")


all.int.use = all.qcd.use[,c(as.character(var.needed),as.character(outcome.needed),as.character(avai.metabs))]
head(all.int.use)
dim(all.int.use)


# inverse normal transformation metabolites by batch

for(i in 1:length(avai.metabs)) {

  metabi = avai.metabs[i]
  all.int.use[,metabi] = int(all.int.use[,metabi],all.int.use[,"batch"])

}

# for actcont and BMIcont - log and then standardize

std = function(trait) {
       return( (trait-mean(trait,na.rm=T))/sd(trait,na.rm=T) )
}

all.int.use$std.log_actcont = std(log(all.int.use$actcont))
all.int.use$std.log_BMIcont = std(log(all.int.use$BMIcont))


#-------------------------------------------------------------------------
# Run association analysis and calculate r2
#-------------------------------------------------------------------------

# emply result table for 7 models

res = as.data.frame(matrix(NA,length(avai.metabs),123))
rptnms = c("beta","se","t","p","var","r2")
names(res) = c("uID","metab_var", "n",
  paste("ageyr","_",rptnms,sep=''),      paste("sex","_",rptnms,sep=''),        paste("std.log_BMIcont","_",rptnms,sep='') , 
  paste("current.smoking","_",rptnms,sep=''), paste("std.log_actcont","_",rptnms,sep=''), 
  paste("redm","_",rptnms,sep=''),       paste("processedm","_",rptnms,sep=''), paste("totalfish","_",rptnms,sep='')   , 
  paste("poultry","_",rptnms,sep=''),    paste("totaldairy","_",rptnms,sep=''), paste("egg","_",rptnms,sep='')         , 
  paste("vege","_",rptnms,sep=''),       paste("nutlegume","_",rptnms,sep=''),  paste("totalfruits","_",rptnms,sep='') , 
  paste("potato","_",rptnms,sep=''),     paste("wholegrain","_",rptnms,sep=''), paste("totalrefined","_",rptnms,sep=''), 
  paste("sugardrinks","_",rptnms,sep=''),paste("coffeetea","_",rptnms,sep=''),  paste("alcohol","_",rptnms,sep='')     )
res$uID = avai.metabs

head(res)
dim(res) # 541 123

# var lista

foods = c("redm","processedm","totalfish","poultry","totaldairy","egg","vege","nutlegume","totalfruits",
    "potato","wholegrain","totalrefined","sugardrinks","coffeetea","alcohol")

covs = c("ageyr","sex","std.log_BMIcont","current.smoking","std.log_actcont")

othervar = c("id","cohort","fast","endpoint","caco")


# check distribution for outliers of foods

check = rbind(
  summary(all.int.use$redm),
  summary(all.int.use$processedm),
  summary(all.int.use$totalfish),
  summary(all.int.use$poultry),
  summary(all.int.use$totaldairy),

  summary(all.int.use$egg),
  summary(all.int.use$vege),
  summary(all.int.use$nutlegume),
  summary(all.int.use$totalfruits),
  summary(all.int.use$potato),

  summary(all.int.use$wholegrain),
  summary(all.int.use$totalrefined),
  summary(all.int.use$sugardrinks),
  summary(all.int.use$coffeetea),
  summary(all.int.use$alcohol))

rownames(check) = foods
print(check) 

# remove participants with missing data // in model, winsorize outliers

all.int.food = all.int.use[which( 
  !is.na(all.int.use$redm) & !is.na(all.int.use$processedm) & !is.na(all.int.use$totalfish) & !is.na(all.int.use$poultry) & 
  !is.na(all.int.use$totaldairy) & !is.na(all.int.use$egg) & !is.na(all.int.use$vege) & !is.na(all.int.use$nutlegume) & 
  !is.na(all.int.use$totalfruits) & !is.na(all.int.use$potato) & !is.na(all.int.use$wholegrain) & 
  !is.na(all.int.use$totalrefined) & !is.na(all.int.use$sugardrinks) & !is.na(all.int.use$coffeetea) & !is.na(all.int.use$alcohol)),]

all.int.food[1:10,1:10]
dim(all.int.food)


# run the regression for all metabolites
# given available metabolites in each batch are diff, models covs will depend on actual data wiht avail metab

for(i in 1:length(avai.metabs)) {

  # --------------
  # metabolite-specific data

  dati = all.int.food[,c(avai.metabs[i],foods,covs,othervar)]
  names(dati)[1] = "metabi"
  dati = dati[!is.na(dati$metabi),]

  dati$sex = ifelse(dati$sex=="women",1,0)
  dati$fast = ifelse(dati$fast=="nonfasting",0,1)

  # --------------
  # define models // + strata(endpoint) + strata(caco), data=dati

  m.base = "metabi ~ ageyr + current.smoking + std.log_BMIcont + std.log_actcont"

  # if more than one sex is in the data, adj sex

  if(length(unique(dati$sex))>1) {
    m.base = paste(m.base, " + sex")
  }

  # if more than one fasting status is in the data, adj fast

  if(length(unique(dati$fast))>1) {
    m.base = paste(m.base, " + fast")
  }

  # if more than one sub-study is in the data, adj endpoint and caco

  if(length(unique(dati$endpoint))>1) {
    m.base = paste(m.base, " + strata(endpoint)")
  }

  # see above

  if(length(unique(dati$caco))>1) {
    m.base = paste(m.base, " + strata(caco)")
  }


  # ------------------------------------------
  # model 2 to include all of foods
  # get model characteristics

  # create working data with winsorization - for all food - assgin values higher than cutoff to the cutoff

  for (k in 1:length(foods)) {
  
    cutoff = quantile(dati[,foods[k]],0.998,na.rm=TRUE) # cut off on original level
    dati[which(dati[,foods[k]]>cutoff),foods[k]] = cutoff

  }

  # get metab var and n
    
    res[i,"metab_var"]  = var(dati$metabi)
    res[i,"n"]  = dim(dati)[1]

  # get model estimate

    model = paste(m.base, " + ", paste(foods,collapse=" + "), sep='')
    
    if (length(unique(dati$sex))>1) {
       getvar = c(covs,foods)
       } else {
              getvar = setdiff(c(covs,foods),"sex")
       }
       
    allstats = coef(summary(lm( as.formula(model), data=dati)))[getvar,]

  for (k in 1:length(getvar)) {
    res[i,paste(getvar[k],"_",c("beta","se","t","p"),sep='')]  = allstats[getvar[k],]
    res[i,paste(getvar[k],"_var",sep='')] = var(dati[,getvar[k]])
  }

  print(i)

}


# calculate r2

getvar = c(covs,foods)

for(i in 1:length(getvar)) {
  res[,paste(getvar[i],"_r2",sep='')]   = res[,paste(getvar[i],"_beta",sep='')]^2  * res[,paste(getvar[i],"_var",sep='')]  / res$metab_var
}


# check results

head(res)
dim(res)


# write results

write.csv(res ,"Metabolite_covs_15foodgroups_assoc.r2_cohort_date.csv" ,row.names=FALSE)









