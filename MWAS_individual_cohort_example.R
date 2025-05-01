#------------------------------------------------------------------------------------------------------
#
# Metabolome wide association analysis with T2D within a cohort
# Example code - name of cohort and cohort-specific data processing were removed/masked
#
#------------------------------------------------------------------------------------------------------

library(survival)

# load QCed, imputated, and batch normalized (inverse-normal-transformed) data

lnames = load(file="~/Use_IncidentT2D.RData")

dim(annotation) # annotation to metabolites in data
dim(dat) # data - row =participants, columns =phenotypes and metabolites

identical(names(dat)[443:923],rownames(annotation)) # check if metabolites' order matched

# exclude the baseline T2D, cancers and CVD

dat = dat[which(is.na(dat$t2dbasediag) & is.na(dat$cancerbasediag) & is.na(dat$cvdbasediag)),]
head(dat)
dim(dat)

# check case numbers: diabetes = incident T2D

table(dat$diabetes)

# emply result table for 7 models

res = data.frame(matrix(NA,dim(annotation)[1],29))
names(res) = c("HMDB","Est1","sem1","P1","FDR1","Est2","sem2","P2","FDR2","Est3","sem3","P3","FDR3",
  "Est4","sem4","P4","FDR4","Est5","sem5","P5","FDR5",
  "Est6","sem6","P6","FDR6","Est7","sem7","P7","FDR7")
res$HMDB = as.character(annotation$HMDB)
res = merge(annotation,res,by="HMDB")
res$HMDB = as.character(res$HMDB)

# run association
# diabetes is incident T2D, 
# ptime_diabetes is time to event
# labcodes/caco are sub-studies and case-control within sub-study

covs = c("ageyr","smoking","alcohol","fast.bloodq","antihluse","antihtuse","fhxdb","phxhbp","phxchol","BMIcont","WHRcont","actcat","ahei2010","sbp","dbp","multivtuse","mnpmh","adj_TChol","adj_Trig","adj_LDLC","diabetes","ptime_diabetes","eGFR")


for(i in 1:dim(res)[1]) {

  d = dat[,c(res$HMDB[i],"diabetes","ptime_diabetes",covs,"caco","labcode")]
  names(d)[1] = "Metab"
  d=d[!is.na(d$Metab),]

  fit1 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp + phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode), data=d)))

  fit2 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp +  phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode) + BMIcont + WHRcont, data=d)))

  fit3 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp + phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode) + actcat + ahei2010, data=d)))

  fit4 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp + phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode) + adj_TChol + adj_Trig + adj_LDLC, data=d)))

  fit5 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp + phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode) + sbp + dbp, data=d)))

  fit6 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp + phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode) + eGFR, data=d)))

  fit7 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Metab + ageyr + as.factor(smoking) + alcohol + fast.bloodq + antihluse + antihtuse + fhxdb + phxhbp +  phxchol + as.factor(mnpmh) + strata(caco) + strata(labcode) + BMIcont + WHRcont + eGFR, data=d)))

  res[i,c("Est1","sem1","P1")] = fit1["Metab",c(1,3,5)]
  res[i,c("Est2","sem2","P2")] = fit2["Metab",c(1,3,5)]
  res[i,c("Est3","sem3","P3")] = fit3["Metab",c(1,3,5)]
  res[i,c("Est4","sem4","P4")] = fit4["Metab",c(1,3,5)]
  res[i,c("Est5","sem5","P5")] = fit5["Metab",c(1,3,5)]

  res[i,c("Est6","sem6","P6")] = fit6["Metab",c(1,3,5)]
  res[i,c("Est7","sem7","P7")] = fit7["Metab",c(1,3,5)]

  print(i)

}

# fdr-adjustment within cohort

res$FDR1 = p.adjust(res$P1, method='fdr', n=length(res$P1))
res$FDR2 = p.adjust(res$P2, method='fdr', n=length(res$P2))
res$FDR3 = p.adjust(res$P3, method='fdr', n=length(res$P3))
res$FDR4 = p.adjust(res$P4, method='fdr', n=length(res$P4))
res$FDR5 = p.adjust(res$P5, method='fdr', n=length(res$P5))

res$FDR6 = p.adjust(res$P6, method='fdr', n=length(res$P6))
res$FDR7 = p.adjust(res$P7, method='fdr', n=length(res$P7))

head(res)
dim(res)


# check on the results in main models 

res[which(res$FDR1<0.05),c("metabolite_name","Est1","FDR1")]
dim(res[which(res$FDR1<0.05),])

res[which(res$FDR2<0.05),c("metabolite_name","Est2","FDR2")]
dim(res[which(res$FDR2<0.05),])

cor(res$Est1,res$Est6,use="complete.obs")
cor(res$Est2,res$Est7,use="complete.obs")


# write results

write.csv(res,"~/MetabolomeWide_T2D_cohort_date.csv",row.names=F)



