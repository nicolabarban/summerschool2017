scores <-read.table("Education_SCORES_AT_ALL_THRESHOLDS.txt", header=T)

sample <- read.table("../data/hapmap_EA.fam", col.names=c("FID","IID", "motherid", "fatherid", "sex", "education" ))


data<- merge(scores, sample, by="IID")

hist(data$pT_1)

summary( lm(education~pT_1, data=data))
summary( lm(education~pT_1, data=data))$r.squared


summary( lm(education~ pT_0.00000005 , data=data))
summary( lm(education~ pT_0.00000005 , data=data))$r.squared
