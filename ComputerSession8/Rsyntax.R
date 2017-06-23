


install.packages(c("interplot", "foreign", "ggplot2"), dependencies=T)
library(foreign)
library(interplot)
library(ggplot2)

mydata <- read.dta("../data/GxEdata.dta")

head(mydata)


### ### ### ### ### ### ### ### ### ### 
### ### ### ### HEIGHT
summary(height)
hist(height)
model_height <- lm(height ~ score_height, data = mydata)
summary(model_height)

model_heightXbirthyear <- lm(height ~ score_height + byear + score_height*byear, data = mydata)
summary(model_heightXbirthyear)

### ### ### ### look at interaction
### ### how does coef change across birth cohorts
ggplot(mydata, aes(x=score_height ,y=height, col=gender))   + 
  					geom_point(alpha = 0.4) + ylim(.93,2.2) 
   				 	stat_smooth(method = "lm", col = "black") 
   
   
ggplot(mydata, aes(x=score_height ,y=height, col=gender))   + 
				geom_point(alpha = 0.4) + ylim(.93,2.2) +
				geom_smooth(method="lm", fill=NA)

#interplot(m = model_heightXbirthyear, var1 = "height", var2 = "score_height") +
#  xlab("Year born") +
 # ylab("Estimated Coefficient for Height genetic score") +
 # theme_bw() +
 # ggtitle("GxE interaction between polygenic score for Height and birth cohort") +
 # theme(plot.title = element_text(face="bold")) 

### ### ### ### look at interaction
### ### how does coef change across birth cohorts
theme_set(theme_bw())
ggplot(mydata, aes(x = score_height, y = height, color = byear_dum, shape = byear_dum)) + 
  stat_smooth(method = 'lm') +
  xlab("Genetic score for Height") +
  ylab("Height") +
  ggtitle("GxE interaction between polygenic score for Height and birth cohort") +
  theme(plot.title = element_text(face="bold")) 

### ### ### ### ### ### ### ### ### ### 
### ### ### ### Education


model_edu <- lm(educyears ~ score_educ, data = mydata)
summary(model_edu)

model_eduXbirthyear <- lm(educyears ~ score_educ + byear + score_educ*byear, data = mydata)
summary(model_eduXbirthyear)

png("images/ea_inter1.png")

interplot(m = model_eduXbirthyear, var1 = "byear", var2 = "score_educ") +
  xlab("Year born") +
 ylab("Estimated Coefficient for Education genetic score") +
 theme_bw() +
  ggtitle("GxE interaction between PGS for EA and birth cohort") +
  theme(plot.title = element_text(face="bold")) 
  dev.off()
  
  
  
theme_set(theme_bw())
ggplot(mydata, aes(x = score_educ, y = educyears, color = byear_dum, shape = byear_dum)) + 
  stat_smooth(method = 'lm') +
  xlab("Genetic score for EA") +
  ylab("Educational Attainment") +
  ggtitle("GxE interaction between polygenic score for EA and birth cohort") +
  theme(plot.title = element_text(face="bold")) 
  
  
### ### ### ### ### ### ### ### ### ### 
### ### ### ### bmi

model_bmi <- lm(bmi ~ score_bmi, data = mydata)
summary(model_bmi)

model_bmiXbirthyear <- lm(bmi ~ score_bmi + byear + score_bmi*byear, data = mydata)
summary(model_bmiXbirthyear)


interplot(m = model_bmiXbirthyear, var1 = "byear", var2 = "score_bmi") +
  xlab("Year born") +
 ylab("Estimated Coefficient for BMI genetic score") +
 theme_bw() +
  ggtitle("GxE interaction between PGS for BMI and birth cohort") +
  theme(plot.title = element_text(face="bold")) 
  dev.off()
  
  
  
theme_set(theme_bw())
ggplot(mydata, aes(x = score_bmi, y = educyears, color = byear_dum, shape = byear_dum)) + 
  stat_smooth(method = 'lm') +
  xlab("Genetic score for BMI") +
  ylab("Educational Attainment") +
  ggtitle("GxE interaction between polygenic score for BMI and birth cohort") +
  theme(plot.title = element_text(face="bold")) 
   

