#Bimodal spatial localization task 
rm(list = ls())
setwd("~/Desktop/NYU/Project2/Experiment code/Stats/unityJudgment")
library(lme4)
library(car)
library(readxl)
library(effects)
library(ggplot2)
library(emmeans)

unity_indvd_data = read_excel('unityJdg_binaryData_overall.xlsx')
unity_indvd_data$SpatialD_abs<-factor(unity_indvd_data$SpatialD_abs, ordered=TRUE)
GLMMmodel_AV <- glmer(UnityJudgment ~ 1 + SpatialD_abs + Condition + Phase + 
                        SpatialD_abs * Condition * Phase  + (1|SubjI), 
                      data = unity_indvd_data, family='binomial', 
                      control = glmerControl(optimizer='bobyqa'))
summary(GLMMmodel_AV, corr = FALSE)
Anova(GLMMmodel_AV)
confint(GLMMmodel_AV)

################################################################################
rm(list = ls())
VE_indvd_data <- read.table("locShifts_indvdTrials_overall.txt", header=TRUE,sep = ",")
A_VE_indvd_data = VE_indvd_data[VE_indvd_data$Modality == 'A',]
V_VE_indvd_data = VE_indvd_data[VE_indvd_data$Modality == 'V',]
#Bimodal spatial localization task (ventriloquism effect, auditory shifts towards the visual stimulus)
A_VE_indvd_data$SpatialD_abs<-factor(A_VE_indvd_data$SpatialD_abs, ordered=TRUE)
A_VE_indvd_data$SpatialD<-factor(A_VE_indvd_data$SpatialD, ordered=TRUE)
A_VE_indvd_data$SubjI<-factor(A_VE_indvd_data$SubjI)
A_VE_indvd_data$Condition<-factor(A_VE_indvd_data$Condition)
A_VE_indvd_data$Phase<-factor(A_VE_indvd_data$Phase)

lmer_resultsA<-lmer(LocShifts_corr~ (Condition*Phase*SpatialD_abs) + (1|SubjI), data=A_VE_indvd_data)
summary(lmer_resultsA)
Anova(lmer_resultsA)
confint(lmer_resultsA)

################################################################################
#Bimodal spatial localization task (ventriloquism effect, auditory shifts towards the visual stimulus)
V_VE_indvd_data$SpatialD_abs<-factor(V_VE_indvd_data$SpatialD_abs, ordered=TRUE)
V_VE_indvd_data$SpatialD<-factor(V_VE_indvd_data$SpatialD, ordered=TRUE)
V_VE_indvd_data$SubjI<-factor(V_VE_indvd_data$SubjI)
V_VE_indvd_data$Condition<-factor(V_VE_indvd_data$Condition)
V_VE_indvd_data$Phase<-factor(V_VE_indvd_data$Phase)

lmer_resultsV<-lmer(LocShifts_corr ~ (Condition*Phase*SpatialD_abs) + (1|SubjI), data=V_VE_indvd_data)
summary(lmer_resultsV)
Anova(lmer_resultsV)
confint(lmer_resultsV)

#post-hoc ttests
V_VE_indvd_data_cong <- V_VE_indvd_data[V_VE_indvd_data$Condition == 'congruent',][c('LocShifts_corr','SubjI')]
V_VE_indvd_data_incong <- V_VE_indvd_data[V_VE_indvd_data$Condition == 'incongruent',][c('LocShifts_corr','SubjI')]
subj_init <- levels(factor(unique((V_VE_indvd_data$SubjI))))
V_VE_data_cong_avg <- rep(NA, length(subj_init))
V_VE_data_incong_avg <- rep(NA, length(subj_init))
for (s in 1:length(subj_init))
{
  cong_data_slc = V_VE_indvd_data_cong[V_VE_indvd_data_cong['SubjI']==subj_init[s],]['LocShifts_corr']
  V_VE_data_cong_avg[s] <- mean(cong_data_slc$LocShifts_corr, na.rm = TRUE)
  incong_data_slc = V_VE_indvd_data_incong[V_VE_indvd_data_incong['SubjI']==subj_init[s],]['LocShifts_corr']
  V_VE_data_incong_avg[s] <- mean(incong_data_slc$LocShifts_corr, na.rm = TRUE)
}
t.test(V_VE_data_cong_avg, V_VE_data_incong_avg, alternative = "two.sided", var.equal = TRUE)

ggplot(data.frame(Effect(c("SpatialD_abs","Condition"),lmer_resultsV)), 
       aes(x=factor(SpatialD_abs, level = c(0,8,16,24)),y=fit,
           color=Condition,group=Condition))+ geom_line()
