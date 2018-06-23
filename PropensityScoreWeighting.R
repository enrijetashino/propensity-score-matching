#----------------------------#
# Propensity Score Weighting #
# Author: Enrijeta Shino     #
# Last Update: Jan 10, 2017  #
#----------------------------#

#-----------------------------------------------------------------------------------------------------#
# Replication code for "Using Propensity Score Matching to Determine the Efficacy of Secondary Career #
# Academies in Raising Educational Aspirations" by Jay W. Rojewski, In Heok Lee, Sinan Gemici         #
#-----------------------------------------------------------------------------------------------------#

#load required survey library
require(survey)


#normalize the base year student sampling weight so that it sums to the sample size
ELS.data.imputed$bystuwt = ELS.data.imputed$bystuwt/mean(ELS.data.imputed$bystuwt)
#bystuwt this is the sampling weight 
#imputed$bystuwt/mean divideing the weight by the mean it is called "the normalization" of the weight

#define design object that describes the sample characteristics
#the variable psu identifies the primary sampling units (cluster ids)
#the variable STRATA_ID identifies the strata ids
#the variable bystuwt identifies the base-year sampling weights for
#respondents of the 2002 and 2004 waves (Base year and 1st follow-up)
surveyDesign = svydesign(ids=~psu, strata=~STRAT_ID, weights=~bystuwt,
                         data = ELS.data.imputed, nest=T)
#svydesign defines the object 
#since the data are stratified we can use multilevel analysis as well, but we're not going to becuase there #is no need for it in this case 

#get weighted percentages of treated and untreated cases
#the treatment variable is "BYS33K" - Ever in career academy from BY Student Questionnaire
(treatmentTable = svymean(~BYS33K, surveyDesign)) #8.8% of cases received treatment
#0.08 of the pop received career academy treatment. 

table(ELS.data.imputed$BYS33K)
mean(ELS.data.imputed$BYS33K)
#-------------------------------------------------

#obtain weights for estimating the ATT for propesnity scores obtained with logistic regression
ELS.data.imputed$weightATT = with(ELS.data.imputed,ifelse(BYS33K==1, 
                                                          1, pScores/(1-pScores))) 
with(ELS.data.imputed, by(weightATT,BYS33K,summary))
#"with" allows to put the name of the variables without adding the name of the dataset
#the minimum should never be 0

#obtain weights for estimating the ATE
ELS.data.imputed$weightATE = with(ELS.data.imputed,ifelse(BYS33K==1, 
                                                          1/pScores, 1/(1-pScores)))
#"ifelse" adds the condition

with(ELS.data.imputed, by(weightATE,BYS33K,summary))
#we run the summary with this command
#this weight seems to be too strong(or extreme) the max for treated is equal to 145.1.
#this can happen due to misspecification, lack of the common support 
#one way to check for extreme weight we have to check the covariate balance first 
#one conseqeuence of the extreme weight due to misspecification is failure to get covariate balance
#and the standard errors get inflated, making it difficult for us to find a treatment effect 

#======================================================

#STRATEGIES TO DEAL WITH EXTREME WEIGHTS

#TRUNCATION OPTION
#check the 99th percentile weight:
with(ELS.data.imputed,quantile(weightATE, 0.99))
#check how many observatios will have truncated weights
with(ELS.data.imputed,sum(weightATE>quantile(weightATE, 0.99)))
#truncate weights at the 99th percentile:
ELS.data.imputed$weightATETruncated = with(ELS.data.imputed,
                                           ifelse(weightATE > quantile(weightATE, 0.99), 
                                                  quantile(weightATE, 0.99),weightATE)) 
#use "ifelse" to include truncation 

#check truncated weights
with(ELS.data.imputed, by(weightATETruncated,BYS33K,summary))


#USE STABILIZAED WEIGHTS
#obtain stabilized weights for the estimation of ATE
#need to calcultae the constants first 
ELS.data.imputed$C=with(ELS.data.imputed,ifelse(BYS33K==1,pScores,1-pScores))
surveyDesign = svydesign(ids=~psu, strata=~STRAT_ID, weights=~bystuwt,
                         data = ELS.data.imputed, nest=T) #recreate survey design

(constants = svyby(~C, by=~BYS33K, design=surveyDesign, FUN=svymean))
#"surveyby" gives the constant 
#double paranthesis (one in the beggining and oen in the end of the command) in the end of the command is a #trick to print the output otherwise it won't print anything

ELS.data.imputed$stabilizedWeightATE = ifelse(ELS.data.imputed$BYS33K==1, 
                                              constants[1,2]/ELS.data.imputed$C, 
                                              constants[2,2]/ELS.data.imputed$C)

#check stabilized weigts for extremeness
with(ELS.data.imputed, by(stabilizedWeightATE,BYS33K,summary))

#compare the distributions of weights and stabilized weights for the ATE
tiff("Chapter3_figure3-1.tif", res=600, compression = "lzw", height=6, width=15, units="in")
with(ELS.data.imputed,
     hist(weightATE[BYS33K==1], density = 10, angle = 45, main="Propensity Score Weights for the ATE",
          breaks=seq(0,150,by=1),
          xlim=c(0,150), ylim=c(0,100), xlab="Shaded = IPTW | Gray = Stabilized IPTW") )
with(ELS.data.imputed,
     hist(stabilizedWeightATE[BYS33K==1], col=gray(0.4,0.25), breaks=seq(0,150,by=1),
          xlim=c(0,150), ylim=c(0,100),add=T) )
dev.off()

#=====================================================
#WEIGHTS WITH PROPENSITY SCORE 
#obtain weights with propensity scores estimated with random forests

#obtain weights for estimating the ATT for propesnity scores obtained with logistic regression
ELS.data.imputed$weightATTRf = with(ELS.data.imputed,ifelse(BYS33K==1, 
                                                            1, pScoresRf/(1-pScoresRf))) 
#pScoresRf this is for the random forest 

with(ELS.data.imputed, by(weightATTRf,BYS33K,summary))

#obtain weights for estimating the ATE
ELS.data.imputed$weightATERf = with(ELS.data.imputed,ifelse(BYS33K==1, 
                                                            1/pScoresRf, 1/(1-pScoresRf)))

with(ELS.data.imputed, by(weightATERf,BYS33K,summary))


#======================================================

#obtain weights with propensity scores estimated with GBM

#obtain weights for estimating the ATT for propesnity scores obtained with logistic regression
ELS.data.imputed$weightATTGBM = with(ELS.data.imputed,ifelse(BYS33K==1, 
                                                             1, pScoresGBM[,1]/(1-pScoresGBM[,1]))) 

#pScoresGBM weight for the generalized boosted model
with(ELS.data.imputed, by(weightATTGBM,BYS33K,summary))

#obtain weights for estimating the ATE
ELS.data.imputed$weightATEGBM = with(ELS.data.imputed,ifelse(BYS33K==1, 
                                                             1/pScoresGBM[,1], 1/(1-pScoresGBM[,1])))


with(ELS.data.imputed, by(weightATEGBM,BYS33K,summary))

#we used only one imputed dataset for all the analysis above

#==========================================================
#obtain weights for multiple imputed datasets
library(mitools)
#add propensity score mean to imputed datasets
allImputations <- update(allImputations, weightATT = ifelse(BYS33K==1, 
                                                            1, pScores/(1-pScores)))

#update gets all the imputed dataset and applies it to them
with(allImputations, by(weightATT,BYS33K,summary))
# this analysis ahs to be doublechecked because the values are exactly the same and that seems to be a #problem. The values should be similar but not equal 

#======================================================
#EVALUATE COVARIATE BALANCE
#create final weights for the base year
#the final weigth is the product of the propensity score weight
#and the sampling weight
ELS.data.imputed$finalWeightBY <- with(ELS.data.imputed,bystuwt*weightATT)

#evaluate covariate balance for ATT
require(twang)
#twang is used gor GBM but it is also good for obtaining tables 
balanceTable <- bal.stat(ELS.data.imputed, vars= covariateNames,
                         treat.var = "BYS33K", 
                         w.all = ELS.data.imputed$finalWeightBY, get.ks=F, #w.all is the final weight
                         sampw = ELS.data.imputed$bystuwt, #sampw is the sampling weight 
                         estimand="ATT", multinom=F)
#Table with results of balance evaluation. The columns are:
# tx.mn - The mean of the treatment group
# tx.sd	- The standard deviation of the treatment group
# ct.mn - The mean of the control group
# ct.sd	- The standard deviation of the control group
# std.eff.sz	- The standardized effect size, (tx.mn-ct.mn)/tx.sd. (this should be less than -1, or -.05)
# stat	- the t-statistic for numeric variables and the chi-square statistic for categorical variables
# p	- the p-value for the test associated with stat (we cannot use this becuase it is affected by the #sample size)

balanceTable <- balanceTable$results 
#calculate variance ratio (variance of the treated divided by the varaince of the control)
balanceTable$varRatio <- with(balanceTable, tx.sd^2/ct.sd^2) 

#summarize the covariate balance quality
(summaryBalance <- rbind(
  std.eff.sz = summary(abs(balanceTable$std.eff.sz)), #standardized effect sizes
  varRatio = summary(balanceTable$varRatio) ))#variance ratios
#the max in the std.eff is equal to .01 which is very good 


#==============================================================
#create a function to return balance summaries for this particular dataset
balanceSummarizer = function(data=ELS.data.imputed, samplingWeight="bystuwt", 
                             PSWeight, treatment="BYS33K",
                             effect="ATT", covariateNames) {
  
  
  finalWeight <- data[,samplingWeight]*data[,PSWeight]
  
  #evaluate covariate balance for ATT
  require(twang)
  balance.table = bal.stat(data, vars= covariateNames,
                           treat.var = treatment, 
                           w.all = finalWeight, get.ks=F,
                           sampw = data[,samplingWeight], 
                           estimand=effect, multinom=F)
  
  balance.table = balance.table$results 
  
  #calculate variance ratio
  balance.table$varRatio <- with(balance.table, tx.sd^2/ct.sd^2) 
  
  #summarize the covariate balance quality
  return(rbind(
    std.eff.sz = summary(abs(balance.table$std.eff.sz)), #standardized effect sizes
    varRatio = summary(balance.table$varRatio) ))#variance ratios
  
  #close function
}

#this loop is a nice and elegant piece but not required for the newbeies in R

#====================================================================
#obtain summaries of covariate balances with different propensity scores

Table3.1 <- rbind(
  logistic = balanceSummarizer(PSWeight="weightATT",covariateNames=covariateNames),
  RF = balanceSummarizer(PSWeight="weightATTRf",covariateNames=covariateNames),
  GBM = balanceSummarizer(PSWeight="weightATTGBM",covariateNames=covariateNames))

Table3.1 = data.frame(Table3.1, 
                      index=rep(c("std.eff.sz","varRatio"),3), method=rep(c("Logistic","RF","GBM"),each=2))
write.csv(Table3.1, file="Table3-1.csv")

Table3.2 <- rbind(
  logistic = balanceSummarizer(PSWeight="weightATE",covariateNames=covariateNames),
  truncated = balanceSummarizer(PSWeight="weightATETruncated",covariateNames=covariateNames),
  stabilized = balanceSummarizer(PSWeight="stabilizedWeightATE",covariateNames=covariateNames),
  RF = balanceSummarizer(PSWeight="weightATERf",covariateNames=covariateNames),
  GBM = balanceSummarizer(PSWeight="weightATEGBM",covariateNames=covariateNames))

Table3.2 = data.frame(Table3.2, 
                      index=rep(c("std.eff.sz","varRatio"),5), 
                      method=rep(c("Logistic","Truncated","Stabilized","RF","GBM"),each=2))
write.csv(Table3.2, file="Table3-2.csv")

#use different methods and compare them to find the one that does better

#===============================================================
#create a final weights for estimating the effect of career academy participation
#on income in the second folow up (2006) 
ELS.data.imputed$finalWeight2006 = with(ELS.data.imputed,F2BYWT*weightATT)


#normalize the base-year to second folow-up (2006) final weight
ELS.data.imputed$finalWeight2006 = ELS.data.imputed$finalWeight2006/mean(ELS.data.imputed$finalWeight2006)


#check distribution of final weights
summary(ELS.data.imputed$finalWeight2006)


#===========================================================
#ESTIMATE THE ATT
#re-create the survey design including the final weight for 2006
surveyDesign2006 <- svydesign(ids=~psu, strata=~STRAT_ID, weights=~finalWeight2006,
                              data = ELS.data.imputed, nest=T)
#we do this to clarify the weight
#this is multilevel data but we don't want to use multilevel modeling, "psu" is the cluster 
#this creates the survey design 

#create replicate weights for bootstrapping
surveyDesign2006Boot = as.svrepdesign(surveyDesign2006, type=c("bootstrap"),replicates=1000)
#this one will create the same design but with a 1000 bootrstaps 
#we could also substitute "bootstrap" with "jackknife" 


#obtain ATT as weighted mean differences.
(weightedMeans=svyby(formula=~F2ERN5P2,by=~BYS33K,design=surveyDesign2006Boot, 
                     FUN=svymean,covmat=TRUE))
(ATT2006 = svycontrast(weightedMeans, contrasts=c(-1,1)))
#svyby funstion to breakdown the estimates by using 
#covmat- surevy covariate matrix 
#svycontrast- to get the difference between means 

pnorm(ATT2006[1]/c(ATT2006)[2])
pnorm(0.23395/0.0845)
pnorm(0.23395/0.0845, lower.tail=F)
#this is the p-value 

#obtain the group variances
(weightedVars=svyby(formula=~F2ERN5P2,by=~BYS33K,design=surveyDesign2006Boot, 
                    FUN=svyvar,covmat=TRUE))
#use the same strategy to get group variances 

#estimate the ATT for 2006 with regression analysis for complex survey data
outcomeModel2006 = svyglm(F2ERN5P2~BYS33K,surveyDesign2006)
summary(outcomeModel2006)
#glm model 
#the two estimators are the same, only SE are different 


#re-estimate the ATT with regression, but this time obtain standard errors with bootstrapping
outcomeModel2006Boot = svyglm(F2ERN5P2~BYS33K,surveyDesign2006Boot)
summary(outcomeModel2006Boot)
#=====================================================================
#Estimate the ATT with multiple imputed datasets
#here we repeat the same but with multiple imputed dataset 

#first, add the final weight to imputed datasets
allImputations <- update(allImputations, finalWeightATT = F2BYWT*weightATT)
#normalize the final weight
allImputations <- update(allImputations, finalWeightATT = finalWeightATT/mean(finalWeightATT))

#create survey design object for multiple imputed dastests
surveyDesign2006MI <- svydesign(ids=~psu, strata=~STRAT_ID, weights=~finalWeightATT,
                                data = allImputations, nest=T)

#estimate the ATT with regression in each imputed dataset
outcomeModel2006MI = with(surveyDesign2006MI, svyglm(F2ERN5P2~BYS33K))

#combine estimates from multiplle imputed datasets
resultsModel2006MI = MIcombine(outcomeModel2006MI)
summary(resultsModel2006MI)
#MIcombine- will combine the estiamtes obtained by imputed datssets 


#=====================================================================
#Doubly-robust (DR) estimation with regression estimation of 
#propensity score weighted means using the propensity score as a covariate
#this is without covariates, even though would be good to include covatiates 

#estimate a regression of the outcome with the propensity score as a covariate
outcomeModelDR = svyglm(F2ERN5P2~BYS33K+pScores+I(pScores^2)+I(pScores^3),surveyDesign2006Boot)

#Estimate ATT with a doubly-robust method
weightedMeansDR = predict(outcomeModelDR,
                          newdata=data.frame(BYS33K=0:1, pScores=0.125), vcov=TRUE, type="response")
(ATT2006DR = svycontrast(weightedMeansDR, c(-1,1)))
#predict- predicts the two means given the model 
#svycontrast- to get the treatment effect 

