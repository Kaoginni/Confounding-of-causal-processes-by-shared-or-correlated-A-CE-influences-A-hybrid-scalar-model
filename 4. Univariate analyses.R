# Script Name: Univariate Sex Differences models incorporating Singetons
# Description: This script tests for quantitative and qualitative sex differences as well as variance-inequality (scalar) in variance components 
# This is based on Neale et al., Multivariate genetic analysis of sex-lim and G x E interaction, Twin Research & Human Genetics, 2006
# And modified from original script provided by Fruhling Rijsdijk.
# The model is specified to deal with the DZ opposite sex pairs 

# Bivariate twin analysis models to estimate causes on (co)variation (ACE) with heterogeneity in ACE paths for males and females
# Phenotypes:	SES, Financial strain, Highest level of education, Western Dietary pattern, Prudent dietary pattern, Depessive symptoms, Anxiety symptoms
# Covariates regressed out: age, sex
# Heterogeneity variable:	sex
# Zygosity variable:		c1c2_sexzyg (1=MZM, 2=DZM, 3=MZF, 4=DZF, 5=DZOS, 6=Singletons, 7=Triplets)

# MODELS (SPECIFIED SEPARATELY FOR EACH OF THE VARIABLES OF INTEREST):

#-----------------------------------------------------------------------
# Model II-1:
# Quantitative sex differences for A, C, and E
# In this model, there are male and female paths for mzm and dzm and mzf and dzf respectively
# For dzos twins, male paths and female paths are combined in one model
# In this model, rAmf and rCmf are constrained to 0.5 and 1 respectively
#-------------------------------------------------------------------------------------------------------------------------------
# Model II-2:
# In this model, we call Model 1 but allow rCmf to be free. THIS IS THE Quantitative model + Qualitative SEX DIF FOR C MODEL
#-------------------------------------------------------------------------------------------------------------------------------
# Model II-3:
# In this model, we call Model 1 but I allow rAmf to be free. THIS IS THE Quantitative model + Qualitative SEX DIF FOR A MODEL
#--------------------------------------------------------------------------------------------------------------------------------
# Model II-4: No qualitative or quantitative differences - Homogeneity CALLING MODEL 1
# In this model, af, cf and ef are constrained to equal am, cm, and em respectively
#----------------------------------------------------------------------------------------------------
# Model II-5: Variance in-equality sex dif model 
# In this model, THERE IS ONLY ONE SET OF A, C AND E PARAMETERS, BUT WE ALLOW THE VARIANCE OF THE 
# FEMALES TO BE ESTIMATED AS A SCALAR * THE VARAINACE OF THE MALES
#---------------------------------------------------------------------------------------------


#Clear workspace
rm(list=ls())

library(mvtnorm)
library(psych)
library(foreign)
require(plyr)
library(OpenMx) 
require(Hmisc)	
library(dplyr)

mxOption(NULL, "Default optimizer", "CSOLNP") # SLSQP is a better optimizer for ordinal data

# -------------------------------------------------------------------------------------------------------------
# 1. Read data into R and prep
# -------------------------------------------------------------------------------------------------------------

ALLdata 	<- read.dta("3. USE - 2017.02.20_COTASS2_Q_and_Clinical_v11.dta")
str(ALLdata)
class(ALLdata) 		# to check that the data is loaded correctly, should return "data.frame"
dim(ALLdata)      	# should be 3969 obs., 852 variables
psych::describe(ALLdata)	# this will show you all the variables in Cotass-2

ALLdata$sex	<-as.numeric(ALLdata$sex)
ALLdata$zyg	<-as.numeric(ALLdata$zyg)
ALLdata$old_zygosity	<-as.numeric(ALLdata$old_zygosity)
ALLdata$sing_or_twin	<-as.numeric(ALLdata$sing_or_twin)
ALLdata$c1c2_sexzyg	<-as.numeric(ALLdata$c1c2_sexzyg)
ALLdata$sexzyg	<-as.numeric(ALLdata$sexzyg)
ALLdata$finstrain	<-as.numeric(ALLdata$finstrain)

psych::describe(ALLdata)	# this will show you all the variables in Cotass-2

#Remove problem participants
ALLdata 	<- ALLdata[ALLdata$res_code!='51941' & ALLdata$res_code!='51942',]  #gender info incorrect
dim(ALLdata)      #should be 3967 obs.

# Generate a new "numtw.noDZOS" (= twin in which males are always first in DZOS pairs, but also including the singletons:
# I already did this in SPSS, but I coded singletons to 3 to include the third of the triplet trio

#ALLdata$TwNo <- ALLdata$numtw_noDZOS; describe(ALLdata$TwNo); describe(ALLdata$numtw_noDZOS)
#ALLdata$TwNo[ALLdata$sing_or_twin=="Singleton"]=1;  
#table(ALLdata$numtw_noDZOS); table(ALLdata$TwNo)

# Select only variables of interest and save them in a smaller file
# In this example apart from all necessary variables we select Diastolic BP, Systolic BP and BMI (just as an example)

subdata1 	<- ALLdata[,c("numclust","TwNo","zyg","sexzyg","age_reported","sex","Rec_Edu","Rec_el_count","Rec_el_tr_count",
				"finstrain","Rec_Tot_Diet_dly","Rec_Unhealthy_Diet_dly","Rec_Healthy_Diet_dly","Dep","Anx")] 

# Rename variables in subdata1 to make it easier. Note CMR stands for Cardio-Metabolic Risk & AM for Anthropometric (generic terms so that you don't 
# need to change the rest of the script if you select different variables 

names(subdata1) <- c("FamID","TwNo","zyg","sexzyg","age","sex","edu","elect","elect_trans",
				"finstrain","totdiet","unhealthdiet","healthdiet","Dep","Anx")

dim (subdata1)
summary(subdata1)
psych::describe(subdata1)

# Because we don't want triplets, remove lines in which TwNo==3; this removes the 3rd sibling in each triplet group
#subdata2 <- subdata1[subdata1$TwNo!=3,]
#nrow(subdata1)
#nrow(subdata2)


# -------------------------------------------------------------------------------------------------------------------
# 2. Regress out covariates, transform (if necessary) and Reshape the data to get pair-wise file for genetic analyses
# -------------------------------------------------------------------------------------------------------------------
# Regress out age, sex (This is standard procedure in twin modeling)
## edu
subdata1$edu_resid <- ((residuals(lm(subdata1$edu ~ subdata1$age + subdata1$sex, na.action="na.exclude")))*1.2)+4.6
psych::describe(subdata1$edu)
psych::describe(subdata1$edu_resid)

# Check distributions
hist(subdata1$edu)
hist(subdata1$edu_resid)

## elect
subdata1$elect_resid <- ((residuals(lm(subdata1$elect ~ subdata1$age + subdata1$sex, na.action="na.exclude")))*1.4)+3.2
psych::describe(subdata1$elect)
psych::describe(subdata1$elect_resid)

# Check distributions
hist(subdata1$elect)
hist(subdata1$elect_resid)

## elect_trans
subdata1$elect_trans_resid <- ((residuals(lm(subdata1$elect_trans ~ subdata1$age + subdata1$sex, na.action="na.exclude")))*0.6)+3.2
psych::describe(subdata1$elect_trans)
psych::describe(subdata1$elect_trans_resid)

# Check distributions
hist(subdata1$elect_trans)
hist(subdata1$elect_trans_resid)

## finstrain
subdata1$finstrain_resid <- residuals(lm(subdata1$finstrain ~ subdata1$age + subdata1$sex, na.action="na.exclude"))
psych::describe(subdata1$finstrain)
psych::describe(subdata1$finstrain_resid)

# Check distributions
hist(subdata1$finstrain)
hist(subdata1$finstrain_resid)

subdata1$finstrain_tr <- ((subdata1$finstrain_resid+2)^(1/3))*5 #cuberoot transformation
psych::describe(subdata1$finstrain_tr)
hist(subdata1$finstrain_tr)

## totdiet
subdata1$totdiet_resid <- (residuals(lm(subdata1$totdiet ~ subdata1$age + subdata1$sex, na.action="na.exclude")))
psych::describe(subdata1$totdiet)
psych::describe(subdata1$totdiet_resid)

# Check distributions
hist(subdata1$totdiet)
hist(subdata1$totdiet_resid)

subdata1$totdiet_tr <- ((subdata1$totdiet_resid + 25)^(2))/200  #square transformation
psych::describe(subdata1$totdiet_tr)
hist(subdata1$totdiet_tr)

## unhealthdiet
subdata1$unhealthdiet_resid <- (residuals(lm(subdata1$unhealthdiet ~ subdata1$age + subdata1$sex, na.action="na.exclude")))
psych::describe(subdata1$unhealthdiet)
psych::describe(subdata1$unhealthdiet_resid)

# Check distributions
hist(subdata1$unhealthdiet)
hist(subdata1$unhealthdiet_resid)

subdata1$unhealthdiet_tr <- ((subdata1$unhealthdiet_resid + 17)^(3))/1500 #cube transformation
psych::describe(subdata1$unhealthdiet_tr)
hist(subdata1$unhealthdiet_tr)

## healthdiet
subdata1$healthdiet_resid <- (residuals(lm(subdata1$healthdiet ~ subdata1$age + subdata1$sex, na.action="na.exclude")))
psych::describe(subdata1$healthdiet)
psych::describe(subdata1$healthdiet_resid)

# Check distributions
hist(subdata1$healthdiet)
hist(subdata1$healthdiet_resid)

subdata1$healthdiet_tr <- (subdata1$healthdiet_resid + 19)/3 #transformation
psych::describe(subdata1$healthdiet_tr)
hist(subdata1$healthdiet_tr)

## Dep
subdata1$Dep_resid <- (residuals(lm(subdata1$Dep ~ subdata1$age + subdata1$sex, na.action="na.exclude")))
psych::describe(subdata1$Dep)
psych::describe(subdata1$Dep_resid)

# Check distributions
hist(subdata1$Dep)
hist(subdata1$Dep_resid)

subdata1$Dep_tr <- ((log(subdata1$Dep_resid + 8)/3)+.5)*5 #cube transformation
psych::describe(subdata1$Dep_tr)
hist(subdata1$Dep_tr)

## Anx
subdata1$Anx_resid <- (residuals(lm(subdata1$Anx ~ subdata1$age + subdata1$sex, na.action="na.exclude")))
psych::describe(subdata1$Anx)
psych::describe(subdata1$Anx_resid)

# Check distributions
hist(subdata1$Anx)
hist(subdata1$Anx_resid)

subdata1$Anx_tr <- log(subdata1$Anx_resid + 2.5)+2 #log transformation
psych::describe(subdata1$Anx_tr)
hist(subdata1$Anx_tr)

summary(subdata1)
psych::describe(subdata1)

# Select only variables of interest and save them in a smaller file

subdata2 	<- subdata1[,c("FamID","TwNo","zyg","sexzyg","age","sex","edu_resid","elect_resid","elect_trans_resid",
				"finstrain_tr","totdiet_tr","unhealthdiet_tr","healthdiet_tr","Dep_tr","Anx_tr")] 

# Rename variables in subdata1 to make it easier. Note CMR stands for Cardio-Metabolic Risk & AM for Anthropometric (generic terms so that you don't 
# need to change the rest of the script if you select different variables 

names(subdata2) <- c("FamID","TwNo","zyg","sexzyg","age","sex","edu","electct","electtranct",
				"finstrain","totdiet","unhealthdiet","healthdiet","Dep","Anx")

summary(subdata2)
psych::describe(subdata2)

# --------------------------------------------------------------------------------------------------------------------------

# Reshape data from long format to wide format for twin analysis:
subdata3 	<- reshape(subdata2, idvar = c("FamID","sexzyg"), timevar = "TwNo", direction = "wide")

# Check data
str(subdata3)
summary(subdata3)
psych::describe(subdata3)
dim(subdata3)

# Rename variables to get rid of the dot '.' in the names
colnames(subdata3) <- c('FamID','sexzyg','zyg3','age3','sex3','edu3','electct3','electtranct3','finstrain3','totdiet3','unhealthdiet3','healthdiet3','Dep3','Anx3',
						     'zyg1','age1','sex1','edu1','electct1','electtranct1','finstrain1','totdiet1','unhealthdiet1','healthdiet1','Dep1','Anx1', 
						     'zyg2','age2','sex2','edu2','electct2','electtranct2','finstrain2','totdiet2','unhealthdiet2','healthdiet2','Dep2','Anx2') 

psych::describe(subdata3)
dim(subdata3)

# -------------------------------------------------------------------------------------------------------------
#Tabulating variable missingness
# -------------------------------------------------------------------------------------------------------------

##Dep
selVars	<- c('Dep1', 'Dep2')
(Depmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$Dep1) & is.na(subdata3$Dep2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$Dep2) & is.na(subdata3$Dep1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(Depmis))
psych::describe(subdata3$Dep3)

##Anx
selVars	<- c('Anx1', 'Anx2')
(Anxmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$Anx1) & is.na(subdata3$Anx2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$Anx2) & is.na(subdata3$Anx1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(Anxmis))
psych::describe(subdata3$Anx3)

##healthdiet
selVars	<- c('healthdiet1', 'healthdiet2')
(healthdietmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$healthdiet1) & is.na(subdata3$healthdiet2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$healthdiet2) & is.na(subdata3$healthdiet1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(healthdietmis))
psych::describe(subdata3$healthdiet3)

##unhealthdiet
selVars	<- c('unhealthdiet1', 'unhealthdiet2')
(unhealthdietmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$unhealthdiet1) & is.na(subdata3$unhealthdiet2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$unhealthdiet2) & is.na(subdata3$unhealthdiet1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(unhealthdietmis))
psych::describe(subdata3$unhealthdiet3)

##totdiet
selVars	<- c('totdiet1', 'totdiet2')
(totdietmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$totdiet1) & is.na(subdata3$totdiet2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$totdiet2) & is.na(subdata3$totdiet1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(totdietmis))
psych::describe(subdata3$totdiet3)

##finstrain
selVars	<- c('finstrain1', 'finstrain2')
(finstrainmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$finstrain1) & is.na(subdata3$finstrain2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$finstrain2) & is.na(subdata3$finstrain1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(finstrainmis))
psych::describe(subdata3$finstrain3)

##electtranct
selVars	<- c('electtranct1', 'electtranct2')
(electtranctmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$electtranct1) & is.na(subdata3$electtranct2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$electtranct2) & is.na(subdata3$electtranct1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(electtranctmis))
psych::describe(subdata3$electtranct3)

##electct
selVars	<- c('electct1', 'electct2')
(electctmis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$electct1) & is.na(subdata3$electct2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$electct2) & is.na(subdata3$electct1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(electctmis))
psych::describe(subdata3$electct3)

##edu
selVars	<- c('edu1', 'edu2')
(edumis	<-matrix(
	c(
	(totalpairs		<-table(subdata3$zyg1)),
	(complete_pairs	<-table(subdata3$zyg1[complete.cases(subdata3[,selVars])])),
	(incomplete1	<-table(subdata3$zyg1[!is.na(subdata3$edu1) & is.na(subdata3$edu2)])),
	(incomplete2	<-table(subdata3$zyg2[!is.na(subdata3$edu2) & is.na(subdata3$edu1)])),
	(incomplete		<-incomplete1 + incomplete2),
	(totalsample	<-(2*complete_pairs) + incomplete)
	),
	ncol=2,
	byrow=TRUE,
	dimnames=list(c("total pairs", "comp pairs","incomp pairs1","incomp pairs2","incomp sing","n (total)"), c("mz","dz"))
	))
(rowtotal		<-rowSums(edumis))
psych::describe(subdata3$edu3)

# -------------------------------------------------------------------------------------------------------------
# 3. Prep Data files for twin analyses
# -------------------------------------------------------------------------------------------------------------

#selVars	<- c('edu1','electct1','electtranct1','finstrain1','totdiet1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electct2','electtranct2','finstrain2','totdiet2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electct3','electtranct3','finstrain3','totdiet3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON


# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)


#**************************************************************************************************************
# (5) Check means and correlations and covariances of the data files
# These will not be used, they are just meant to get a feel for the data and to help decide on starting values
#************************************************************************************************************
#*******************************************************************
#  Twin Means, Covariances & Correlations
#-------------------------------------------------------------------
# Means
(colMeans(mzmData, na.rm=T))
(colMeans(dzmData, na.rm=T))
(colMeans(mzfData, na.rm=T))
(colMeans(dzfData, na.rm=T))
(colMeans(dzoData, na.rm=T))

# Variance-Covariance Matrices
(cov(mzmData, use="complete"))
(cov(dzmData, use="complete"))
(cov(mzfData, use="complete"))
(cov(dzfData, use="complete"))
(cov(dzoData, use="complete"))

# Correlations
(cor(mzmData, use="complete"))
(cor(dzmData, use="complete"))
(cor(mzfData, use="complete"))
(cor(dzfData, use="complete"))
(cor(dzoData, use="complete"))


#------------------------------------------------------------------------------------------
# (II-1) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
#	 for males and females separately using same-sex twin pair & DZ-OS PAIRS
# Education
#------------------------------------------------------------------------------------------
nv		<- 1				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor		<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars		<- c('edu1', 'edu2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('edu3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables

Stpathf	<-c(.4) 		
(Stmean	<-4.6)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.1, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE1Model	<-mxModel('HetACE1', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE1Fit    <- mxRun(HetACE1Model, intervals=F)
(HetACE1Summ  <- summary(HetACE1Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE1Fit)
mxEval (MZM.cm2, HetACE1Fit)
mxEval (MZM.em2, HetACE1Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE1Fit)
mxEval (MZF.cf2, HetACE1Fit)
mxEval (MZF.ef2, HetACE1Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE1cModel	<-mxModel(HetACE1Fit, name="HetACE1c")
HetACE1cModel	<-omxSetParameters(HetACE1cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE1cFit		<-mxTryHard(HetACE1cModel, intervals=F)
(HetACE1cSumm	<-summary(HetACE1cFit))

mxEval(MZM.Vm, HetACE1cFit)
mxEval(MZM.hm2, HetACE1cFit)
mxEval(MZM.cm2, HetACE1cFit)
mxEval(MZM.em2, HetACE1cFit)

mxEval(MZF.Vf, HetACE1cFit)
mxEval(MZF.hf2, HetACE1cFit)
mxEval(MZF.cf2, HetACE1cFit)
mxEval(MZF.ef2, HetACE1cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE1aModel	<-mxModel(HetACE1Fit, name="HetACE1a")
HetACE1aModel	<-omxSetParameters(HetACE1aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE1aFit		<-mxTryHard(HetACE1aModel, intervals=F)
(HetACE1aSumm	<-summary(HetACE1aFit, verbose = F))

mxEval(MZM.Vm, HetACE1aFit)
mxEval(MZM.hm2, HetACE1aFit)
mxEval(MZM.cm2, HetACE1aFit)
mxEval(MZM.em2, HetACE1aFit)

mxEval(MZF.Vf, HetACE1aFit)
mxEval(MZF.hf2, HetACE1aFit)
mxEval(MZF.cf2, HetACE1aFit)
mxEval(MZF.ef2, HetACE1aFit)

#---------------------------------------------------------------------
# Model II-4 - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE1Model	<-mxModel(HetACE1Fit, name="HomACE1")
HomACE1Model	<-omxSetParameters(HomACE1Model, 	labels=	c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE1Model	<-omxSetParameters(HomACE1Model, 	labels=	c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE1Model	<-omxSetParameters(HomACE1Model, 	labels=	c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE1Model	<-omxSetParameters(HomACE1Model, 	labels=	c("mf1"), free=T, values=c(4.6), newlabels=c("mm1"))
HomACE1Model	<-omxAssignFirstParameters(HomACE1Model) 
HomACE1Fit		<-mxTryHard(HomACE1Model, intervals=F)
(HomACE1Summ	<-summary(HomACE1Fit))


#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE1Fit)
mxEval(MZM.hm2, HomACE1Fit)
mxEval(MZM.cm2, HomACE1Fit)
mxEval(MZM.em2, HomACE1Fit)

mxEval(MZF.Vf, HomACE1Fit)
mxEval(MZF.hf2, HomACE1Fit)
mxEval(MZF.cf2, HomACE1Fit)
mxEval(MZF.ef2, HomACE1Fit)


#--------------------------------------------------------------------------------------------------
# Model II-5a - Scalar (variance inequality) Model - Education
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models
# --------------------------------------------------------------------------------------------------

nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('edu1', 'edu2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('edu3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

(LabAf	<- paste("apf",1:nv,sep=""))		
(LabCf	<- paste("cpf",1:nv,sep=""))		
(LabEf	<- paste("epf",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE1Model	<-mxModel("ScACE1", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model II-5
ScACE1Fit	<-mxTryHard(ScACE1Model, intervals=F)
(ScACE1Sum 	<-summary(ScACE1Fit))

mxEval (MZM.h2m, ScACE1Fit)
mxEval (MZM.c2m, ScACE1Fit)
mxEval (MZM.e2m, ScACE1Fit)
mxEval (MZM.am, ScACE1Fit)
mxEval (MZM.cm, ScACE1Fit)
mxEval (MZM.em, ScACE1Fit)

mxEval (MZF.h2f, ScACE1Fit)
mxEval (MZF.c2f, ScACE1Fit)
mxEval (MZF.e2f, ScACE1Fit)
mxEval (MZF.af, ScACE1Fit)
mxEval (MZF.cf, ScACE1Fit)
mxEval (MZF.ef, ScACE1Fit)
mxEval (MZF.ScF, ScACE1Fit)

#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE1cFit, HetACE1Fit)
mxCompare(HetACE1aFit, HetACE1Fit)
mxCompare(HetACE1Fit, HomACE1Fit)
mxCompare(HetACE1Fit, ScACE1Fit)
##
##

#------------------------------------------------------------------------------------------
# (II-1) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
#	 for males and females separately using same-sex twin pair & DZ-OS PAIRS
# electtranct
#------------------------------------------------------------------------------------------
nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('electtranct1', 'electtranct2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('electtranct3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E latent Factor Correlation Matrices (for the correlation between variables)
#(LabsRAM		<- paste("Ram",1,sep=""))		# Males (Twins & singletons)
#(LabsRCM		<- paste("Rcm",1,sep=""))		# Males (Twins & singletons)
#(LabsREM		<- paste("Rem",1,sep=""))		# Males (Twins & singletons)
#(LabsRAF		<- paste("Raf",1,sep=""))	# Females (Twins & singletons)
#(LabsRCF		<- paste("Rcf",1,sep=""))		# Females (Twins & singletons)
#(LabsREF		<- paste("Ref",1,sep=""))		# Females (Twins & singletons)

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables
#Stcorm	<-c(.4)

Stpathf	<-c(.4) 		
#Stcorf	<-c(.4)
(Stmean	<-3.2)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.1, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE2Model	<-mxModel('HetACE2', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE2Fit    <- mxRun(HetACE2Model, intervals=F)
(HetACE2Summ  <- summary(HetACE2Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE2Fit)
mxEval (MZM.cm2, HetACE2Fit)
mxEval (MZM.em2, HetACE2Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE2Fit)
mxEval (MZF.cf2, HetACE2Fit)
mxEval (MZF.ef2, HetACE2Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE2cModel	<-mxModel(HetACE2Fit, name="HetACE2c")
HetACE2cModel	<-omxSetParameters(HetACE2cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE2cFit		<-mxTryHard(HetACE2cModel, intervals=F)
(HetACE2cSumm	<-summary(HetACE2cFit))

mxEval(MZM.Vm, HetACE2cFit)
mxEval(MZM.hm2, HetACE2cFit)
mxEval(MZM.cm2, HetACE2cFit)
mxEval(MZM.em2, HetACE2cFit)

mxEval(MZF.Vf, HetACE2cFit)
mxEval(MZF.hf2, HetACE2cFit)
mxEval(MZF.cf2, HetACE2cFit)
mxEval(MZF.ef2, HetACE2cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE2aModel	<-mxModel(HetACE2Fit, name="HetACE2a")
HetACE2aModel	<-omxSetParameters(HetACE2aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE2aFit		<-mxTryHard(HetACE2aModel, intervals=F)
(HetACE2aSumm	<-summary(HetACE2aFit, verbose = F))

mxEval(MZM.Vm, HetACE2aFit)
mxEval(MZM.hm2, HetACE2aFit)
mxEval(MZM.cm2, HetACE2aFit)
mxEval(MZM.em2, HetACE2aFit)

mxEval(MZF.Vf, HetACE2aFit)
mxEval(MZF.hf2, HetACE2aFit)
mxEval(MZF.cf2, HetACE2aFit)
mxEval(MZF.ef2, HetACE2aFit)

#---------------------------------------------------------------------
# Model II-4b - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE2Model	<-mxModel(HetACE2Fit, name="HomACE2")
HomACE2Model	<-omxSetParameters(HomACE2Model, 	labels=	c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE2Model	<-omxSetParameters(HomACE2Model, 	labels=	c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE2Model	<-omxSetParameters(HomACE2Model, 	labels=	c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE2Model	<-omxSetParameters(HomACE2Model, 	labels=	c("mf1"), free=T, values=c(4.6), newlabels=c("mm1"))
HomACE2Model	<-omxAssignFirstParameters(HomACE2Model) 
HomACE2Fit		<-mxTryHard(HomACE2Model, intervals=F)
(HomACE2Summ	<-summary(HomACE2Fit))


#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE2Fit)
mxEval(MZM.hm2, HomACE2Fit)
mxEval(MZM.cm2, HomACE2Fit)
mxEval(MZM.em2, HomACE2Fit)

mxEval(MZF.Vf, HomACE2Fit)
mxEval(MZF.hf2, HomACE2Fit)
mxEval(MZF.cf2, HomACE2Fit)
mxEval(MZF.ef2, HomACE2Fit)


#--------------------------------------------------------------------------------------------------
# Model II-5a - Scalar (variance inequality) Model - Education
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models
# --------------------------------------------------------------------------------------------------

nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('electtranct1', 'electtranct2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('electtranct3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE2Model	<-mxModel("ScACE2", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model III-5
ScACE2Fit	<-mxTryHard(ScACE2Model, intervals=F)
(ScACE2Sum 	<-summary(ScACE2Fit))

mxEval (MZM.h2m, ScACE2Fit)
mxEval (MZM.c2m, ScACE2Fit)
mxEval (MZM.e2m, ScACE2Fit)
mxEval (MZM.am, ScACE2Fit)
mxEval (MZM.cm, ScACE2Fit)
mxEval (MZM.em, ScACE2Fit)

mxEval (MZF.h2f, ScACE2Fit)
mxEval (MZF.c2f, ScACE2Fit)
mxEval (MZF.e2f, ScACE2Fit)
mxEval (MZF.af, ScACE2Fit)
mxEval (MZF.cf, ScACE2Fit)
mxEval (MZF.ef, ScACE2Fit)
mxEval (MZF.ScF, ScACE2Fit)


#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE2cFit, HetACE2Fit)
mxCompare(HetACE2aFit, HetACE2Fit)
mxCompare(HetACE2Fit, HomACE2Fit)
mxCompare(HetACE2Fit, ScACE2Fit)

##
##
##

#------------------------------------------------------------------------------------------
# (II-1c) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
#	 for males and females separately using same-sex twin pair & DZ-OS PAIRS
# Finstrain
#------------------------------------------------------------------------------------------
nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('finstrain1', 'finstrain2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('finstrain3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E latent Factor Correlation Matrices (for the correlation between variables)
#(LabsRAM		<- paste("Ram",1,sep=""))		# Males (Twins & singletons)
#(LabsRCM		<- paste("Rcm",1,sep=""))		# Males (Twins & singletons)
#(LabsREM		<- paste("Rem",1,sep=""))		# Males (Twins & singletons)
#(LabsRAF		<- paste("Raf",1,sep=""))	# Females (Twins & singletons)
#(LabsRCF		<- paste("Rcf",1,sep=""))		# Females (Twins & singletons)
#(LabsREF		<- paste("Ref",1,sep=""))		# Females (Twins & singletons)

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables
#Stcorm	<-c(.4)

Stpathf	<-c(.4) 		
#Stcorf	<-c(.4)
(Stmean	<-3.2)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.1, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE3Model	<-mxModel('HetACE3', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE3Fit    <- mxRun(HetACE3Model, intervals=F)
(HetACE3Summ  <- summary(HetACE3Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE3Fit)
mxEval (MZM.cm2, HetACE3Fit)
mxEval (MZM.em2, HetACE3Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE3Fit)
mxEval (MZF.cf2, HetACE3Fit)
mxEval (MZF.ef2, HetACE3Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE3cModel	<-mxModel(HetACE3Fit, name="HetACE3c")
HetACE3cModel	<-omxSetParameters(HetACE3cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE3cFit		<-mxTryHard(HetACE3cModel, intervals=F)
(HetACE3cSumm	<-summary(HetACE3cFit))

mxEval(MZM.Vm, HetACE3cFit)
mxEval(MZM.hm2, HetACE3cFit)
mxEval(MZM.cm2, HetACE3cFit)
mxEval(MZM.em2, HetACE3cFit)

mxEval(MZF.Vf, HetACE3cFit)
mxEval(MZF.hf2, HetACE3cFit)
mxEval(MZF.cf2, HetACE3cFit)
mxEval(MZF.ef2, HetACE3cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE3aModel	<-mxModel(HetACE3Fit, name="HetACE3a")
HetACE3aModel	<-omxSetParameters(HetACE3aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE3aFit		<-mxTryHard(HetACE3aModel, intervals=F)
(HetACE3aSumm	<-summary(HetACE3aFit, verbose = F))

mxEval(MZM.Vm, HetACE3aFit)
mxEval(MZM.hm2, HetACE3aFit)
mxEval(MZM.cm2, HetACE3aFit)
mxEval(MZM.em2, HetACE3aFit)

mxEval(MZF.Vf, HetACE3aFit)
mxEval(MZF.hf2, HetACE3aFit)
mxEval(MZF.cf2, HetACE3aFit)
mxEval(MZF.ef2, HetACE3aFit)

#---------------------------------------------------------------------
# Model II-4b - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE3Model	<-mxModel(HetACE3Fit, name="HomACE3")
HomACE3Model	<-omxSetParameters(HomACE3Model, 	labels=	c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE3Model	<-omxSetParameters(HomACE3Model, 	labels=	c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE3Model	<-omxSetParameters(HomACE3Model, 	labels=	c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE3Model	<-omxSetParameters(HomACE3Model, 	labels=	c("mf1"), free=T, values=c(4.6), newlabels=c("mm1"))
HomACE3Model	<-omxAssignFirstParameters(HomACE3Model) 
HomACE3Fit		<-mxTryHard(HomACE3Model, intervals=F)
(HomACE3Summ	<-summary(HomACE3Fit))


#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE3Fit)
mxEval(MZM.hm2, HomACE3Fit)
mxEval(MZM.cm2, HomACE3Fit)
mxEval(MZM.em2, HomACE3Fit)

mxEval(MZF.Vf, HomACE3Fit)
mxEval(MZF.hf2, HomACE3Fit)
mxEval(MZF.cf2, HomACE3Fit)
mxEval(MZF.ef2, HomACE3Fit)


#--------------------------------------------------------------------------------------------------
# Model II-5c - Scalar (variance inequality) Model - Financial strain
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models
# --------------------------------------------------------------------------------------------------

nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('finstrain1', 'finstrain2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('finstrain3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

(LabAf	<- paste("apm",1:nv,sep=""))		
(LabCf	<- paste("cpm",1:nv,sep=""))		
(LabEf	<- paste("epm",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEf, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE3Model	<-mxModel("ScACE3", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model III-5
ScACE3Fit	<-mxTryHard(ScACE3Model, intervals=F)
(ScACE3Sum 	<-summary(ScACE3Fit))

mxEval (MZM.h2m, ScACE3Fit)
mxEval (MZM.c2m, ScACE3Fit)
mxEval (MZM.e2m, ScACE3Fit)
mxEval (MZM.am, ScACE3Fit)
mxEval (MZM.cm, ScACE3Fit)
mxEval (MZM.em, ScACE3Fit)

mxEval (MZF.h2f, ScACE3Fit)
mxEval (MZF.c2f, ScACE3Fit)
mxEval (MZF.e2f, ScACE3Fit)
mxEval (MZF.af, ScACE3Fit)
mxEval (MZF.cf, ScACE3Fit)
mxEval (MZF.ef, ScACE3Fit)
mxEval (MZF.ScF, ScACE3Fit)


#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE3cFit, HetACE3Fit)
mxCompare(HetACE3aFit, HetACE3Fit)
mxCompare(HetACE3Fit, HomACE3Fit)
mxCompare(HetACE3Fit, ScACE3Fit)

##
##

#------------------------------------------------------------------------------------------
# (II-1d) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
#	 for males and females separately using same-sex twin pair & DZ-OS PAIRS
# Unhealthdiet
#------------------------------------------------------------------------------------------
nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('unhealthdiet1', 'unhealthdiet2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('unhealthdiet3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E latent Factor Correlation Matrices (for the correlation between variables)
#(LabsRAM		<- paste("Ram",1,sep=""))		# Males (Twins & singletons)
#(LabsRCM		<- paste("Rcm",1,sep=""))		# Males (Twins & singletons)
#(LabsREM		<- paste("Rem",1,sep=""))		# Males (Twins & singletons)
#(LabsRAF		<- paste("Raf",1,sep=""))	# Females (Twins & singletons)
#(LabsRCF		<- paste("Rcf",1,sep=""))		# Females (Twins & singletons)
#(LabsREF		<- paste("Ref",1,sep=""))		# Females (Twins & singletons)

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables
#Stcorm	<-c(.4)

Stpathf	<-c(.4) 		
#Stcorf	<-c(.4)
(Stmean	<-3.2)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE4Model	<-mxModel('HetACE4', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE4Fit    <- mxRun(HetACE4Model, intervals=F)
(HetACE4Summ  <- summary(HetACE4Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE4Fit)
mxEval (MZM.cm2, HetACE4Fit)
mxEval (MZM.em2, HetACE4Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE4Fit)
mxEval (MZF.cf2, HetACE4Fit)
mxEval (MZF.ef2, HetACE4Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE4cModel	<-mxModel(HetACE4Fit, name="HetACE4c")
HetACE4cModel	<-omxSetParameters(HetACE4cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE4cFit		<-mxTryHard(HetACE4cModel, intervals=F)
(HetACE4cSumm	<-summary(HetACE4cFit))

mxEval(MZM.Vm, HetACE4cFit)
mxEval(MZM.hm2, HetACE4cFit)
mxEval(MZM.cm2, HetACE4cFit)
mxEval(MZM.em2, HetACE4cFit)

mxEval(MZF.Vf, HetACE4cFit)
mxEval(MZF.hf2, HetACE4cFit)
mxEval(MZF.cf2, HetACE4cFit)
mxEval(MZF.ef2, HetACE4cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE4aModel	<-mxModel(HetACE4Fit, name="HetACE4a")
HetACE4aModel	<-omxSetParameters(HetACE4aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE4aFit		<-mxTryHard(HetACE4aModel, intervals=F)
(HetACE4aSumm	<-summary(HetACE4aFit, verbose = F))

mxEval(MZM.Vm, HetACE4aFit)
mxEval(MZM.hm2, HetACE4aFit)
mxEval(MZM.cm2, HetACE4aFit)
mxEval(MZM.em2, HetACE4aFit)

mxEval(MZF.Vf, HetACE4aFit)
mxEval(MZF.hf2, HetACE4aFit)
mxEval(MZF.cf2, HetACE4aFit)
mxEval(MZF.ef2, HetACE4aFit)

#---------------------------------------------------------------------
# Model II-4b - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE4Model	<-mxModel(HetACE4Fit, name="HomACE4")
HomACE4Model	<-omxSetParameters(HomACE4Model, 	labels=	c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE4Model	<-omxSetParameters(HomACE4Model, 	labels=	c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE4Model	<-omxSetParameters(HomACE4Model, 	labels=	c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE4Model	<-omxSetParameters(HomACE4Model, 	labels=	c("mf1"), free=T, values=c(4.0), newlabels=c("mm1"))
HomACE4Model	<-omxAssignFirstParameters(HomACE4Model) 
HomACE4Fit		<-mxTryHard(HomACE4Model, intervals=F)
(HomACE4Summ	<-summary(HomACE4Fit))

#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE4Fit)
mxEval(MZM.hm2, HomACE4Fit)
mxEval(MZM.cm2, HomACE4Fit)
mxEval(MZM.em2, HomACE4Fit)

mxEval(MZF.Vf, HomACE4Fit)
mxEval(MZF.hf2, HomACE4Fit)
mxEval(MZF.cf2, HomACE4Fit)
mxEval(MZF.ef2, HomACE4Fit)


#--------------------------------------------------------------------------------------------------
# Model IV-5c - Scalar (variance inequality) Model - Financial strain
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models
# --------------------------------------------------------------------------------------------------

nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('unhealthdiet1', 'unhealthdiet2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('unhealthdiet3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

(LabAf	<- paste("apm",1:nv,sep=""))		
(LabCf	<- paste("cpm",1:nv,sep=""))		
(LabEf	<- paste("epm",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEf, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE4Model	<-mxModel("ScACE4", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model II-5
ScACE4Fit	<-mxTryHard(ScACE4Model, intervals=F)
(ScACE4Sum 	<-summary(ScACE1Fit))

mxEval (MZM.h2m, ScACE4Fit)
mxEval (MZM.c2m, ScACE4Fit)
mxEval (MZM.e2m, ScACE4Fit)
mxEval (MZM.am, ScACE4Fit)
mxEval (MZM.cm, ScACE4Fit)
mxEval (MZM.em, ScACE4Fit)

mxEval (MZF.h2f, ScACE4Fit)
mxEval (MZF.c2f, ScACE4Fit)
mxEval (MZF.e2f, ScACE4Fit)
mxEval (MZF.af, ScACE4Fit)
mxEval (MZF.cf, ScACE4Fit)
mxEval (MZF.ef, ScACE4Fit)
mxEval (MZF.ScF, ScACE4Fit)

#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE4cFit, HetACE4Fit)
mxCompare(HetACE4aFit, HetACE4Fit)
mxCompare(HetACE4Fit, HomACE4Fit)
mxCompare(HetACE4Fit, ScACE4Fit)

##
##

#------------------------------------------------------------------------------------------
# (II-1e) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
# for males and females separately using same-sex twin pair & DZ-OS PAIRS
# Health diet
#------------------------------------------------------------------------------------------
nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('healthdiet1', 'healthdiet2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('healthdiet3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E latent Factor Correlation Matrices (for the correlation between variables)
#(LabsRAM		<- paste("Ram",1,sep=""))		# Males (Twins & singletons)
#(LabsRCM		<- paste("Rcm",1,sep=""))		# Males (Twins & singletons)
#(LabsREM		<- paste("Rem",1,sep=""))		# Males (Twins & singletons)
#(LabsRAF		<- paste("Raf",1,sep=""))	# Females (Twins & singletons)
#(LabsRCF		<- paste("Rcf",1,sep=""))		# Females (Twins & singletons)
#(LabsREF		<- paste("Ref",1,sep=""))		# Females (Twins & singletons)

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables
#Stcorm	<-c(.4)

Stpathf	<-c(.4) 		
#Stcorf	<-c(.4)
(Stmean	<-3.2)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE5Model	<-mxModel('HetACE5', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE5Fit    <- mxRun(HetACE5Model, intervals=F)
(HetACE5Summ  <- summary(HetACE5Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE5Fit)
mxEval (MZM.cm2, HetACE5Fit)
mxEval (MZM.em2, HetACE5Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE5Fit)
mxEval (MZF.cf2, HetACE5Fit)
mxEval (MZF.ef2, HetACE5Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE5cModel	<-mxModel(HetACE5Fit, name="HetACE5c")
HetACE5cModel	<-omxSetParameters(HetACE5cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE5cFit		<-mxTryHard(HetACE5cModel, intervals=F)
(HetACE5cSumm	<-summary(HetACE5cFit))

mxEval(MZM.Vm, HetACE5cFit)
mxEval(MZM.hm2, HetACE5cFit)
mxEval(MZM.cm2, HetACE5cFit)
mxEval(MZM.em2, HetACE5cFit)

mxEval(MZF.Vf, HetACE5cFit)
mxEval(MZF.hf2, HetACE5cFit)
mxEval(MZF.cf2, HetACE5cFit)
mxEval(MZF.ef2, HetACE5cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE5aModel	<-mxModel(HetACE5Fit, name="HetACE5a")
HetACE5aModel	<-omxSetParameters(HetACE5aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE5aFit		<-mxTryHard(HetACE5aModel, intervals=F)
(HetACE5aSumm	<-summary(HetACE5aFit, verbose = F))

mxEval(MZM.Vm, HetACE5aFit)
mxEval(MZM.hm2, HetACE5aFit)
mxEval(MZM.cm2, HetACE5aFit)
mxEval(MZM.em2, HetACE5aFit)

mxEval(MZF.Vf, HetACE5aFit)
mxEval(MZF.hf2, HetACE5aFit)
mxEval(MZF.cf2, HetACE5aFit)
mxEval(MZF.ef2, HetACE5aFit)

#---------------------------------------------------------------------
# Model II-4b - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE5Model	<-mxModel(HetACE5Fit, name="HomACE5")
HomACE5Model	<-omxSetParameters(HomACE5Model, 	labels=c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE5Model	<-omxSetParameters(HomACE5Model, 	labels=c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE5Model	<-omxSetParameters(HomACE5Model, 	labels=c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE5Model	<-omxSetParameters(HomACE5Model, 	labels=c("mf1"), free=T, values=c(4.0), newlabels=c("mm1"))
HomACE5Model	<-omxAssignFirstParameters(HomACE5Model) 
HomACE5Fit		<-mxTryHard(HomACE5Model, intervals=F)
(HomACE5Summ	<-summary(HomACE5Fit))

#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE5Fit)
mxEval(MZM.hm2, HomACE5Fit)
mxEval(MZM.cm2, HomACE5Fit)
mxEval(MZM.em2, HomACE5Fit)

mxEval(MZF.Vf, HomACE5Fit)
mxEval(MZF.hf2, HomACE5Fit)
mxEval(MZF.cf2, HomACE5Fit)
mxEval(MZF.ef2, HomACE5Fit)


#--------------------------------------------------------------------------------------------------
# Model II-5e - Scalar (variance inequality) Model - Healthy diet
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models
# --------------------------------------------------------------------------------------------------

nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('healthdiet1', 'healthdiet2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('healthdiet3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

(LabAf	<- paste("apm",1:nv,sep=""))		
(LabCf	<- paste("cpm",1:nv,sep=""))		
(LabEf	<- paste("epm",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEf, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE5Model	<-mxModel("ScACE5", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model V-5
ScACE5Fit	<-mxTryHard(ScACE5Model, intervals=F)
(ScACE5Sum 	<-summary(ScACE5Fit))

mxEval (MZM.h2m, ScACE5Fit)
mxEval (MZM.c2m, ScACE5Fit)
mxEval (MZM.e2m, ScACE5Fit)
mxEval (MZM.am, ScACE5Fit)
mxEval (MZM.cm, ScACE1Fit)
mxEval (MZM.em, ScACE5Fit)

mxEval (MZF.h2f, ScACE5Fit)
mxEval (MZF.c2f, ScACE5Fit)
mxEval (MZF.e2f, ScACE5Fit)
mxEval (MZF.af, ScACE5Fit)
mxEval (MZF.cf, ScACE5Fit)
mxEval (MZF.ef, ScACE5Fit)
mxEval (MZF.ScF, ScACE5Fit)

#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE5cFit, HetACE5Fit)
mxCompare(HetACE5aFit, HetACE5Fit)
mxCompare(HetACE5Fit, HomACE5Fit)
mxCompare(HetACE5Fit, ScACE5Fit)

##
##

#------------------------------------------------------------------------------------------
# (II-1f) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
# for males and females separately using same-sex twin pair & DZ-OS PAIRS
# Depressive symptoms
#------------------------------------------------------------------------------------------
nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('Dep1', 'Dep2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('Dep3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E latent Factor Correlation Matrices (for the correlation between variables)
#(LabsRAM		<- paste("Ram",1,sep=""))		# Males (Twins & singletons)
#(LabsRCM		<- paste("Rcm",1,sep=""))		# Males (Twins & singletons)
#(LabsREM		<- paste("Rem",1,sep=""))		# Males (Twins & singletons)
#(LabsRAF		<- paste("Raf",1,sep=""))	# Females (Twins & singletons)
#(LabsRCF		<- paste("Rcf",1,sep=""))		# Females (Twins & singletons)
#(LabsREF		<- paste("Ref",1,sep=""))		# Females (Twins & singletons)

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables
#Stcorm	<-c(.4)

Stpathf	<-c(.4) 		
#Stcorf	<-c(.4)
(Stmean	<-3.2)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE6Model	<-mxModel('HetACE6', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE6Fit    <- mxRun(HetACE6Model, intervals=F)
(HetACE6Summ  <- summary(HetACE6Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE6Fit)
mxEval (MZM.cm2, HetACE6Fit)
mxEval (MZM.em2, HetACE6Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE6Fit)
mxEval (MZF.cf2, HetACE6Fit)
mxEval (MZF.ef2, HetACE6Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE6cModel	<-mxModel(HetACE6Fit, name="HetACE6c")
HetACE6cModel	<-omxSetParameters(HetACE6cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE6cFit		<-mxTryHard(HetACE6cModel, intervals=F)
(HetACE6cSumm	<-summary(HetACE6cFit))

mxEval(MZM.Vm, HetACE6cFit)
mxEval(MZM.hm2, HetACE6cFit)
mxEval(MZM.cm2, HetACE6cFit)
mxEval(MZM.em2, HetACE6cFit)

mxEval(MZF.Vf, HetACE6cFit)
mxEval(MZF.hf2, HetACE6cFit)
mxEval(MZF.cf2, HetACE6cFit)
mxEval(MZF.ef2, HetACE6cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE6aModel	<-mxModel(HetACE6Fit, name="HetACE6a")
HetACE6aModel	<-omxSetParameters(HetACE6aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE6aFit		<-mxTryHard(HetACE6aModel, intervals=F)
(HetACE6aSumm	<-summary(HetACE6aFit, verbose = F))

mxEval(MZM.Vm, HetACE6aFit)
mxEval(MZM.hm2, HetACE6aFit)
mxEval(MZM.cm2, HetACE6aFit)
mxEval(MZM.em2, HetACE6aFit)

mxEval(MZF.Vf, HetACE6aFit)
mxEval(MZF.hf2, HetACE6aFit)
mxEval(MZF.cf2, HetACE6aFit)
mxEval(MZF.ef2, HetACE6aFit)

#---------------------------------------------------------------------
# Model II-4b - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE6Model	<-mxModel(HetACE6Fit, name="HomACE6")
HomACE6Model	<-omxSetParameters(HomACE6Model, 	labels=c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE6Model	<-omxSetParameters(HomACE6Model, 	labels=c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE6Model	<-omxSetParameters(HomACE6Model, 	labels=c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE6Model	<-omxSetParameters(HomACE6Model, 	labels=c("mf1"), free=T, values=c(4.0), newlabels=c("mm1"))
HomACE6Model	<-omxAssignFirstParameters(HomACE6Model) 
HomACE6Fit		<-mxTryHard(HomACE6Model, intervals=F)
(HomACE6Summ	<-summary(HomACE6Fit))

#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE6Fit)
mxEval(MZM.hm2, HomACE6Fit)
mxEval(MZM.cm2, HomACE6Fit)
mxEval(MZM.em2, HomACE6Fit)

mxEval(MZF.Vf, HomACE6Fit)
mxEval(MZF.hf2, HomACE6Fit)
mxEval(MZF.cf2, HomACE6Fit)
mxEval(MZF.ef2, HomACE6Fit)


#--------------------------------------------------------------------------------------------------
# Model II-5f - Scalar (variance inequality) Model - Depressive symptoms
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models 
# --------------------------------------------------------------------------------------------------

nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('Dep1', 'Dep2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('Dep3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

(LabAf	<- paste("apm",1:nv,sep=""))		
(LabCf	<- paste("cpm",1:nv,sep=""))		
(LabEf	<- paste("epm",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEf, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE6Model	<-mxModel("ScACE6", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model VI-5
ScACE6Fit	<-mxTryHard(ScACE6Model, intervals=F)
(ScACE6Sum 	<-summary(ScACE6Fit))

mxEval (MZM.h2m, ScACE6Fit)
mxEval (MZM.c2m, ScACE6Fit)
mxEval (MZM.e2m, ScACE6Fit)
mxEval (MZM.am, ScACE6Fit)
mxEval (MZM.cm, ScACE6Fit)
mxEval (MZM.em, ScACE6Fit)

mxEval (MZF.h2f, ScACE6Fit)
mxEval (MZF.c2f, ScACE6Fit)
mxEval (MZF.e2f, ScACE6Fit)
mxEval (MZF.af, ScACE6Fit)
mxEval (MZF.cf, ScACE6Fit)
mxEval (MZF.ef, ScACE6Fit)
mxEval (MZF.ScF, ScACE6Fit)


#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE6cFit, HetACE6Fit)
mxCompare(HetACE6aFit, HetACE6Fit)
mxCompare(HetACE6Fit, HomACE6Fit)
mxCompare(HetACE6Fit, ScACE6Fit)

##
##

#------------------------------------------------------------------------------------------
# (II-1f) Specify the Quantitative Heterogeneity ACE Model (univariate), modelling a set of ACE parameters 
# for males and females separately using same-sex twin pair & DZ-OS PAIRS
# Depressive symptoms
#------------------------------------------------------------------------------------------
nv			<- 1				# number of variables for a twin = 1 in Univariate
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower		<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor			<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('Anx1', 'Anx2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('Anx3')		# THE VARIABLES FOR A SINGLETON

#selVars	<- c('edu1','electtranct1','finstrain1','unhealthdiet1','healthdiet1','Dep1','Anx1',
#			'edu2','electtranct2','finstrain2','unhealthdiet2','healthdiet2','Dep2','Anx2')	# THE VARIABLES FOR TWIN PAIRS
#selVarsS 	<- c('edu3','electtranct3','finstrain3','unhealthdiet3','healthdiet3','Dep3','Anx3')				# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData		<- subset(subdata3, sexzyg == 6, selVarsS)
sfData		<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))				# male-male pairs
(LabMFF	<- c(LabMF,LabMF))				# female-female pairs
(LabMMF	<- c(LabMM,LabMF))				# male-female pairs

# Create Labels for the A, C and E latent Factor Correlation Matrices (for the correlation between variables)
#(LabsRAM		<- paste("Ram",1,sep=""))		# Males (Twins & singletons)
#(LabsRCM		<- paste("Rcm",1,sep=""))		# Males (Twins & singletons)
#(LabsREM		<- paste("Rem",1,sep=""))		# Males (Twins & singletons)
#(LabsRAF		<- paste("Raf",1,sep=""))	# Females (Twins & singletons)
#(LabsRCF		<- paste("Rcf",1,sep=""))		# Females (Twins & singletons)
#(LabsREF		<- paste("Ref",1,sep=""))		# Females (Twins & singletons)

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables)
(LabAm	<- paste("amp",1:nv,sep=""))		# Males
(LabCm	<- paste("cmp",1:nv,sep=""))		# Males
(LabEm	<- paste("emp",1:nv,sep=""))		# Males
(LabAf	<- paste("afp",1:nv,sep=""))		# Females
(LabCf	<- paste("cfp",1:nv,sep=""))		# Females
(LabEf	<- paste("efp",1:nv,sep=""))		# Females

# Create objects for Start values
Stpathm	<-c(.4) 			# change here according to the variances of the variables
#Stcorm	<-c(.4)

Stpathf	<-c(.4) 		
#Stcorf	<-c(.4)
(Stmean	<-3.2)

# ---------------------------------------------------------------------------------------------------------------------
# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

# Define matrices a, c, and e to store a, c, and e path coefficients
pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.7, label=LabEm, name="em")
pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.3, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=.6, label=LabEf, name="ef")

# Setting the A & C correlations between males and females from the DZ opposite-sex pairs 
# In this 1st model, rA = .5 (free = F), rC = 1 (free = F), thus both are fixed to the expected values when there is no qualitative sex dif
rAmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels = c("rao11"), free=F, values=.5, lbound=-.5, ubound=.5, name="Rao" )
rCmf  <- mxMatrix( type="Full", nrow=nv, ncol=nv, labels =c("rco11") , free=F, values=1, lbound=-1, ubound=1,name="Rco" )

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")
covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

# Algebra to compute standardized variance components
covM		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
covF		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAm		<- mxAlgebra( expression=Am/Vm, name="hm2")
StCm		<- mxAlgebra( expression=Cm/Vm, name="cm2")
StEm		<- mxAlgebra( expression=Em/Vm, name="em2")
StAf		<- mxAlgebra( expression=Af/Vf, name="hf2")
StCf		<- mxAlgebra( expression=Cf/Vf, name="cf2")
StEf		<- mxAlgebra( expression=Ef/Vf, name="ef2")

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, Am+Cm),		
                                           cbind(Am+Cm, Am+Cm+Em))  		, 						name="expCovMZM")
covMZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, Af+Cf),		
                                           cbind(Af+Cf, Af+Cf+Ef))  		, 						name="expCovMZF")
covDZM	<- mxAlgebra( expression= rbind  ( 	cbind(Am+Cm+Em, 0.5%x%Am+Cm),	
                                           cbind(0.5%x%Am+Cm, Am+Cm+Em))		, 					name="expCovDZM")
covDZF	<- mxAlgebra( expression= rbind  ( 	cbind(Af+Cf+Ef, 0.5%x%Af+Cf),	
                                           cbind(0.5%x%Af+Cf, Af+Cf+Ef))		, 					name="expCovDZF")
# corrected model below!! May 2020
covDZO <- mxAlgebra( expression= rbind  (     cbind( Am+Cm+Em, 				am%*%Rao%*%af + cm%*%Rco%*%cf ), 
                                              cbind( t(am%*%Rao%*%af) + t(cm%*%Rco%*%cf), 	Af+Cf+Ef) ), 	name="expCovDZO")                 

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, name="expCovMSing")

covFs		<- mxAlgebra( expression= Af+Cf+Ef 	, name="expCovFSing")

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDZO	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups
objmzm  	<- mxExpectationNormal( covariance="expCovMZM", means="expMeanMM", dimnames=selVars)
objdzm	<- mxExpectationNormal( covariance="expCovDZM", means="expMeanMM", dimnames=selVars)
objmzf	<- mxExpectationNormal( covariance="expCovMZF", means="expMeanFF", dimnames=selVars)
objdzf  	<- mxExpectationNormal( covariance="expCovDZF", means="expMeanFF", dimnames=selVars)
objdzo  	<- mxExpectationNormal( covariance="expCovDZO", means="expMeanMF", dimnames=selVars)
objsM	  	<- mxExpectationNormal( covariance="expCovMSing", means="expMeanM", dimnames=selVarsS)
objsF	  	<- mxExpectationNormal( covariance="expCovFSing", means="expMeanF", dimnames=selVarsS)

fitFunction 	<- mxFitFunctionML()

# Combine Groups
parsm		<- list( MeanM, pathAm, pathCm, pathEm, covAm, covCm, covEm, covM, StAm, StCm, StEm, fitFunction )
parsf		<- list( MeanF, pathAf, pathCf, pathEf, covAf, covCf, covEf, covF, StAf, StCf, StEf, fitFunction )
modelMZM	<- mxModel(parsm, MeanMM, covMZM, dataMZM, objmzm, name="MZM")
modelDZM	<- mxModel(parsm, MeanMM, covDZM, dataDZM, objdzm, name="DZM")
modelMZF	<- mxModel(parsf, MeanFF, covMZF, dataMZF, objmzf, name="MZF")
modelDZF	<- mxModel(parsf, MeanFF, covDZF, dataDZF, objdzf, name="DZF")
modelDZO	<- mxModel(parsm, parsf, MeanMF, rAmf, rCmf, covDZO, dataDZO, objdzo, name="DZO")
modelsM	<- mxModel(parsm, covMs, dataSM, objsM, name="SingM")
modelsF	<- mxModel(parsf, covFs, dataSF, objsF, name="SingF")

minus2ll	<- mxAlgebra(expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective + SingM.objective + SingF.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")

ciM  		<- mxCI (c ('MZM.hm2[1,1]','MZM.cm2[1,1]','MZM.em2[1,1]','DZO.rAmf [1,1]','DZO.rCmf [1,1]') )		# h2, c2, e2 males
ciF		<- mxCI (c ('MZF.hf2[1,1]','MZF.cf2[1,1]','MZF.ef2[1,1]') )		# h2, c2, e2 females

HetACE7Model	<-mxModel('HetACE7', modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, modelsM, modelsF, minus2ll, obj,ciM, ciF) 

# ------------------------------------------------------------------------------
# RUN Quantitative Heterogeneity ACE Model 
HetACE7Fit    <- mxRun(HetACE7Model, intervals=F)
(HetACE7Summ  <- summary(HetACE7Fit, verbose = F))

#Generate parameter estimates in a more convenient way, shorter

mxEval (MZM.hm2, HetACE7Fit)
mxEval (MZM.cm2, HetACE7Fit)
mxEval (MZM.em2, HetACE7Fit)

#Generate parameter estimates for females

mxEval (MZF.hf2, HetACE7Fit)
mxEval (MZF.cf2, HetACE7Fit)
mxEval (MZF.ef2, HetACE7Fit)

#---------------------------------------------------------------------
# Run Model II-2 Allow rCmf to be estimated
# Adds qualitative sex dif for C to the first model
#---------------------------------------------------------------------
HetACE7cModel	<-mxModel(HetACE7Fit, name="HetACE7c")
HetACE7cModel	<-omxSetParameters(HetACE7cModel, labels =c("rco11"), free=T, values=1, lbound=-1, ubound=1)
HetACE7cFit		<-mxTryHard(HetACE7cModel, intervals=F)
(HetACE7cSumm	<-summary(HetACE7cFit))

mxEval(MZM.Vm, HetACE7cFit)
mxEval(MZM.hm2, HetACE7cFit)
mxEval(MZM.cm2, HetACE7cFit)
mxEval(MZM.em2, HetACE7cFit)

mxEval(MZF.Vf, HetACE7cFit)
mxEval(MZF.hf2, HetACE7cFit)
mxEval(MZF.cf2, HetACE7cFit)
mxEval(MZF.ef2, HetACE7cFit)

#---------------------------------------------------------------------
# Run Model II-3 Allow only rAmf to be estimated
# Adds qualitative sex dif for A to the first model
#---------------------------------------------------------------------

HetACE7aModel	<-mxModel(HetACE7Fit, name="HetACE7a")
HetACE7aModel	<-omxSetParameters(HetACE7aModel, labels=  c("rao11"), free=T, values=0.5, lbound=-.5, ubound=.5)

HetACE7aFit		<-mxTryHard(HetACE7aModel, intervals=F)
(HetACE7aSumm	<-summary(HetACE7aFit, verbose = F))

mxEval(MZM.Vm, HetACE7aFit)
mxEval(MZM.hm2, HetACE7aFit)
mxEval(MZM.cm2, HetACE7aFit)
mxEval(MZM.em2, HetACE7aFit)

mxEval(MZF.Vf, HetACE7aFit)
mxEval(MZF.hf2, HetACE7aFit)
mxEval(MZF.cf2, HetACE7aFit)
mxEval(MZF.ef2, HetACE7aFit)

#---------------------------------------------------------------------
# Model II-4b - Specify and run the ACE homogeneity model for h2, c2 and e2 only
#---------------------------------------------------------------------

HomACE7Model	<-mxModel(HetACE7Fit, name="HomACE7")
HomACE7Model	<-omxSetParameters(HomACE7Model, 	labels=c("amp1"), free=T, values=c(.8), newlabels=c("afp1"))
HomACE7Model	<-omxSetParameters(HomACE7Model, 	labels=c("cmp1"), free=T, values=c(-.4), newlabels=c("cfp1"))
HomACE7Model	<-omxSetParameters(HomACE7Model, 	labels=c("efp1"), free=T, values=c(.6), newlabels=c("emp1"))
HomACE7Model	<-omxSetParameters(HomACE7Model, 	labels=c("mf1"), free=T, values=c(4.0), newlabels=c("mm1"))
HomACE7Model	<-omxAssignFirstParameters(HomACE7Model) 
HomACE7Fit		<-mxTryHard(HomACE7Model, intervals=F)
(HomACE7Summ	<-summary(HomACE7Fit))

#Generate parameter estimates (check that hf2, cf2 and ef2 will be the same as hm2, cm2 and em2)
mxEval(MZM.Vm, HomACE7Fit)
mxEval(MZM.hm2, HomACE7Fit)
mxEval(MZM.cm2, HomACE7Fit)
mxEval(MZM.em2, HomACE7Fit)

mxEval(MZF.Vf, HomACE7Fit)
mxEval(MZF.hf2, HomACE7Fit)
mxEval(MZF.cf2, HomACE7Fit)
mxEval(MZF.ef2, HomACE7Fit)


#--------------------------------------------------------------------------------------------------
# Model II-6f - Scalar (variance inequality) Model - Anxiety symptoms
# Modelling one set of ACE parameters for Males and Females 
# but specifying a scalar to multiply e.g. the female variances
# NOTE: This model only makes sense if you know from the Univariate Results that BOTH variables
# show a variance difference between males and females. If only one shows variance differences
# we need more complicated 'hybrid' models 
# --------------------------------------------------------------------------------------------------

nv		<- 1				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower	<- ntv*(ntv+1)/2 		# number of free elements in a lower matrix ntv*ntv
ncor		<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

selVars	<- c('Anx1', 'Anx2')	# THE VARIABLES FOR TWIN PAIRS
selVarsS 	<- c('Anx3')		# THE VARIABLES FOR A SINGLETON

# Subset data
# This will generate 5 sub data files for the 5 sex-by-zygosity groups
mzmData 	<- subdata3[subdata3$sexzyg == 1, selVars]
dzmData 	<- subdata3[subdata3$sexzyg == 2, selVars]
mzfData 	<- subdata3[subdata3$sexzyg == 3, selVars]
dzfData 	<- subdata3[subdata3$sexzyg == 4, selVars]
dzoData 	<- subdata3[subdata3$sexzyg == 5, selVars]

# Subset singletons
# This will generate 2 sub data files for the male adn female singelton groups
smData	<- subset(subdata3, sexzyg == 6, selVarsS)
sfData	<- subset(subdata3, sexzyg == 7, selVarsS)

# Check sub-data files
psych::describe(mzmData)
psych::describe(dzmData)
psych::describe(mzfData)
psych::describe(dzfData)
psych::describe(dzoData)
psych::describe(smData)
psych::describe(sfData)

# FIRST, CREATE OBJECTS WITH Labels for the Means, SD, and correlation to ease specification in the body of the model

# Create Labels for the Means Matrices
(LabMM	<- paste("mm",1:nv,sep=""))		# Male singleton
(LabMF	<- paste("mf",1:nv,sep=""))		# Female singleton
(LabMMM	<- c(LabMM,LabMM))			# male-male pairs
(LabMFF	<- c(LabMF,LabMF))			# female-female pairs
(LabMMF	<- c(LabMM,LabMF))			# male-female pairs

# Create Labels for the A, C and E Standard Deviation Matrices (paths from Latent Factors to variables) no sex dif
(LabAm	<- paste("apm",1:nv,sep=""))		
(LabCm	<- paste("cpm",1:nv,sep=""))		
(LabEm	<- paste("epm",1:nv,sep=""))		

(LabAf	<- paste("apm",1:nv,sep=""))		
(LabCf	<- paste("cpm",1:nv,sep=""))		
(LabEf	<- paste("epm",1:nv,sep=""))		

# Create Labels for the Scalar Matrix
(LabS		<- paste("s",1:nv,sep=""))		
(LabSS	<- c(LabS,LabS))
(LabNA	<- c("NA"))
(Labdos	<- c(LabNA,LabS))
(PatF		<- c(F))	
(PatT		<- c(T))
(Pat		<- c(PatF,PatT))

# Create objects for Start values
Stpath	<-c(.6) 			# change here according to the variances of the variables
Stcor		<-c(.6)
(Stmean	<-4.6)
StScalar	<-c(.8)				# change here to >1 if female variance is larger than males or if it is mixed
(DumOnes	<-c(1))
(StSdos	<-c(DumOnes,StScalar))		# Scalar for Males is just 1

# -------------------------------------------------------------------------------------------------- 
# Define matrices for the means

# Matrices to estimate Means for Male/Female pairs, OS pairs and singletons
MeanMM	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMM, 	name="expMeanMM" )
MeanM		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMM, 	name="expMeanM" )
MeanFF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMFF, 	name="expMeanFF" )
MeanF		<-mxMatrix( "Full", 1, nv, free=T, 		values=c(Stmean), 		labels=LabMF, 	name="expMeanF" )
MeanMF	<-mxMatrix( "Full", 1, ntv, free=T, 	values=c(Stmean,Stmean), 	labels=LabMMF, 	name="expMeanMF" )

pathAm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAm, name="am")
pathCm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCm, name="cm")
pathEm	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEm, name="em")

pathAf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabAf, name="af")
pathCf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabCf, name="cf")
pathEf	<- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=Stpath, label=LabEf, name="ef")

# Note, these are the correlations between the latent factors (not the variables)

# Algebra to compute the variance components
# Matrices generated to hold A, C and E computed Variance Components
covAm		<- mxAlgebra( expression=am %*% t(am), name="Am")
covCm		<- mxAlgebra( expression=cm %*% t(cm), name="Cm")
covEm		<- mxAlgebra( expression=em %*% t(em), name="Em")

covAf		<- mxAlgebra( expression=af %*% t(af), name="Af")
covCf		<- mxAlgebra( expression=cf %*% t(cf), name="Cf")
covEf		<- mxAlgebra( expression=ef %*% t(ef), name="Ef")

covAmf	<- mxAlgebra( expression=am %*% t(af), name="Amf")
covCmf	<- mxAlgebra( expression=cm %*% t(cf), name="Cmf")

# Algebra to compute standardized variance components
covm		<- mxAlgebra( expression=Am+Cm+Em, name="Vm")
StAm		<- mxAlgebra( expression=Am/Vm, name="h2m")
StCm		<- mxAlgebra( expression=Cm/Vm, name="c2m")
StEm		<- mxAlgebra( expression=Em/Vm, name="e2m")

covf		<- mxAlgebra( expression=Af+Cf+Ef, name="Vf")
StAf		<- mxAlgebra( expression=Af/Vf, name="h2f")
StCf		<- mxAlgebra( expression=Cf/Vf, name="c2f")
StEf		<- mxAlgebra( expression=Ef/Vf, name="e2f")

ScalarF	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=T, values=c(StScalar,StScalar), label=LabSS, name="ScF")
ScalarOS	<- mxMatrix( type="Diag", nrow=2, ncol=2, free=Pat, values=StSdos, label=Labdos, name="ScOS")
ScalarFs	<- mxMatrix( type="Diag", nrow=1, ncol=1, free=T, values=StScalar, label=LabS, name="ScFs")

# Algebra for expected variance/covariance matrix in the 4 groups
covMZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, Am+Cm),	cbind(Am+Cm, Am+Cm+Em)), name="ExpCovMZM")
covDZM	<-mxAlgebra( expression= rbind  ( cbind(Am+Cm+Em, .5%x%Am+Cm),	cbind(.5%x%Am+Cm, Am+Cm+Em)	), name="ExpCovDZM")

covMZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, Af+Cf),		cbind(Af+Cf, Af+Cf+Ef))),	name="ExpCovMZF")
covDZF	<-mxAlgebra( expression= ScF %&% (rbind  ( cbind(Af+Cf+Ef, .5%x%Af+Cf),	cbind(.5%x%Af+Cf, Af+Cf+Ef))),	name="ExpCovDZF")

covDOS	<-mxAlgebra( expression= ScOS %&% (rbind  ( cbind(Am+Cm+Em, .5%x%Amf+Cmf),	cbind(.5%x%Amf+Cmf, Af+Cf+Ef)	)),	name="ExpCovDOS")

covMs		<- mxAlgebra( expression= Am+Cm+Em 	, 			name="expCovMSing")
covFs		<- mxAlgebra( expression= ScFs %&% (Af+Cf+Ef) 	, 	name="expCovFSing")    

# Data objects for Multiple Groups
dataMZM	<-mxData(mzmData, type="raw")
dataMZF	<-mxData(mzfData, type="raw")
dataDZM	<-mxData(dzmData, type="raw")
dataDZF	<-mxData(dzfData, type="raw")
dataDOS	<-mxData(dzoData, type="raw")
dataSM	<-mxData(smData, type="raw")
dataSF	<-mxData(sfData, type="raw")

# Objective objects for Multiple Groups     
objMZM	<- mxExpectationNormal(covariance="ExpCovMZM", means="expMeanMM", dimnames=selVars)
objDZM	<- mxExpectationNormal(covariance="ExpCovDZM", means="expMeanMM", dimnames=selVars)
objMZF	<- mxExpectationNormal(covariance="ExpCovMZF", means="expMeanFF", dimnames=selVars)
objDZF	<- mxExpectationNormal(covariance="ExpCovDZF", means="expMeanFF", dimnames=selVars)
objDOS	<- mxExpectationNormal(covariance="ExpCovDOS", means="expMeanMF", dimnames=selVars)
objsM		<- mxExpectationNormal(covariance="expCovMSing", means="expMeanM", dimnames=selVarsS )
objsF		<- mxExpectationNormal(covariance="expCovFSing", means="expMeanF", dimnames=selVarsS )

fitFunction <- mxFitFunctionML()

# Combine Groups
parsm	 	<-list(pathAm, pathCm, pathEm, covAm, covCm, covEm, covm, StAm, StCm, StEm, fitFunction)
parsf	 	<-list(pathAf, pathCf, pathEf, covAf, covCf, covEf, covf, StAf, StCf, StEf, fitFunction)
modelMZM 	<-mxModel(parsm, MeanM, MeanMM, covMZM, dataMZM, objMZM,  name="MZM")
modelDZM 	<-mxModel(parsm, MeanM, MeanMM, covDZM, dataDZM, objDZM,  name="DZM")
modelMZF 	<-mxModel(parsf, MeanF, MeanFF, covMZF, dataMZF, objMZF, ScalarF, name="MZF")
modelDZF 	<-mxModel(parsf, MeanF, MeanFF, covDZF, dataDZF, objDZF, ScalarF, name="DZF")
modelDOS	<-mxModel(parsm, parsf, MeanM, MeanF, MeanMF, covAmf, covCmf, covDOS, dataDOS, objDOS, ScalarOS, name="DOS")
modelsM	<-mxModel(parsm, MeanM, covMs, dataSM, objsM, name="SM")
modelsF	<-mxModel(parsf, MeanF, covFs, dataSF, objsF, ScalarFs, name="SF")

minus2ll  <-mxAlgebra(expression=	MZM.objective + DZM.objective + MZF.objective + DZF.objective 
                      			+ DOS.objective + SM.objective + SF.objective, name="m2LL")
obj		<-mxFitFunctionAlgebra("m2LL")
ci		<-mxCI (c ('MZM.h2m[1,1]','MZM.c2m[1,1]','MZM.e2m[1,1]','MZF.h2f[1,1]','MZF.c2f[1,1]','MZF.e2f[1,1]','SF.ScFs[1,1]') )
ScACE7Model	<-mxModel("ScACE7", parsm, parsf, modelMZM, modelDZM, modelMZF, modelDZF, modelDOS, modelsM, modelsF, minus2ll, obj, ci) 

# -------------------------------------------------------------------------------
# RUN Scalar ACE Model VII-5
ScACE7Fit	<-mxTryHard(ScACE7Model, intervals=F)
(ScACE7Sum 	<-summary(ScACE7Fit))

mxEval (MZM.h2m, ScACE7Fit)
mxEval (MZM.c2m, ScACE7Fit)
mxEval (MZM.e2m, ScACE7Fit)
mxEval (MZM.am, ScACE7Fit)
mxEval (MZM.cm, ScACE7Fit)
mxEval (MZM.em, ScACE7Fit)

mxEval (MZF.h2f, ScACE7Fit)
mxEval (MZF.c2f, ScACE7Fit)
mxEval (MZF.e2f, ScACE7Fit)
mxEval (MZF.af, ScACE7Fit)
mxEval (MZF.cf, ScACE7Fit)
mxEval (MZF.ef, ScACE7Fit)
mxEval (MZF.ScF, ScACE7Fit)

#----------------------------------------------------------------
# Print Comparative Fit Statistics between models
#----------------------------------------------------------------
mxCompare(HetACE7cFit, HetACE7Fit)
mxCompare(HetACE7aFit, HetACE7Fit)
mxCompare(HetACE7Fit, HomACE7Fit)
mxCompare(HetACE7Fit, ScACE7Fit)



