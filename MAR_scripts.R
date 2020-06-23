##################################################################
# R-script for data analysis and visualization pertaining to Dexter
# et. al manuscript "Modeling the trophic impacts of invasive 
# zooplankton in a highly invaded river"
# Last updated 23 June 2020
##################################################################
#Load required packages
library(Hmisc)
library(MARSS)
library(MAR1)

#Load data
MAR_input <- read.csv("C:/Users/Ericd/Dropbox/Eric Work/PhD/MAR/data/MAR_input.csv")

#Aggregate autotrophic taxa into a single variable
MAR_input$autotrophs<- MAR_input$Diatom+MAR_input$Flagellate + MAR_input$Green_algae + MAR_input$Cyanobact

#Format data for model input using MAR1 package
mar1 <- prepare.data(data=MAR_input,increment="month",z.method="deseason",fill.gap = 1,log = FALSE,order="ymd")

#Create t+1 lagged variables for C matrix
mar1$Corbicula.lag <- Lag(mar1$Corbicula, +1)
mar1$Forb_Cope_lag <- Lag(mar1$Forbesi_Cope, +1)
mar1$Forb_Naup_lag <- Lag(mar1$Forbesi_Naup, +1)

#Drop first row (not needed)
mar1<-mar1[-1,]

################################################################################
#M1 - Full model
###############################################################################

##Define a matrix of variates (the elements of the B Matrix)
variate.list<- c(12,14,18,19,20,25,28)
Var.mat <- cbind(mar1[,c(variate.list)])
names(Var.mat) <- colnames(mar1[c(variate.list)])
Var.mat <- t(Var.mat)

#Now set up the matrix of covariates.
covariate.list <-c(5,29,30,31)
Cov.mat <- cbind(mar1[,c(covariate.list)])                                   
names(Cov.mat) <- colnames(mar1[c(covariate.list)])
Cov.mat <- t(Cov.mat)

##Now Define the model - 
# In this case we define the structure of the terms in the Process and Observation equations of the 
# state space MAR.  If we make the observation process zero (A, D and R from equation 1.1 in the MARSS 
# users manual) then we are left with only the process equation - the "vanilla" MAR1.
# C matrix is unconstrained, so we can estimate co-variate effects.
cntl.list=list(safe=TRUE, maxit=1500, allow.degen=TRUE)
UU="zero"
BB="unconstrained"
CC="unconstrained"
QQ="diagonal and unequal"

#For Observation Equation Only
RR="zero"
ZZ="identity"
AA="zero"
DD="zero"
dd="zero"

CRD.model.1=list(Z=ZZ, U=UU, Q=QQ, R=RR, B=BB, C=CC, d=dd, D=DD, A=AA)
CRD.model.1$c=Cov.mat

#Run the model
CRD.MARSSC1 <- MARSS(as.matrix(Var.mat),model=CRD.model.1, control=cntl.list, silent=2)
m1<-CRD.MARSSC1

#Get bootstrapped CI's for params
CRD.MARSSC1.CI <- MARSSparamCIs(CRD.MARSSC1) 

#Now pull the results out of the results object so we can plot them
BBB <-  matrix(CRD.MARSSC1$coef[1:length(variate.list)^2],nrow=length(variate.list), ncol=length(variate.list))
colnames(BBB) <- rownames(BBB) <- colnames(mar1[c(variate.list)])  ## Define the matrix of co-variates (the elements of the C Matrix)

cov.index.start <- length(CRD.MARSSC1$coef)+1 - length(covariate.list)*length(variate.list)
cov.index.end <- length(CRD.MARSSC1$coef)
CCC <-  matrix(CRD.MARSSC1$coef[cov.index.start:cov.index.end],nrow=length(variate.list), ncol=length(covariate.list))
rownames(CCC) <- colnames(mar1[c(variate.list)])    ## Define the matrix of co-variates (the elements of the C Matrix)
colnames(CCC) <- colnames(mar1[c(covariate.list)])

##############################################################################
#Extract the bootstrapped CI's for plotting  - Pull the data then make them a matrix, then add them to the list of data to be plotted.
#First do this for the B matrix

for (i in 1:length(variate.list)^2) {
  if ((CRD.MARSSC1.CI$par.lowCI$B[i]>0)|(CRD.MARSSC1.CI$par.upCI$B[i]<0) ) 
  {CRD.MARSSC1$bootstrap$B[i] <- CRD.MARSSC1$coef[i]}
  else {CRD.MARSSC1$bootstrap$B[i] <- 0}
}

CRD.MARSSC1$bootstrap$B <- matrix(CRD.MARSSC1$bootstrap$B[1:length(variate.list)^2], nrow=length(variate.list), ncol=length(variate.list))
str(CRD.MARSSC1$bootstrap$B)
rownames(CRD.MARSSC1$bootstrap$B) <- colnames(CRD.MARSSC1$bootstrap$B) <- colnames(mar1[c(variate.list)])

#Now do the same thing for the C matrix
Cindex <-(length(covariate.list)*length(variate.list))
Cstart <-length(CRD.MARSSC1$coef)-Cindex
for (i in 1:Cindex) {
  if ((CRD.MARSSC1.CI$par.lowCI$U[i]>0)|(CRD.MARSSC1.CI$par.upCI$U[i]<0) ) 
  {CRD.MARSSC1$bootstrap$C[i] <- CRD.MARSSC1$coef[i+Cstart]}
  else {CRD.MARSSC1$bootstrap$C[i] <- 0}
}
CRD.MARSSC1$bootstrap$C <- matrix(CRD.MARSSC1$bootstrap$C[1:(length(covariate.list)*length(variate.list))], nrow=length(variate.list), ncol=length(covariate.list))

rownames(CRD.MARSSC1$bootstrap$C) <- colnames(mar1[c(variate.list)])
colnames(CRD.MARSSC1$bootstrap$C) <- colnames(mar1[c(covariate.list)])

#Now create the MAR plot object and plot it.
myres <- matrix(0.5, nrow=length(variate.list),ncol=length(variate.list)+length(covariate.list)) #No restrictions

CRD.MARSSC1.plot <-
  list(restrictions.set=myres,
       bestfit=list(B=BBB,C=CCC), bootstrap=CRD.MARSSC1$bootstrap)
class(CRD.MARSSC1.plot) <- "MAR"

plot(CRD.MARSSC1.plot,legend=F)

#Export as 8X5 PDF
############################################################################
#Custom plots
####################################################################
library(ggplot2)
library(reshape)

Clow <- CRD.MARSSC1.CI$par.lowCI$U
Chigh <-  CRD.MARSSC1.CI$par.upCI$U
Ccoef <- CRD.MARSSC1$coef[64:91]
noSig <-ifelse(Clow < 0 & Chigh > 0, 0, Ccoef)
temp <- data.frame(Clow,Chigh,Ccoef,noSig)
temp2 <- temp[which(temp$noSig != 0),]
xgroup<-as.factor(c("Temp","Temp","P. forbesi","P. forbesi","Nauplii","Nauplii"))
ygroup<-as.factor(c("Bosmina","Daphnia","Daphnia","cyclo","Ciliate","Autotroph"))
temp2 <-cbind(xgroup,ygroup,temp2)
ggplot(data=temp2, aes(x=xgroup, y=Ccoef,group=ygroup)) +
  geom_errorbar(aes(ymin=temp2$Clow, ymax=temp2$Chigh),width=.5,position=position_dodge(.9)  ) +
  geom_bar(stat="identity", position="dodge", fill="darkgrey",color="black",width=1)+ ylim(-0.5,0.5) +
  geom_text(aes(label=ygroup), position=position_dodge(width=.9), size=4,angle = 90)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size = 18))
  #save as 12x5 PDF

lowerCI<-CRD.MARSSC1.CI$par.lowCI$U
upperCI<-CRD.MARSSC1.CI$par.upCI$U
temp<-melt(as.matrix(CCC))
temp2<-data.frame(temp,lowerCI,upperCI)
specA<-temp$X1
specB<-temp$X2
midCI<-temp$value
df<-data.frame(specA,specB,lowerCI,upperCI,midCI)
df<-df[which(df$midCI != 0),] #Remove interactions where CI includes zero
df<-df[which(specA != specB),]

lowerCI<-melt(as.matrix(lowerB[,1:7]))
upperCI<-melt(as.matrix(upperB[,1:7]))
upperCI<-upperCI$value
midCI<-melt(as.matrix(midB[,1:7]))
midCI<-midCI$value
#includeZero<-between(0,lowerCI$value,upperCI)
df<-data.frame(lowerCI,upperCI,midCI)
df<-df[which(df$midCI != 0),] #Remove interactions where CI includes zero
df<-df[which(specA != specB),]

specA<-as.character(df$X1)
specB<-as.character(df$X2)
df<-df[which(specA != specB),]

ggplot(data=df, aes(x=specB, y=midCI,group=specA)) +
  geom_errorbar(aes(ymin=df$lowerCI, ymax=df$upperCI),width=.5,position=position_dodge(.9)  ) +
  geom_bar(stat="identity", position="dodge", fill="darkgrey",color="black",width=1)+ ylim(-0.5,0.5) +
  geom_text(aes(label=specB), position=position_dodge(width=.9), size=4,angle = 90)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size = 18))
#save as 12x5 PDF

##################################################################
#################################################################################
