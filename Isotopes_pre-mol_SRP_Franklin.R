## SIBER CODE UPDATED FROM: 
##https://cran.r-project.org/web/packages/SIBER/vignettes/siber-comparing-populations.html

rm(list=ls())
graphics.off()

set.seed(1)


library(SIBER)
library(ggplot2)
library(dplyr)


# load data
mydata<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Isotopos/PPA_GLS.csv",header=T)

mydata2<-dplyr::filter(mydata, tissue %in% c("feather"))

mydata3 <- mydata2 %>%
  rename(iso1 = d13C,
         iso2 = d15N,
         group = sex)%>%
  mutate(community = 1)%>%
  select(iso1,iso2, group,community)

mydata3 <- mydata3 %>%
  arrange(group)

# create the siber object
siber.example<- createSiberObject(mydata3)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 3, lwd = 3)
group.ellipses.args  <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "cyan", "violet")

par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
                x.limits = c(-26,-18),
                y.limits = c(7,16)
)

legend("bottomright",
       legend = as.character(paste("sex",unique(mydata3$group))),
       pch = 19,
       col = 1:length(unique(mydata3$group)),
       bty = "n")


# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siber.example) 
print(community.ML)

group.ML <- groupMetricsML(siber.example)
group.ML
write.csv(group.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3


# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)

SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = c("F","M"), 
                 ylims = c(0,8),
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)
data



cr.p <- c(0.95) # vector of quantiles

SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

SEA.B.credibles
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)


# ---------------------------------------------------------------------
# Compare the posterior distributions
# ---------------------------------------------------------------------

#females vs Males
Pg1.lt.g2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.lt.g2) #[1] 0.80425


# ---------------------------------------------------------------------
# Overlap between ellipses
# ---------------------------------------------------------------------

group.ML 

#FEMALES
overlap.g1.g2 <- maxLikOverlap("1.F", "1.M", siber.example, p = 0.95)


#the proportion of females that overlaps with males
prop.of.first1.2 <- as.numeric(overlap.g1.g2["overlap"] / overlap.g1.g2["area.1"])
print(prop.of.first1.2)#[1] 0.378955


#the proportion of males that overlaps with females
prop.of.second1.2 <- as.numeric(overlap.g1.g2["overlap"] / overlap.g1.g2["area.2"])
print(prop.of.second1.2)#[1] 0.552773


#the proportion of females and males that overlap with each other
prop.of.both <- as.numeric(overlap.g1.g2["overlap"] / (overlap.g1.g2["area.1"] + overlap.g1.g2["area.2"]))
print(prop.of.both)#[1] 0.2248254


bayes.overlap.G1.G2 <- bayesianOverlap("1.F", "1.M", ellipses.posterior, 
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
print(bayes.overlap.G1.G2)


################## Comparisons between sex isotopic values #########

library(onewaytests)

library(ggplot2)
library(lme4)
library(lattice)
library(nlme)
library(geepack)
library(MuMIn)
library(sjPlot)


#### Normality ####

hist(mydata3$iso1)
hist(mydata3$iso2)


lm_d13C.M<-lm(iso1~group,data=mydata3, na.action=na.fail)
plot(lm_d13C.F)#qqplot ok

lm_d15.M<-lm(iso2~group,data=mydata3, na.action=na.fail)
plot(lm_d15.M)#qqplot ok


## We model the heterocedacea using the varIdent function

vf2 <- varIdent(form= ~ 1 | group)


#We run a model for each isotope separately.
#Carbono:

M.gls1<-gls(iso1~group,weights = vf2, data=mydata3, method="ML", na.action=na.fail)

summary(M.gls1)
dredge(M.gls1)

anova(M.gls1, test = "F")
tab_model(M.gls1)

sum(mydata3$iso1 & mydata3$group == "F")
mean(mydata3$iso1[mydata3$group == "F"])
sd(mydata3$iso1[mydata3$group == "F"])
range(mydata3$iso1[mydata3$group == "F"])

sum(mydata3$iso1 & mydata3$group == "M")
mean(mydata3$iso1[mydata3$group == "M"])
sd(mydata3$iso1[mydata3$group == "M"])
range(mydata3$iso1[mydata3$group == "M"])


#Nitrogen:

M.gls2<-gls(iso2~group,weights = vf2, data=mydata3, method="ML", na.action=na.fail)

summary(M.gls2)
dredge(M.gls2)

anova(M.gls2, test = "F")
tab_model(M.gls2)

sum(mydata3$iso2 & mydata3$group == "F")
mean(mydata3$iso2[mydata3$group == "F"])
sd(mydata3$iso2[mydata3$group == "F"])
range(mydata3$iso2[mydata3$group == "F"])

sum(mydata3$iso2 & mydata3$group == "M")
mean(mydata3$iso2[mydata3$group == "M"])
sd(mydata3$iso2[mydata3$group == "M"])
range(mydata3$iso2[mydata3$group == "M"])

# posthoc test for multiple comparisons 
library(rstatix)

games_howell_test(data=mydata3, iso1~group, conf.level = 0.95, detailed = FALSE)
games_howell_test(data=mydata3, iso2~group, conf.level = 0.95, detailed = FALSE)


################### Trophic Position ##########################
# tRophicPosition Rockhopper penguins, code addapted from https://github.com/clquezada/tRophicPosition

rm(list=ls())
graphics.off()

#install the latest version 
library(usethis)
library(devtools)
library(tRophicPosition)
library(rjags)

library(SIBER)

#set you working directory
setwd("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Isotopos/trophic_position")

######################################### FEMALES ################################################

# 2020 pre-molt 

#Loading data from a csv/txt file
mydata1 <- read.table("PPA_Female_GLS_TP_cherel.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE)
#mydata1 <- read.csv("C:/Users/User/Desktop/PPA_project/Isotopos/TP_files/Females/PPA_Female_2014_TP.csv", header=TRUE, sep="\t", dec=".", strip.white=TRUE)


# We need to separate the baselines (if there are 2)
# First we select d15N and d13C for baseline 1
dNb1 <- mydata1$d15N[which(mydata1$FG == "Baseline" & mydata1$Spp == "Eufa_CA")]
dCb1 <- mydata1$d13C[which(mydata1$FG == "Baseline" & mydata1$Spp == "Eufa_CA")]

# Then we select d15N and d13C for baseline 2
dNb2 <- mydata1$d15N[which(mydata1$FG == "Baseline" & mydata1$Spp == "Eufa_PFZ")]
dCb2 <- mydata1$d13C[which(mydata1$FG == "Baseline" & mydata1$Spp == "Eufa_PFZ")]

dNb1
dCb1
dNb2
dCb2

#Now we want to extract the target species. This is done with a similar code, 
#but changing the Species value:
dN_sp1 <- mydata1$d15N[which(mydata1$Spp == "Rockhopper_chrysocome_F")]
dC_sp1 <- mydata1$d13C[which(mydata1$Spp == "Rockhopper_chrysocome_F")]

dN_sp1 
dC_sp1

#2 baselines
sp1 <- list("dNb1" = dNb1, "dCb1" = dCb1, "dNb2" = dNb2,
            "dCb2" = dCb2, "dNc" = dN_sp1, "dCc" = dC_sp1)
#Generating TEFs
TDF_CHEREL <- simulateTDF(nN = 10, meanN = 4.4, sdN = 0.5,
                        nC = 10, meanC = 0.1, sdC = 0.5)
TDF_CHEREL

sp1$deltaN <- TDF_CHEREL$deltaN
sp1$deltaC <- TDF_CHEREL$deltaC

str(sp1)

#Defining the Bayesian model for trophic position

model.string <- jagsTwoBaselinesFull()


#Initializing the model and sampling the posterior estimates

IsotopeData <- lapply(sp1, na.omit)


screenIsotopeData(isotopeData = sp1, density = "both",
                  consumer = "Rockhopper_chrysocome_F", b1 = "Eufa_CA",
                  b2 = "Eufa_PFZ", legend = c(1.15, 1.15), title = NULL)

model.sp1 <- TPmodel(data = sp1,
                     model.string = model.string,
                     n.chains = 4,
                     n.adapt = 200000)
model.sp1

samples.sp1 <- posteriorTP(model.sp1)


# summarize the posterior data
summary(samples.sp1)
# plot the posterior data:
plot(samples.sp1)
# save the posterior samples of trophic position (only first chain)
sp1.TP <- as.data.frame(samples.sp1[[1]][,"TP"])$var1
# combine it with the posterior samples from the second chain
sp1.TP <- c(sp1.TP,as.data.frame(samples.sp1[[2]][,"TP"])$var1)
# check the length of the variable created
length(sp1.TP)


Species <- c(rep("Rockhopper_chrysocome_F", length(sp1.TP)))
# and then we combine it with the posterior samples of the species trophic position
df1 <- data.frame(sp1.TP, Species)


colnames(df1) <- c("TP", "Species")

# check everything worked well
summary(df1)

# plot it:

trophicDensityPlot(df1)
trophicDensityPlot(df1, quantiles = TRUE)

# use siberDensityPlot() function from SIBER packege:
# First show the posterior samples into a data frame
dfBB <- data.frame(sp1.TP)
# change the name of each column
colnames(dfBB) <- c("Rockhopper_chrysocome_F")
# plot the data
SIBER::siberDensityPlot(dfBB, xlab = "Species",
                        ylab = "Trophic Position")

################################## MALES #################################3


mydata2 <- read.table("PPA_Male_GLS_TP_cherel.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE)

dNb1 <- mydata2$d15N[which(mydata2$FG == "Baseline" & mydata2$Spp == "Eufa_CA")]
dCb1 <- mydata2$d13C[which(mydata2$FG == "Baseline" & mydata2$Spp == "Eufa_CA")]

# Baseline 2
dNb2 <- mydata2$d15N[which(mydata2$FG == "Baseline" & mydata2$Spp == "Eufa_PFZ")]
dCb2 <- mydata2$d13C[which(mydata2$FG == "Baseline" & mydata2$Spp == "Eufa_PFZ")]


dN_sp2 <- mydata2$d15N[which(mydata2$Spp == "Rockhopper_chrysocome_M")]
dC_sp2 <- mydata2$d13C[which(mydata2$Spp == "Rockhopper_chrysocome_M")]

#2 baselines
sp2 <- list("dNb1" = dNb1, "dCb1" = dCb1, "dNb2" = dNb2,
            "dCb2" = dCb2, "dNc" = dN_sp2, "dCc" = dC_sp2)
#TDF
sp2$deltaN <- TDF_CHEREL$deltaN
sp2$deltaC <- TDF_CHEREL$deltaC

str(sp2)

model.string <- jagsTwoBaselinesFull()


IsotopeData <- lapply(sp2, na.omit)


screenIsotopeData(isotopeData = sp2, density = "both",
                  consumer = "Rockhopper_chrysocome_M", b1 = "Eufa_CA",
                  b2 = "Eufa_PFZ", legend = c(1.15, 1.15), title = NULL)


# choose the model:

model.sp2 <- TPmodel(data = sp2,
                     model.string = model.string,
                     n.chains = 4,
                     n.adapt = 200000)
# posteriorTP() will monitor by default TP and muDeltaN.

samples.sp2 <- posteriorTP(model.sp2)

summary(samples.sp2)
plot(samples.sp2)

# save the posterior samples of trophic position (only first chain)
sp2.TP <- as.data.frame(samples.sp2[[1]][,"TP"])$var1
# combine it with the posterior samples from the second chain
sp2.TP <- c(sp2.TP,as.data.frame(samples.sp2[[2]][,"TP"])$var1)
# check the length of the variable created
length(sp2.TP)

# create the data frame:
# With this code we make the Species variable, including the species name
Species <- c(rep("Rockhopper_chrysocome_M", length(sp2.TP)))
# and then we combine it with the posterior samples of the species trophic position
df2 <- data.frame(sp2.TP, Species)

# We need the data frame to include the variable names TP and Species, 
# so we need to change the names accordingly:
# Changing the variable names of the dataframe "df"

colnames(df2) <- c("TP", "Species")

# check everything worked well
summary(df2)

# plot it:
trophicDensityPlot(df2, quantiles = TRUE)

# use siberDensityPlot() function from SIBER packege:
# First show the posterior samples into a data frame
dfBB <- data.frame(sp2.TP)
# change the name of each column
colnames(dfBB) <- c("Rockhopper_chrysocome")
# plot the data
SIBER::siberDensityPlot(dfBB, xlab = "Species",
                        ylab = "Trophic Position")

####################Comparing two or more trophic positions##########################
# To compare the TP of 2 species use the function compareTwoDistributions(). 
# This function receives 2 posterior distributions and a test, and returns the
# probability of ocurring that comparison, in a Bayesian way. 

# posterior trophic position:

summary(sp1.TP)
summary(sp2.TP)

sd(sp1.TP)
sd(sp2.TP)


# We want to know whether there is any statistical support for differences in TP 
# between the species:
# We can inquire wheter "the posterior distribution of trophic position of 
# sp 1 is higher than the posterior distribution of sp 2".

# Female vs Male
compareTwoDistributions(sp1.TP,sp2.TP, test = ">=")
# [1] 0.473
#[1] 0.497 Cherel

# build a data frame combining all TP and creating a Species factor

# Males Feathers PPA
TP <- c(sp1.TP,sp2.TP)
Species <- c(rep("Female", length(sp1.TP)),
             rep("Male", length(sp1.TP)))
df_all_3 <- data.frame(TP, Species)


# use the function for plotting trophic position grouped by the Species factor
# and without quantiles (the default)

trophicDensityPlot(df_all_3)

p <- trophicDensityPlot(df_all_3)


library(ggplot2)
p + scale_y_continuous(limits = c(0, 2), breaks = seq(0, 5, by = 1), labels = seq(0, 5, by = 1))

# use siberDensityPlot() function from SIBER packege:
# First combine the two posterior samples into a data frame

# Females Feather PPA

df_all4 <- data.frame(sp1.TP,sp2.TP)
colnames(df_all4) <- c("Female", "Male")

SIBER::siberDensityPlot(df_all4, xlab = "Sex",
                        ylab = "Trophic Position",ylim=c(1,5),prn = TRUE)









