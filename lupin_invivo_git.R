multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

setwd("/path/to/in_vivo/Stage_CamilleAllanic/Analyses/")
#setwd('//trs-nas1/caallanic/Mes Documents/Camille INRA/Th?se/Projet LUPINEMA/R?sultats/Ch?vres/')
require(nlme)
require(ade4)
require(GGally)
# require(lmer)
# require(lmerTest)
require(Hmisc)
require(car)
require(lsmeans)
require(RVAideMemoire)
library("FactoMineR")
library("factoextra")
theme_set(theme_bw())
### Time points
tp = c('J0','J18','J21','J24','J30')

##======================================================================================
## GOATS 
##======================================================================================

dat = read.csv(file='Goat_res.csv', header=T, sep=";")
head(dat)
dat$Jour_Prlv=as.character(dat$Jour_Prlv)
dat$Jour_Prlv[dat$Jour_Prlv=='J31']='J30'
goat = dat[dat$Infestation==2 
         & dat$Jour_Prlv %in% tp,
         c("ID_ANIMAL",'Ht','FEC','PL','TB','TP','scs','Jour_Prlv','Lot','Lignee')]
goat$Jour_Prlv=factor(goat$Jour_Prlv)
goat$Group = 'Lup-Inf'
goat$Group[goat$Lot=='ConcInf'] = 'Conc-Inf'
goat$Group[goat$Lot=='ConcNInf'] = 'Conc-Ninf'
goat$Group[goat$Lot=='LupNInf'] = 'Lup-Ninf'
goat$Group = as.factor(goat$Group)
goat$Day = gsub('J','D',goat$Jour_Prlv)

## Prepatent period
ppg = dat[!is.na(dat$Tps_exc) & dat$Infestation==2 ,c('ID_ANIMAL','Tps_exc','Lot'),]

##======================================================================================
## EWES 
##======================================================================================
rm(dat)
dat=read.csv(file='Ewe_res.csv', header=T, sep=";")
dat$Infestation = factor(dat$Infestation)

ewe = dat[which(dat$Infestation==2 & dat$Jour_Prlv %in% tp),
           c("ID_ANIMAL","Lot", "Infestation", "Jour_Prlv","Ht", "FEC","Pepsino", "Poids", "GMQtot", "GMQinst")]
ewe$Jour_Prlv=factor(ewe$Jour_Prlv)
ewe$species='ewe'
ewe$Group = 'Lup-Inf'
ewe$Group[ewe$Lot=='ConcInf'] = 'Conc-Inf'
ewe$Group[ewe$Lot=='ConcNInf'] = 'Conc-Ninf'
ewe$Group[ewe$Lot=='LupNInf'] = 'Lup-Ninf'
ewe$Group = as.factor(ewe$Group)
ewe$Day = gsub('J','D',ewe$Jour_Prlv)

## Prepatent period
ppe = dat[!is.na(dat$Tps_exc) & dat$Infestation==2 ,c('ID_ANIMAL','Tps_exc','Lot'),]

##======================================================================================
## PLOT FEC, HT, GW or MILK for each group
##======================================================================================
require(RColorBrewer)
#require(ggridges)

infewe = ewe[ewe$Lot %in% c('LupInf','ConcInf'),]
infewe$Lot = factor(infewe$Lot)
infewe$Group = factor(infewe$Group)
vecol = c('#e66101','#fdb863','#5e3c99','#b2abd2')
scl = c(0,50,500,5000,25000)

p1 = ggplot(infewe,aes(x = Day,y = FEC+1, gp = Group, fill = Group)) +
  geom_violin() +
  #geom_point(position=position_dodge(.8)) +
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(.8), binwidth = 100) + 
  scale_fill_manual(values = vecol[c(1,3)]) + 
  scale_y_log10(limits=c(1,35000),breaks = c(1,50,500,5000,25000),labels= scl) +
  theme(text = element_text(size=14),legend.position='none') +
  xlab('Days post-infection') + ylab('Faecal Egg Count (eggs/g)') +ggtitle('a')

p2 = ggplot(infewe,aes(x = Day,y = Ht, gp = Group, fill = Group)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(.6), binwidth = 0.8) + 
  scale_fill_manual(values = vecol[c(1,3)]) + 
  theme(text = element_text(size=14),legend.position='none') +
  xlab('Days post-infection') + ylab('Haematocrit (%)') +ggtitle('c')
ewe$stat = sapply(strsplit(as.character(ewe$Group),'-'),function(x) x[2])

p3 = ggplot(ewe[!is.na(ewe$GMQtot),],aes(x = Group,y = GMQtot, gp = Day, fill = Group)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(.6), binwidth = 30) + 
  scale_fill_manual(values = vecol) + 
  theme(text = element_text(size=14),legend.position = 'bottom')+ 
  scale_y_continuous(limits = c(-300,200),breaks=seq(-300,200,100)) +
  xlab('Days post-infection') + ylab('Average daily gain (g/day)') + ggtitle('e')

infgoat = goat[goat$Lot %in% c('LupInf','ConcInf'),]
infgoat$Lot = factor(infgoat$Lot)
infgoat$Group = 'Lup-Inf'
infgoat$Group[infgoat$Lot=='ConcInf'] = 'Conc-Inf'
infgoat$Group = as.factor(infgoat$Group)
infgoat$Day = gsub('J','D',infgoat$Jour_Prlv)
scl = c(0,50,500,1500,4500)
vecol=c('#a6611a','#dfc27d','#80cdc1','#018571')

p4 = ggplot(infgoat,aes(x = Day,y = FEC+1, gp = Group, fill = Group)) +
  geom_violin() +
  #geom_point(position=position_dodge(.8)) +
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(.8), binwidth = 100) + 
  scale_fill_manual(values = vecol[c(1,3)]) + 
  scale_y_log10(limits=c(1,6400),breaks = c(1,50,500,1500,4500),labels= scl) +
  theme(text = element_text(size=14),legend.position='none') +
  xlab('Days post-infection') + ylab('Faecal Egg Count (eggs/g)') + ggtitle ('b')  

p5 = ggplot(infgoat,aes(x = Day,y = Ht, gp = Group, fill = Group)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(.6), binwidth = 0.6) + 
  scale_fill_manual(values = vecol[c(1,3)]) + 
  theme(text = element_text(size=14),legend.position='none') +
  xlab('Days post-infection') + ylab('Haematocrit (%)') + ggtitle('d')

##-- Create variable Overall Milk Volume
a = aggregate(PL ~ ID_ANIMAL,FUN = sum,data=goat[!is.na(goat$PL),])
ab = unique(merge(a,goat[,c('ID_ANIMAL','Group')],by='ID_ANIMAL'))
p6 = ggplot(ab,aes(x = Group,y = PL/1000, gp = Group, fill = Group)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(.9), binwidth = .4) + 
  scale_fill_manual(values = vecol) + 
  theme(text = element_text(size=14),legend.position = 'bottom') + 
  scale_y_continuous(limits = c(2,8),breaks = seq(2,8,2)) +
  xlab('Days post-infection') + ylab('Milk volume (L)') + ggtitle('f')

pdf(file='./Figure5.pdf',width=14,height=8)
multiplot(p1,p2,p3,p4,p5,p6,cols=2)
dev.off()

#======================================================================================
###== Model testing EWES
#======================================================================================
require(lme4)

infewe$Lot = factor(infewe$Lot)
infewe$Day = factor(infewe$Day)
infewe$ID_ANIMAL = factor(infewe$ID_ANIMAL)
infewe$ind = match(infewe$ID_ANIMAL,levels(infewe$ID_ANIMAL))

### FEC
# moyennes et ?cart-types FEC suivant le lot et les jours de pr?l?vement
a = aggregate(FEC~Jour_Prlv, FUN=mean, data=infewe)
a$std = aggregate(FEC~Jour_Prlv, FUN=sd, data=infewe)
a

# get rid of D0 as 0 count in each group
fecewe = infewe[infewe$Day!='D0',]
fecewe$Day = factor(fecewe$Day)

# Normality test
shapiro.test(log(50+fecewe$FEC))
# data:  log(50 + infewe$FEC)
# W = 0.77764, p-value = 9.704e-11

shapiro.test(sqrt(sqrt(fecewe$FEC)))
# data:  sqrt(sqrt(fecewe$FEC))
# W = 0.8096, p-value = 8.584e-10

## Model failed to converge
m0 = glmer.nb(FEC ~ Lot*Day + (1|ind), 
              data = fecewe)

## Use transformed FEC data instead
m0 = lme(sqrt(sqrt(FEC)) ~ Lot*Day,
         random =~ 1|ind, 
         data = fecewe)

summary(m0)
# Random effects:
#   Formula: ~1 | ind
#         (Intercept) Residual
# StdDev:    2.316413 1.779074

# Fixed effects: sqrt(sqrt(FEC)) ~ Lot * Day 
#                      Value Std.Error DF   t-value p-value
# (Intercept)       0.429404 0.8431525 66  0.509284  0.6123
# LotLupInf        -0.234375 1.1923977 22 -0.196558  0.8460
# DayD21            1.529762 0.7263040 66  2.106228  0.0390
# DayD24            3.519288 0.7263040 66  4.845475  0.0000
# DayD30            5.194100 0.7263040 66  7.151414  0.0000
# LotLupInf:DayD21 -0.098683 1.0271490 66 -0.096074  0.9238
# LotLupInf:DayD24  1.481139 1.0271490 66  1.441991  0.1540
# LotLupInf:DayD30  3.287710 1.0271490 66  3.200811  0.0021

#3.287710/sd(sqrt(sqrt(fecewe$FEC)))

Anova(m0, type = "III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: sqrt(sqrt(FEC))
#               Chisq Df Pr(>Chisq)    
# (Intercept)  0.2594  1   0.610553    
# Lot          0.0386  1   0.844174    
# Day         58.6661  3  1.133e-12 ***
# Lot:Day     14.3312  3   0.002487 ** 

# Check for normality of residuals / OK
shapiro.test(residuals(m0))
# W = 0.98002, p-value = 0.1509

# QQplot résidus
qqnorm(m0) # OK
plot(m0) # OK
rm(fecewe)

### Ht
htewe = infewe[,c('ind','Lot','Day','Ht','Group')]

vecht0 = htewe$Ht[htewe$Day=='D0']
htewe$Ht0 = vecht0[match(htewe$ind,htewe$ind[htewe$Day=='D0'])]
htewe = htewe[htewe$Day!='D0',]
htewe$Day = factor(htewe$Day)
htewe$Group = factor(htewe$Group)
# Normality test
shapiro.test(htewe$Ht)
# data:  htewe$Ht
# W = 0.96458, p-value = 0.01072

## Use transformed FEC data instead
m1 = lme(Ht ~ Lot*Day + Ht0,
         random =~ 1|ind, 
         data = htewe)

summary(m1)
# Random effects:
#   Formula: ~1 | ind
#         (Intercept) Residual
# StdDev:    3.192032 1.119868
# 
# Fixed effects: Ht ~ Lot * Day + Ht0 
#                      Value Std.Error DF   t-value p-value
# (Intercept)      17.317451  8.716620 66  1.986716  0.0511
# LotLupInf        -2.775650  1.696149 21 -1.636442  0.1166
# DayD21           -1.916667  0.457184 66 -4.192331  0.0001
# DayD24           -1.083333  0.457184 66 -2.369578  0.0207
# DayD30           -1.833333  0.457184 66 -4.010055  0.0002
# Ht0               0.503922  0.241162 21  2.089557  0.0490
# LotLupInf:DayD21  0.333333  0.646556 66  0.515552  0.6079
# LotLupInf:DayD24 -0.500000  0.646556 66 -0.773328  0.4421
# LotLupInf:DayD30 -1.750000  0.646556 66 -2.706649  0.0086

Anova(m1, type = "III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: Ht
#               Chisq Df Pr(>Chisq)    
# (Intercept)  3.9470  1   0.046954 *  
# Lot          2.6779  1   0.101747    
# Day         22.6589  3  4.756e-05 ***
# Ht0          4.3662  1   0.036658 *  
# Lot:Day     11.9857  3   0.007432 ** 

# Check for normality of residuals / OK
shapiro.test(residuals(m1))
# W = 0.98938, p-value = 0.6426

# QQplot résidus
qqnorm(m1) # OK
plot(m1) # OK
rm(htewe)

### Weight
adg = ewe[!is.na(ewe$GMQtot),]
adg$ID_ANIMAL = factor(adg$ID_ANIMAL)
adg$ind = match(adg$ID_ANIMAL,levels(adg$ID_ANIMAL))
adg$Group = factor(adg$Group,levels = c("Conc-Ninf","Conc-Inf","Lup-Ninf","Lup-Inf"))
shapiro.test(adg$GMQtot)
# data:  adg$GMQtot
# W = 0.84994, p-value = 2.166e-05

m3 = lm(GMQtot ~ Group ,
         data = adg)

summary.aov(m3)
#               Df Sum Sq Mean Sq F value Pr(>F)  
#   Group        3 142709   47570   3.724  0.018 *
#   Residuals   44 562053   12774  

summary(m3)
# Call:
#   lm(formula = GMQtot ~ Group, data = adg)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -281.825  -36.469    9.796   65.954  242.375 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     51.325     32.627   1.573   0.1229  
# GroupConc-Inf    4.308     46.141   0.093   0.9260  
# GroupLup-Ninf    7.058     46.141   0.153   0.8791  
# GroupLup-Inf  -122.000     46.141  -2.644   0.0113 *


### Prepatent period
shapiro.test(ppe$Tps_exc[ppe$Lot=='ConcInf'])
# data:  ppe$Tps_exc[ppe$Lot == "ConcInf"]
# W = 0.85737, p-value = 0.04532
shapiro.test(ppe$Tps_exc[ppe$Lot=='LupInf'])
# data:  ppe$Tps_exc[ppe$Lot == "LupInf"]
# W = 0.8294, p-value = 0.02063
var.test(Tps_exc ~ Lot,data=ppe)
# data:  Tps_exc by Lot
# F = 2.1707, num df = 11, denom df = 11, p-value = 0.2144
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.6249049 7.5404687
# sample estimates:
#   ratio of variances 
# 2.170732 
t.test(Tps_exc ~ Lot,data=ppe,var.equal = TRUE)
# data:  Tps_exc by Lot
# t = 1.0077, df = 22, p-value = 0.3246
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.587153  4.587153
# sample estimates:
#   mean in group ConcInf  mean in group LupInf 
# 24.75                 23.25  

#======================================================================================
###== Model testing GOATS
#======================================================================================

infgoat$Lot = factor(infgoat$Lot)
infgoat$Day = factor(infgoat$Day)
infgoat$ID_ANIMAL = factor(infgoat$ID_ANIMAL)
infgoat$ind = match(infgoat$ID_ANIMAL,levels(infgoat$ID_ANIMAL))

### FEC
# moyennes et ?cart-types FEC suivant le lot et les jours de pr?l?vement
a = aggregate(FEC~Lot+Jour_Prlv, FUN=mean, data=infgoat)
a$std = aggregate(FEC~Lot+Jour_Prlv, FUN=sd, data=infgoat)

# get rid of D0 as 0 count in each group
fecgoat = infgoat[infgoat$Day!='D0',]
fecgoat$Day = factor(fecgoat$Day)

# Normality test
shapiro.test(log(50+fecgoat$FEC))
# data:  log(50 + infgoat$FEC)
# W = 0.77044, p-value = 6.114e-11

shapiro.test(sqrt(sqrt(fecgoat$FEC)))
# data:  sqrt(sqrt(fecgoat$FEC))
# W = 0.8096, p-value = 8.589e-10

## Use transformed FEC data instead
m0g = lme(sqrt(sqrt(FEC)) ~ Lot*Day,
         random =~ 1|ind, 
         data = fecgoat)

summary(m0g)
# Random effects:
#   Formula: ~1 | ind
#         (Intercept) Residual
# StdDev:   0.8296296 1.219841
# 
# Fixed effects: sqrt(sqrt(FEC)) ~ Lot * Day 
#                      Value Std.Error DF   t-value p-value
# (Intercept)       0.359028 0.4258617 66  0.843063  0.4022
# LotLupInf        -0.359028 0.6022594 22 -0.596135  0.5572
# DayD21            0.184805 0.4979981 66  0.371096  0.7118
# DayD24            2.291319 0.4979981 66  4.601059  0.0000
# DayD30            4.974419 0.4979981 66  9.988830  0.0000
# LotLupInf:DayD21  0.804903 0.7042757 66  1.142880  0.2572
# LotLupInf:DayD24  1.006336 0.7042757 66  1.428894  0.1578
# LotLupInf:DayD30  0.990397 0.7042757 66  1.406264  0.1643

Anova(m0g, type = "III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: sqrt(sqrt(FEC))
#                 Chisq Df Pr(>Chisq)    
# (Intercept)   0.7108  1     0.3992    
# Lot           0.3554  1     0.5511    
# Day         130.2528  3     <2e-16 ***
# Lot:Day       2.7386  3     0.4337   

# Check for normality of residuals / OK
shapiro.test(residuals(m0))
# W = 0.98002, p-value = 0.1509

# QQplot résidus
qqnorm(m0) # OK
plot(m0) # OK
rm(fecgoat)

### Ht
htgoat = infgoat[,c('ind','Lot','Day','Ht','Group')]

vecht0 = htgoat$Ht[htgoat$Day=='D0']
htgoat$Ht0 = vecht0[match(htgoat$ind,htgoat$ind[htgoat$Day=='D0'])]
htgoat = htgoat[htgoat$Day!='D0',]
htgoat$Day = factor(htgoat$Day)
htgoat$Group = factor(htgoat$Group)
# Normality test
shapiro.test(htgoat$Ht)
# data:  htgoat$Ht
# W = 0.9806, p-value = 0.08092

## Use transformed FEC data instead
m1g = lme(Ht ~ Lot*Day + Ht0,
         random =~ 1|ind, 
         data = htgoat)

summary(m1g)
# Random effects:
#   Formula: ~1 | ind
# (Intercept) Residual
# StdDev:    1.467343 1.084207
# 
# Fixed effects: Ht ~ Lot * Day + Ht0 
#                      Value Std.Error DF   t-value p-value
# (Intercept)       5.532589  5.478098 66  1.009947  0.3162
# LotLupInf         0.300679  0.939176 21  0.320152  0.7520
# DayD21            0.000000  0.442626 66  0.000000  1.0000
# DayD24            0.166667  0.442626 66  0.376541  0.7077
# DayD30           -1.500000  0.442626 66 -3.388868  0.0012
# Ht0               0.831505  0.214533 21  3.875877  0.0009
# LotLupInf:DayD21 -0.916667  0.625967 66 -1.464400  0.1478
# LotLupInf:DayD24 -0.333333  0.625967 66 -0.532509  0.5962
# LotLupInf:DayD30  2.416667  0.625967 66  3.860692  0.0003

Anova(m1g, type = "III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: Ht
#               Chisq Df Pr(>Chisq)    
# (Intercept)  1.0200  1  0.3125206    
# Lot          0.1025  1  0.7488531    
# Day         18.7154  3  0.0003131 ***
# Ht0         15.0224  1  0.0001062 ***
# Lot:Day     32.9291  3  3.333e-07 *** 

# Check for normality of residuals / OK
shapiro.test(residuals(m1))
# W = 0.96066, p-value = 0.00143

# QQplot résidus
qqnorm(m1) # OK
plot(m1) # OK
rm(htgoat)

### Milk
shapiro.test(ab$PL)
#W = 0.97527, p-value = 0.4

ab$Group = factor(ab$Group,levels = c("Conc-Ninf","Conc-Inf","Lup-Ninf","Lup-Inf"))

m4 = lm(PL ~ Group ,
        data = ab)

summary.aov(m4)
#             Df   Sum Sq Mean Sq F value Pr(>F)
# Group        3  1015999  338666   0.299  0.826
# Residuals   44 49764754 1131017               

summary(m4)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1968.7  -689.7   155.7   760.1  2065.1 
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)     5064.9      307.0  16.498   <2e-16 ***
#  GroupConc-Inf   203.25     434.17   0.468    0.642    
#  GroupLup-Ninf  -207.67     434.17  -0.478    0.635    
#  GroupLup-Inf    -19.25     434.17  -0.044    0.965    

### Prepatent period
shapiro.test(ppg$Tps_exc[ppg$Lot=='ConcInf'])
# data:  ppg$Tps_exc[ppg$Lot == "ConcInf"]
# W = 0.84492, p-value = 0.0318
shapiro.test(ppg$Tps_exc[ppg$Lot=='LupInf'])
# data:  ppg$Tps_exc[ppg$Lot == "LupInf"]
# W = 0.74519, p-value = 0.002381
var.test(Tps_exc ~ Lot,data=ppg)
# data:  Tps_exc by Lot
# F = 1.7107, num df = 11, denom df = 11, p-value = 0.3868
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.4924788 5.9425375
# sample estimates:
#   ratio of variances 
# 1.710723 

t.test(Tps_exc ~ Lot,data=ppg,var.equal = TRUE)
# data:  Tps_exc by Lot
# t = 0.70417, df = 22, p-value = 0.4887
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.269303  4.602637
# sample estimates:
#   mean in group ConcInf  mean in group LupInf 
# 25.33333              24.16667 

