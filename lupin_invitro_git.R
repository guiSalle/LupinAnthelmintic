#Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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

#==13.04.2017================================================================
# LUPINE SEED ANTHELMINTIC EFFECT ANALYSIS
#============================================================================
require(lattice)
require(nlme)
require(ggplot2)
require(drc)
library(multcomp)
require(lsmeans)
require(lawstat)
draft='.'
### colors used: 
#a50026
#d73027
#f46d43
#fdae61
#fee08b
#ffffbf
#d9ef8b
#a6d96a
#66bd63
#1a9850
#006837

######====================================
#= FRACTION test
######====================================
setwd('./Fractions/')
theme_set(theme_bw())
f= read.csv(file='./LMIAOctobre2016.csv',header=T,sep=';')
f$mig = f$NoL3

#f = aggregate(mig ~ Strain + Fraction, data = frac,FUN = mean)
f$inhib=0
f$inhib[f$Strain=='KOK'] = 100*f$mig[f$Strain=='KOK']/mean(f$mig[f$Fraction=='TM' & f$Strain=='KOK'])
f$inhib[f$Strain=='WEY'] = 100*f$mig[f$Strain=='WEY']/mean(f$mig[f$Fraction=='TM' & f$Strain=='WEY'])

fres = aggregate(inhib ~ Strain + Fraction,data=f,FUN=mean)
fres$std = aggregate(inhib ~ Strain + Fraction,data=f,FUN=sd)[,3]
fres$souche = 'Isolat sensible'
fres$souche[fres$Strain=='KOK'] = 'Isolat multi-résistant'
fres$inhib[fres$inhib<0]=0
levels(fres$Fraction) = c('Alcaloides',"F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","LEV","TM")
fres


######====================================
#= VARIETY SCREENING
######====================================
setwd('./HFevrier/')

mig = read.csv("mig_tot.csv",sep=";")

##-- Screening based on 28.05.2014 results
mig = mig[mig$Date=='28-mai-14',]
mig$Mol = factor(mig$Mol)

levels(mig$Mol)
# [1] "CLOVIS_(EB)"  "E063_(EB)"    "E063_(FA)"    "EGY014_(EB)"  "EGY014_(FA)"  "EGY100_(EB)"  "EGY100_(FA)"  "ENERGY_(EB)" 
# [9] "LANG061_(EB)" "LANG061_(FA)" "LANG172_(EB)" "LL049_(EB)"   "LL151_(EB)"   "LL151_(FA)"   "LM261_(EB)"   "LM261_(FA)"  
# [17] "ORUS_(EB)"    "Temoin" 

extr = mig[mig$Mol!='Temoin',]
extr$Mol2 = unlist(strsplit(as.character(extr$Mol),'_'))[which(seq(1:(2*(dim(extr)[1])))%%2!=0)]
extr$ext =unlist(strsplit(as.character(extr$Mol),'_'))[which(seq(1:(2*(dim(extr)[1])))%%2==0)]
fa = extr[extr$ext=="(FA)",] #& extr$Souche=='Weybridge',]

extr = extr[extr$ext=="(EB)",] # & extr$Souche=='Weybridge',]

###---- Crude extracts
tmSH=mig[row.names(mig) %in% c(1456,1457,1458),] ###-- Control for Alkaloid-rich seeds
tmSL=mig[row.names(mig) %in% c(1495,1496,1497),] ###-- Control for Low-Alkaloid seeds
tmRH=mig[row.names(mig) %in% c(1435,1436,1437),] ###-- Control for Alkaloid-rich seeds
tmRL=mig[row.names(mig) %in% c(1477,1478,1479),] ###-- Control for Low-Alkaloid seeds

tmSH$Mol2='tmSH'
tmSH$ext='tm'
tmSL$Mol2='tmSL'
tmSL$ext='tm'
tmRH$Mol2='tmRH'
tmRH$ext='tm'
tmRL$Mol2='tmRL'
tmRL$ext='tm'

extr=rbind(extr,tmSH)
extr=rbind(extr,tmSL)
extr=rbind(extr,tmRH)
extr=rbind(extr,tmRL)

extr$Date=NULL
#extr$Souche=NULL
extr$Mol=extr$Mol2
extr$Mol2=NULL
extr = extr[,c(1,2,8,9)]
extr$group = 'Low Alkaloid Content'
extr$group[extr$Mol %in% c('E063','LANG061','PAC036','EGY014','EGY100','LM261','LL151')]='High Alkaloid Content'

extr$pmig = -1

extr$pmig[(extr$group=='High Alkaloid Content' & extr$Souche=='Weybridge')|extr$group=='tmSH'] = extr$Mig[(extr$group=='High Alkaloid Content' & extr$Souche=='Weybridge')|extr$group=='tmSH']/mean(extr$Mig[extr$Mol=='tmSH'])
extr$pmig[(extr$group=='Low Alkaloid Content'& extr$Souche=='Weybridge')|extr$group=='tmSL'] = extr$Mig[(extr$group=='Low Alkaloid Content'& extr$Souche=='Weybridge')|extr$group=='tmSL']/mean(extr$Mig[extr$Mol=='tmSL'])
extr$pmig[(extr$group=='High Alkaloid Content' & extr$Souche=='Kokstad')|extr$group=='tmRH'] = extr$Mig[(extr$group=='High Alkaloid Content' & extr$Souche=='Kokstad')|extr$group=='tmRH']/mean(extr$Mig[extr$Mol=='tmRH'])
extr$pmig[(extr$group=='Low Alkaloid Content'& extr$Souche=='Kokstad')|extr$group=='tmRL'] = extr$Mig[(extr$group=='Low Alkaloid Content'& extr$Souche=='Kokstad')|extr$group=='tmRL']/mean(extr$Mig[extr$Mol=='tmRL'])

extr$Mol[extr$Mol=='tmSH']='Control'
extr$Mol[extr$Mol=='tmSL']='Control'
extr$Mol[extr$Mol=='tmRH']='Control'
extr$Mol[extr$Mol=='tmRL']='Control'

extr$Mol = factor(extr$Mol)

#-- Modeling: Mean migration of controls as covariates
#df = extr[extr$Mol !='Control',]
df = extr
df$Mol=relevel(df$Mol, ref="Control") 
df$Mol = factor(df$Mol)
df$control = 0
df$control[df$Souche=='Weybridge'] = mean(df$pmig[df$Souche=='Weybridge'])
df$control[df$Souche=='Kokstad'] = mean(df$pmig[df$Souche=='Kokstad'])
df$pinhib = 1 - df$pmig
df$pinhib[df$pinhib<0] = 0

#df = df[row.names(df)!=1446,] ##-- discard one LL151 replicate

ggplot(df,aes(Mig,fill=Mol)) +
  geom_density(alpha=.1) +
  facet_wrap(~ Mol)

# bartlett.test(pinhib ~ Mol,data=df) ###- Homogeneous if no control
# bartlett.test(pinhib ~ Souche,data=df) ###- Homogeneous var for each strain
###df = df[df$Mol!='Control',]

a = lm(pinhib ~ Mol*Souche , data = df)
anova(a)
summary(a)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                 0.05282    0.03449   1.531 0.131530    
#   MolCLOVIS                   0.56690    0.05974   9.489 4.23e-13 ***
#   MolE063                     0.72820    0.05974  12.189  < 2e-16 ***
#   MolEGY014                   0.79390    0.05974  13.289  < 2e-16 ***
#   MolEGY100                   0.76470    0.05974  12.800  < 2e-16 ***
#   MolENERGY                   0.70070    0.05974  11.729  < 2e-16 ***
#   MolLANG061                  0.65521    0.05974  10.968 2.35e-15 ***
#   MolLANG172                  0.54577    0.05974   9.136 1.52e-12 ***
#   MolLL049                    0.67958    0.05974  11.375 5.88e-16 ***
#   MolLL151                    0.58222    0.05974   9.746 1.69e-13 ***
#   MolLM261                    0.78660    0.05974  13.167  < 2e-16 ***
#   MolORUS                     0.53169    0.05974   8.900 3.58e-12 ***
#   SoucheWeybridge            -0.02759    0.04878  -0.566 0.574026    
#   MolCLOVIS:SoucheWeybridge  -0.25268    0.08449  -2.991 0.004185 ** 
#   MolE063:SoucheWeybridge    -0.22190    0.08449  -2.627 0.011202 *  
#   MolEGY014:SoucheWeybridge   0.14033    0.08449   1.661 0.102508    
#   MolEGY100:SoucheWeybridge  -0.12777    0.08449  -1.512 0.136287    
#   MolENERGY:SoucheWeybridge  -0.32226    0.08449  -3.814 0.000353 ***
#   MolLANG061:SoucheWeybridge  0.01325    0.08449   0.157 0.875944    
#   MolLANG172:SoucheWeybridge -0.48691    0.08449  -5.763 4.08e-07 ***
#   MolLL049:SoucheWeybridge   -0.25068    0.08449  -2.967 0.004472 ** 
#   MolLL151:SoucheWeybridge   -0.09844    0.08449  -1.165 0.249074    
#   MolLM261:SoucheWeybridge    0.11610    0.08449   1.374 0.175061    
#   MolORUS:SoucheWeybridge    -0.27710    0.08449  -3.280 0.001822 **

###--- Figure 1 : co-efficacy on Susceptible and Resistant strains
dfres = df[df$Souche=='Kokstad',]
dfsus = df[df$Souche=='Weybridge',]
tabsus = aggregate(pinhib ~ Mol ,data=dfsus,FUN=mean)
tabsus$stds = aggregate(pinhib ~ Mol ,data=dfsus,FUN=sd)[,2]
tabres = aggregate(pinhib ~ Mol ,data=dfres,FUN=mean)
tabres$stdr = aggregate(pinhib ~ Mol ,data=dfres,FUN=sd)[,2]
tab = merge(tabsus,tabres,by='Mol')
tab = tab[tab$Mol!='Control',]
tab = merge(unique(df[,c('Mol','group')]),tab,by='Mol')
tab = tab[order(tab$group,tab$Mol),]
tab$Mol = factor(tab$Mol, levels=tab$Mol)
tab$com = 'no'
tab$com[tab$Mol %in% c('CLOVIS','ENERGY','ORUS')]='yes'

pdf(file=paste0(draft,'Figure1.pdf'),width=14,height=8)
ggplot(tab,aes(x=pinhib.x,y=pinhib.y,col=com)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=tab$pinhib.y-tab$stdr,ymax=tab$pinhib.y + tab$stdr,width=0.02)) +
  geom_errorbarh(aes(xmin=tab$pinhib.x-tab$stds,xmax=tab$pinhib.x + tab$stds,height=0.02)) +
  ggrepel::geom_text_repel(segment.alpha = 0.2,size=5,box.padding = unit(1.5, "lines"),aes(x=pinhib.x,y=pinhib.y,label=Mol)) +
  scale_color_manual(values=c('black','forestgreen')) +
  facet_wrap(~ group) +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1)) +
  scale_x_continuous(breaks = c(0,0.2,0.4,.6,.8,1),limits = c(-0.02,1)) +
  ylab('% inhibition for the resistant isolate') +
  xlab('% inhibition for the susceptible isolate') +
  ggtitle('') + 
  theme(axis.title = element_text(size=18),legend.position = 'none',
        strip.text.x = element_text(size=18),
        axis.text.x = element_text(size=16,hjust = 1),
        axis.text.y = element_text(size=16,hjust = 1))
dev.off()

####### Suppl. TABLE 3
stab1 = aggregate(pinhib ~ Mol + Souche, data = df,FUN =mean)
stab1$std = aggregate(pinhib ~ Mol + Souche, data = df,FUN =sd)[,3]
stab1x = unique(merge(stab1,df[,c('Mol','group')],by='Mol'))

write.csv(stab1x,file=paste0(draft,'supplementary_table3.csv'),row.names = F,quote=F)

##-- Test low vs. high alkaloidic content
b = lm(pinhib ~ group , data=df)
summary(b)
# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.77675    0.03244  23.945  < 2e-16 ***
# groupLow Alkaloid Content -0.27114    0.04392  -6.173 5.11e-08 ***
lsmeans(b,'group')
# group                    lsmean         SE df  lower.CL  upper.CL
# High Alkaloid Content 0.7767476 0.03243935 64 0.7119425 0.8415526
# Low Alkaloid Content  0.5056042 0.02961294 64 0.4464455 0.5647628

##-- Test low vs. high alkaloidic content
dft = df[df$Mol!='Control' & df$group=='Low Alkaloid Content',]
dft$gp = 'other'
dft$gp[dft$Mol=='ENERGY'|dft$Mol=='LL049'] = 'ENLL'
b = lm(pinhib ~ gp + Souche , data=dft)
summary(b)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  0.58593    0.08658   6.768  2.4e-09 ***
#   gpother     -0.05985    0.09412  -0.636    0.527    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###---- Alkaloidic Fraction

fa$Souche=factor(fa$Souche)

tmS=mig[row.names(mig) %in% c(1414,1415,1416),] ###-- Control for susceptible isolate
tmR=mig[row.names(mig) %in% c(1393,1394,1395),] ###-- Control for resistant isolate
tmS$Mol2='tmS'
tmS$ext='tm'
tmR$Mol2='tmR'
tmR$ext='tm'

fa=rbind(fa,tmS)
fa=rbind(fa,tmR)

fa$Date=NULL
#fa$Souche=NULL
fa$Mol=fa$Mol2
fa$Mol2=NULL
fa = fa[,c(1,2,9)]
fa$pinhib=-1
fa$pinhib[(fa$Souche=='Weybridge')|fa$Mol=='tmS'] = 1 - fa$Mig[(fa$Souche=='Weybridge')|fa$Mol=='tmS']/mean(fa$Mig[fa$Mol=='tmS'])
fa$pinhib[(fa$Souche=='Kokstad')|fa$Mol=='tmR'] = 1 - fa$Mig[(fa$Souche=='Kokstad')|fa$Mol=='tmR']/mean(fa$Mig[fa$Mol=='tmR'])
fa$Mol[fa$Mol=='tmS']='Control'
fa$Mol[fa$Mol=='tmR']='Control'
fa$Mol = factor(fa$Mol)

#df2 = fa[fa$Mol !='Control',]
df2 = fa
df2$Mol = factor(df2$Mol)
df2$control = 0
df2$control[df2$Souche=='Weybridge'] = mean(df2$pinhib[df2$Souche=='Weybridge'])
df2$control[df2$Souche=='Kokstad'] = mean(df2$pinhib[df2$Souche=='Kokstad'])

# ggplot(df2,aes(pinhib,fill=Mol)) +
#   geom_density(alpha=.1) +
#   facet_wrap(~ Mol)

bartlett.test(pinhib ~ Mol,data=df2) ###- Homogeneous variance for Mol
bartlett.test(pinhib ~ Souche,data=df2) ###- Homogeneous var for each strain

p1B = ggplot(tab2,aes(x=Mol,y=lsmean,fill=Mol)) +
  geom_errorbar(aes(ymin=(tab2$lower.CL),ymax=(tab2$upper.CL),width=0.1)) +
  geom_bar(stat='identity') +
  theme_classic() +
  scale_fill_manual(values=c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf')) +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1)) +
  ylab('% inhibition relative to control') +
  xlab('Lupine breed') +
  ggtitle('B') +
  theme(axis.title = element_text(size=16),legend.position = 'none',
        strip.text.x = element_text(size=14),
        axis.text.x = element_text(size=14,angle = 45,hjust = 1),
        axis.text.y = element_text(size=14,hjust = 1))

####### Suppl. TABLE 4
stab2 = aggregate(pinhib ~ Mol + Souche, data = df2,FUN =mean)
stab2$std = aggregate(pinhib ~ Mol + Souche, data = df2,FUN =sd)[,3]

write.csv(stab2,file=paste0(draft,'supplementary_table4.csv'),row.names = F,quote=F)


#####--- Keep ENERGY = best commercial breed // E063 based on alkaloid fraction effect


######====================================
#= LUPINE BREED x ISOLATE effect
######====================================

setwd('./ODubois/res/')

#--- Metafile including every LMIA tests
lmia1=read.csv(file="./LMIA_M_EA_FNA/test de migration - e063 et ener (macer-fa-fna) sur wey et kok 120615.csv",header=T,sep=";")
lmia2=read.csv(file="./LMIA_M_EA_FNA/test de migration - e063 et ener (macer-fa-fna) sur wey et kok 160615.csv",header=T,sep=";")
lmia3=read.csv(file="./LMIA_M_EA_FNA/Test de migration (WEY, KOK, TELA) FA, FAm, FNA, macération RASSEMBLEMENT DES DEUX MANIPS.csv",header=T,sep=";")
lmia2$plq=lmia2$plq+4
lmiat=rbind(lmia3,lmia2)

options(digits=3)

# #----------------------------------------------------------------------------------------------------------------------------
# #--- Data description & QC check
# #----------------------------------------------------------------------------------------------------------------------------
# par(mfrow=c(2,2))
# 
# #---- Control larvae migration level (based on 100 Larvae)
# 
# tmlmia1=lmia1[lmia1$var=='tm' & lmia1$frac=='neg',]   #-- XP 12/06/2015
# aggregate(l3 ~ souche + frac, data=tmlmia1,FUN=mean)
# # souche frac       l3
# # 1    kok  neg 24.00000
# # 2    wey  neg 45.66667
# boxplot(l3 ~ souche, data=tmlmia1,ylim=c(0,100),main='Neg Control Migration - 12.06.2015')
# lmia1[lmia1$souche=="wey" & lmia1$var=='tm' & lmia1$frac=='neg',]
# 
# tmlmia2=lmia2[lmia2$var=='tm' & lmia2$frac=='neg',]   #-- XP 16/06/2015
# aggregate(l3 ~ souche + frac, data=tmlmia2,FUN=mean)
# # souche frac       l3
# # 1    kok  neg 38.83333
# # 2    wey  neg 51.83333
# boxplot(l3 ~ souche, data=tmlmia2,ylim=c(0,100),main='Neg Control Migration - 16.06.2015')
# lmia2[lmia2$souche=="wey" & lmia2$var=='tm' & lmia2$frac=='neg',]
# 
# tmlmia3=lmia3[lmia3$var=='tm'& lmia3$frac=='neg',]   #-- XP Août 2015
# aggregate(l3 ~ souche + frac, data=tmlmia3,FUN=mean)
# # souche frac    l3
# # 1    kok  neg 43.75
# # 2   tela  neg 66.00
# # 3    wey  neg 59.75
# boxplot(l3 ~ souche, data=tmlmia3,ylim=c(0,100),main='Neg Control Migration - 08.2015')
# lmia3[lmia3$souche=="wey" & lmia3$var=='tm' & lmia3$frac=='neg',]
# 
# #--- Density plots
# plot(density(tmlmia1$l3))
# plot(density(tmlmia2$l3))
# plot(density(tmlmia3$l3))

#--------- consider the only LMIA 3 results as not enough larvae migrated in other XP -----------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#--- ANALYSIS on LMIA 3 results only
#----------------------------------------------------------------------------------------------------------------------------

#--- Compute % of migration

lmia3$pmig=NA

mcont.1=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==1 & lmia3$souche=="wey")])
mcont.2=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==2 & lmia3$souche=="wey")])
mcont.3=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==3 & lmia3$souche=="wey")])
mcont.4=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==4 & lmia3$souche=="wey")])

mcont.5=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==1 & lmia3$souche=="kok")])
mcont.6=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==2 & lmia3$souche=="kok")])		
mcont.7=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==3 & lmia3$souche=="kok")])
mcont.8=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==4 & lmia3$souche=="kok")])		

mcont.9=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==1 & lmia3$souche=="tela")])
mcont.10=mean(lmia3$l3[which(lmia3$frac=="neg" & lmia3$plq==2 & lmia3$souche=="tela")])

lmia3$pmig[which(lmia3$plq=="1" & lmia3$souche=="wey")]= round(lmia3$l3[lmia3$plq=="1" & lmia3$souche=="wey"]/mcont.1,digits=2)
lmia3$pmig[which(lmia3$plq=="2" & lmia3$souche=="wey")]= round(lmia3$l3[lmia3$plq=="2" & lmia3$souche=="wey"]/mcont.2,digits=2)
lmia3$pmig[which(lmia3$plq=="3" & lmia3$souche=="wey")]= round(lmia3$l3[lmia3$plq=="3" & lmia3$souche=="wey"]/mcont.3,digits=2)
lmia3$pmig[which(lmia3$plq=="4" & lmia3$souche=="wey")]= round(lmia3$l3[lmia3$plq=="4" & lmia3$souche=="wey"]/mcont.4,digits=2)

lmia3$pmig[which(lmia3$plq=="1" & lmia3$souche=="kok")]= round(lmia3$l3[lmia3$plq=="1" & lmia3$souche=="kok"]/mcont.5,digits=2)
lmia3$pmig[which(lmia3$plq=="2" & lmia3$souche=="kok")]= round(lmia3$l3[lmia3$plq=="2" & lmia3$souche=="kok"]/mcont.6,digits=2)
lmia3$pmig[which(lmia3$plq=="3" & lmia3$souche=="kok")]= round(lmia3$l3[lmia3$plq=="3" & lmia3$souche=="kok"]/mcont.7,digits=2)
lmia3$pmig[which(lmia3$plq=="4" & lmia3$souche=="kok")]= round(lmia3$l3[lmia3$plq=="4" & lmia3$souche=="kok"]/mcont.8,digits=2)

lmia3$pmig[which(lmia3$plq=="1" & lmia3$souche=="tela")]= round(lmia3$l3[lmia3$plq=="1" & lmia3$souche=="tela"]/mcont.9,digits=2)
lmia3$pmig[which(lmia3$plq=="2" & lmia3$souche=="tela")]= round(lmia3$l3[lmia3$plq=="2" & lmia3$souche=="tela"]/mcont.10,digits=2)

####--- QC CHECK: % migration KOKSTAD / RESISTANT ON PLATE 3 AND 4 ONLY
lmia3[lmia3$souche=="wey" & lmia3$var=="tm" & (lmia3$frac=='neg'),]
#     souche var frac plq rep l3 pmig
# 28     wey  tm  neg   1   1 71 1.02
# 29     wey  tm  neg   1   2 66 0.95
# 30     wey  tm  neg   1   3 72 1.03
# 118    wey  tm  neg   2   1 71 1.08
# 119    wey  tm  neg   2   2 62 0.94
# 120    wey  tm  neg   2   3 64 0.97
# 205    wey  tm  neg   3   1 47 0.99
# 206    wey  tm  neg   3   2 49 1.03
# 207    wey  tm  neg   3   3 47 0.99
# 259    wey  tm  neg   4   1 59 1.05
# 260    wey  tm  neg   4   2 52 0.93
# 261    wey  tm  neg   4   3 57 1.02

lmia3[lmia3$souche=="kok" & lmia3$var=="tm" & (lmia3$frac=='pos' | lmia3$frac=='lev'),]
#     souche var frac plq rep l3 pmig
# 55     kok  tm  pos   1   1  0 0.00
# 56     kok  tm  pos   1   2  0 0.00
# 57     kok  tm  pos   1   3  0 0.00
# 145    kok  tm  pos   2   1  2 0.05
# 146    kok  tm  pos   2   2  0 0.00
# 147    kok  tm  pos   2   3  0 0.00
# 226    kok  tm  lev   3   1 31 0.78
# 227    kok  tm  lev   3   2 29 0.73
# 228    kok  tm  lev   3   3 28 0.71
# 280    kok  tm  lev   4   1 27 0.71
# 281    kok  tm  lev   4   2 23 0.61
# 282    kok  tm  lev   4   3 24 0.63

####################--------------------------------------------------------------------
##---- KEEP PLATES 1 AND 2 FOR WEYBRIDGE // PLATES 3 AND 4 FOR KOKSTAD
##-- Split dataset by strain 
##-- Do not include control wells
####################-------------------------------------------------------------------------

SUS = lmia3[lmia3$souche=='wey' & lmia3$plq==1 |
              lmia3$souche=='wey' & lmia3$plq==2,]
SUS$frac[SUS$frac=='pos']='lev'
RES = lmia3[lmia3$souche=='kok' & lmia3$plq==3 |
              lmia3$souche=='kok' & lmia3$plq==4,]
RES=RES[RES$frac!="fam" & RES$frac !="nico",]
SUS=SUS[SUS$frac!="fam" & SUS$frac !="nico",]

RES$mol=factor(paste(RES$var,RES$frac,sep='.'))
SUS$mol=factor(paste(SUS$var,SUS$frac,sep='.'))

RES$solution='bla'
RES$solution[RES$mol=='e063.fa'] ='E063.Alkaloids'
RES$solution[RES$mol=='ener.fa'] ='ENERGY.Alkaloids'
RES$solution[RES$mol=='e063.fna'] ='E063.Non-alkaloids'
RES$solution[RES$mol=='ener.fna'] ='ENERGY.Non-alkaloids'
RES$solution[RES$mol=='e063.macer'] ='E063-Tot. Extract'
RES$solution[RES$mol=='ener.macer'] ='ENERGY-Tot. Extract'
RES$solution[RES$mol=='tm.lev'] ='Levamisole'
RES$solution[RES$mol=='tm.neg'] ='Water'
SUS$solution='bla'
SUS$solution[SUS$mol=='e063.fa'] ='E063.Alkaloids'
SUS$solution[SUS$mol=='ener.fa'] ='ENERGY.Alkaloids'
SUS$solution[SUS$mol=='e063.fna'] ='E063.Non-alkaloids'
SUS$solution[SUS$mol=='ener.fna'] ='ENERGY.Non-alkaloids'
SUS$solution[SUS$mol=='e063.macer'] ='E063-Tot. Extract'
SUS$solution[SUS$mol=='ener.macer'] ='ENERGY-Tot. Extract'
SUS$solution[SUS$mol=='tm.lev'] ='Levamisole'
SUS$solution[SUS$mol=='tm.neg'] ='Water'

SUS$solution=factor(SUS$solution)
RES$solution=factor(RES$solution)

RES$isolate = 'Resistant H.contortus'
SUS$isolate = 'Susceptible H.contortus'

#--- DATA DESCRIPTION
table(RES$mol)
table(SUS$mol)

########---------------------- Test and graph original data
# require(fitdistrplus)
# plot(fitdist(data = RES$l3, distr = "pois", method = "mle"))
# plot(fitdist(data = RES$l3, distr = "nbinom", method = "mle"))
# plot(fitdist(data = RES$l3, distr = "norm", method = "mle"))
# 
# qqcomp(fitdist(data = RES$l3, distr = "pois", method = "mle"),main="QQ-Plot - Poisson")
# qqcomp(fitdist(data = RES$l3, distr = "nbinom", method = "mle"),main="QQ-Plot - Nbinom")
# qqcomp(fitdist(data = RES$l3, distr = "norm", method = "mle"),main="QQ-Plot - Gaussian")


#########--------------------- PLOT INHIBITION FRACTION
tmp=rbind(RES,SUS)
tmp$isolate=factor(tmp$isolate)
# 
# ggplot(tmp,aes(x=solution,y=1-pmig,fill=solution))+
#   geom_boxplot() +
#   facet_wrap(~ isolate,ncol=1) +
#   theme_bw() +
#   ylab('% inhibition') +
#   xlab('Solution') +
#   theme(legend.pos='none',
#         axis.title = element_text(size=14),
#         strip.text.x = element_text(size=14),
#         axis.text = element_text(size=12,angle = 45,hjust = 1))

#########--------------------- MODELING MIGRATION - SUSCEPTIBLE ISOLATE
SUS$solution = factor(SUS$solution,levels=c('Water', 'Levamisole',
                                            'E063-Tot. Extract','E063.Alkaloids','E063.Non-alkaloids',
                                            'ENERGY-Tot. Extract','ENERGY.Alkaloids','ENERGY.Non-alkaloids'))

m = lme(1-pmig ~ solution ,
        random=~1|plq, data=SUS)

summary(m)
# Random effects:
#   Formula: ~1 | plq
# (Intercept) Residual
# StdDev:      0.0398    0.101
# 
# Fixed effects: pmig ~ solution 
#                               Value Std.Error DF t-value p-value
# (Intercept)                   0.998    0.0498 39   20.06  0.0000
# solutionLevamisole            0.993    0.0581 39   17.11  0.0000
# solutionE063-Aqu. Extract     0.828    0.0581 39   14.27  0.0000
# solutionE063.Alkaloids        0.410    0.0581 39    7.06  0.0000
# solutionE063.Non-alkaloids    0.375    0.0581 39    6.46  0.0000
# solutionENERGY-Aqu. Extract   0.425    0.0581 39    7.32  0.0000
# solutionENERGY.Alkaloids      0.413    0.0581 39    7.12  0.0000
# solutionENERGY.Non-alkaloids -0.125    0.0581 39   -2.15  0.0376

anova.lme(m)
#             numDF denDF F-value p-value
# (Intercept)     1    39   173.1  <.0001
# solution        7    39    82.3  <.0001
summary(glht(m,mcp(solution="Tukey")))
#                                                   Estimate Std. Error z value Pr(>|z|)    
#   Levamisole - Water == 0                          0.99333    0.05805   17.11   <0.001 ***
#   E063-Aqu. Extract - Water == 0                   0.82833    0.05805   14.27   <0.001 ***
#   E063.Alkaloids - Water == 0                      0.41000    0.05805    7.06   <0.001 ***
#   E063.Non-alkaloids - Water == 0                  0.37500    0.05805    6.46   <0.001 ***
#   ENERGY-Aqu. Extract - Water == 0                 0.42500    0.05805    7.32   <0.001 ***
#   ENERGY.Alkaloids - Water == 0                    0.41333    0.05805    7.12   <0.001 ***
#   ENERGY.Non-alkaloids - Water == 0               -0.12500    0.05805   -2.15    0.381    
#   E063-Aqu. Extract - Levamisole == 0             -0.16500    0.05805   -2.84    0.085 .  
#   E063.Alkaloids - Levamisole == 0                -0.58333    0.05805  -10.05   <0.001 ***
#   E063.Non-alkaloids - Levamisole == 0            -0.61833    0.05805  -10.65   <0.001 ***
#   ENERGY-Aqu. Extract - Levamisole == 0           -0.56833    0.05805   -9.79   <0.001 ***
#   ENERGY.Alkaloids - Levamisole == 0              -0.58000    0.05805   -9.99   <0.001 ***
#   ENERGY.Non-alkaloids - Levamisole == 0          -1.11833    0.05805  -19.26   <0.001 ***
#   E063.Alkaloids - E063-Aqu. Extract == 0         -0.41833    0.05805   -7.21   <0.001 ***
#   E063.Non-alkaloids - E063-Aqu. Extract == 0     -0.45333    0.05805   -7.81   <0.001 ***
#   ENERGY-Aqu. Extract - E063-Aqu. Extract == 0    -0.40333    0.05805   -6.95   <0.001 ***
#   ENERGY.Alkaloids - E063-Aqu. Extract == 0       -0.41500    0.05805   -7.15   <0.001 ***
#   ENERGY.Non-alkaloids - E063-Aqu. Extract == 0   -0.95333    0.05805  -16.42   <0.001 ***
#   E063.Non-alkaloids - E063.Alkaloids == 0        -0.03500    0.05805   -0.60    0.999    
#   ENERGY-Aqu. Extract - E063.Alkaloids == 0        0.01500    0.05805    0.26    1.000    
#   ENERGY.Alkaloids - E063.Alkaloids == 0           0.00333    0.05805    0.06    1.000    
#   ENERGY.Non-alkaloids - E063.Alkaloids == 0      -0.53500    0.05805   -9.22   <0.001 ***
#   ENERGY-Aqu. Extract - E063.Non-alkaloids == 0    0.05000    0.05805    0.86    0.989    
#   ENERGY.Alkaloids - E063.Non-alkaloids == 0       0.03833    0.05805    0.66    0.998    
#   ENERGY.Non-alkaloids - E063.Non-alkaloids == 0  -0.50000    0.05805   -8.61   <0.001 ***
#   ENERGY.Alkaloids - ENERGY-Aqu. Extract == 0     -0.01167    0.05805   -0.20    1.000    
#   ENERGY.Non-alkaloids - ENERGY-Aqu. Extract == 0 -0.55000    0.05805   -9.47   <0.001 ***
#   ENERGY.Non-alkaloids - ENERGY.Alkaloids == 0    -0.53833    0.05805   -9.27   <0.001 ***

a=summary(lsmeans(m,'solution'))[,c('lsmean','SE','lower.CL','upper.CL')]
# solution             lsmean     SE df lower.CL upper.CL
# Water                 0.998 0.0498  1  0.36591    1.631
# Levamisole            0.005 0.0498  1 -0.62742    0.637
# E063-Aqu. Extract     0.170 0.0498  1 -0.46242    0.802
# E063.Alkaloids        0.588 0.0498  1 -0.04409    1.221
# E063.Non-alkaloids    0.623 0.0498  1 -0.00909    1.256
# ENERGY-Aqu. Extract   0.573 0.0498  1 -0.05909    1.206
# ENERGY.Alkaloids      0.585 0.0498  1 -0.04742    1.217
# ENERGY.Non-alkaloids  1.123 0.0498  1  0.49091    1.756
a$solution=levels(SUS$solution)
a$std = aggregate((1-pmig) ~ solution, data=SUS,FUN=sd)[,2]
a
p1a = ggplot(a,aes(x=solution,y=lsmean,fill=solution))+
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('#a50026','#d73027','#f46d43',
                             '#66bd63','#1a9850','#006837',
                             '#542788','#92c5de'))+
  theme_classic() +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(-0.31,1)) +
  geom_errorbar(aes(ymin=(a$lsmean-a$SE),ymax=(a$lsmean + a$SE),width=0.1)) +
  ylab('% inhibition relative to control') +
  xlab('') +
  ggtitle('a')+
  theme(legend.pos='none',
        plot.title = element_text(size=20),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45,hjust = 1))

#########--------------------- MODELING MIGRATION - RESISTANT ISOLATE
RES$solution = factor(RES$solution,levels=c('Water', 'Levamisole',
                                            'E063-Tot. Extract','E063.Alkaloids','E063.Non-alkaloids',
                                            'ENERGY-Tot. Extract','ENERGY.Alkaloids','ENERGY.Non-alkaloids'))

m1 = lme(1-pmig ~ solution ,
        random=~1|plq, data=RES)

summary(m1)

# Fixed effects: pmig ~ solution 
#                               Value Std.Error DF t-value p-value
# (Intercept)                   1.002    0.0462 39   21.67  0.0000
# solutionLevamisole            0.307    0.0654 39    4.69  0.0000
# solutionE063-Aqu. Extract     0.285    0.0654 39    4.36  0.0001
# solutionE063.Alkaloids        0.403    0.0654 39    6.17  0.0000
# solutionE063.Non-alkaloids    0.127    0.0654 39    1.94  0.0599
# solutionENERGY-Aqu. Extract   0.488    0.0654 39    7.47  0.0000
# solutionENERGY.Alkaloids      0.550    0.0654 39    8.41  0.0000
# solutionENERGY.Non-alkaloids -0.098    0.0654 39   -1.50  0.1406

summary(glht(m1,mcp(solution="Tukey")))
# Estimate Std. Error z value Pr(>|z|)    
# Levamisole - Water == 0                           0.3067     0.0654    4.69   <0.001 ***
# E063-Aqu. Extract - Water == 0                    0.2850     0.0654    4.36   <0.001 ***
# E063.Alkaloids - Water == 0                       0.4033     0.0654    6.17   <0.001 ***
# E063.Non-alkaloids - Water == 0                   0.1267     0.0654    1.94   0.5247    
# ENERGY-Aqu. Extract - Water == 0                  0.4883     0.0654    7.47   <0.001 ***
# ENERGY.Alkaloids - Water == 0                     0.5500     0.0654    8.41   <0.001 ***
# ENERGY.Non-alkaloids - Water == 0                -0.0983     0.0654   -1.50   0.8056    
# E063-Aqu. Extract - Levamisole == 0              -0.0217     0.0654   -0.33   1.0000    
# E063.Alkaloids - Levamisole == 0                  0.0967     0.0654    1.48   0.8190    
# E063.Non-alkaloids - Levamisole == 0             -0.1800     0.0654   -2.75   0.1071    
# ENERGY-Aqu. Extract - Levamisole == 0             0.1817     0.0654    2.78   0.1001    
# ENERGY.Alkaloids - Levamisole == 0                0.2433     0.0654    3.72   0.0051 ** 
# ENERGY.Non-alkaloids - Levamisole == 0           -0.4050     0.0654   -6.20   <0.001 ***
# E063.Alkaloids - E063-Aqu. Extract == 0           0.1183     0.0654    1.81   0.6131    
# E063.Non-alkaloids - E063-Aqu. Extract == 0      -0.1583     0.0654   -2.42   0.2305    
# ENERGY-Aqu. Extract - E063-Aqu. Extract == 0      0.2033     0.0654    3.11   0.0394 *  
# ENERGY.Alkaloids - E063-Aqu. Extract == 0         0.2650     0.0654    4.05   0.0013 ** 
# ENERGY.Non-alkaloids - E063-Aqu. Extract == 0    -0.3833     0.0654   -5.86   <0.001 ***
# E063.Non-alkaloids - E063.Alkaloids == 0         -0.2767     0.0654   -4.23   <0.001 ***
# ENERGY-Aqu. Extract - E063.Alkaloids == 0         0.0850     0.0654    1.30   0.8991    
# ENERGY.Alkaloids - E063.Alkaloids == 0            0.1467     0.0654    2.24   0.3255    
# ENERGY.Non-alkaloids - E063.Alkaloids == 0       -0.5017     0.0654   -7.67   <0.001 ***
# ENERGY-Aqu. Extract - E063.Non-alkaloids == 0     0.3617     0.0654    5.53   <0.001 ***
# ENERGY.Alkaloids - E063.Non-alkaloids == 0        0.4233     0.0654    6.48   <0.001 ***
# ENERGY.Non-alkaloids - E063.Non-alkaloids == 0   -0.2250     0.0654   -3.44   0.0137 *  
# ENERGY.Alkaloids - ENERGY-Aqu. Extract == 0       0.0617     0.0654    0.94   0.9818    
# ENERGY.Non-alkaloids - ENERGY-Aqu. Extract == 0  -0.5867     0.0654   -8.97   <0.001 ***
# ENERGY.Non-alkaloids - ENERGY.Alkaloids == 0     -0.6483     0.0654   -9.92   <0.001 ***
#   
b=summary(lsmeans(m1,'solution'))[,c('lsmean','SE','lower.CL','upper.CL')]
b$solution=levels(RES$solution)
b$std = aggregate(1-pmig ~ solution, data=RES,FUN=sd)[,2]
b
# lsmean     SE lower.CL upper.CL solution                std
# -0.00167 0.0462   -0.589    0.586 Water                0.0842
# 0.30500 0.0462   -0.282    0.892 Levamisole           0.0638
# 0.28333 0.0462   -0.304    0.871 E063-Aqu. Extract    0.0728
# 0.40167 0.0462   -0.186    0.989 E063.Alkaloids       0.0906
# 0.12500 0.0462   -0.462    0.712 E063.Non-alkaloids   0.0641
# 0.48667 0.0462   -0.101    1.074 ENERGY-Aqu. Extract  0.1555
# 0.54833 0.0462   -0.039    1.136 ENERGY.Alkaloids     0.0906
# -0.10000 0.0462   -0.687    0.487 ENERGY.Non-alkaloids 0.2034

p1b = ggplot(b,aes(x=solution,y=lsmean,fill=solution))+
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('#a50026','#d73027','#f46d43',
                             '#66bd63','#1a9850','#006837',
                             '#542788','#92c5de'))+
  theme_classic() +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(-0.31,1)) +
  geom_errorbar(aes(ymin=(b$lsmean-b$SE),ymax=(b$lsmean + b$SE),width=0.1)) +
  ylab('% inhibition relative to control') +
  xlab('') +
  ggtitle('b')+
  theme(legend.pos='none',
        plot.title = element_text(size=20),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45,hjust = 1))

#===========================================================================================================================
# Résultats LMIA Teladorsagia circumcincta
#===========================================================================================================================
tela = lmia3[which(lmia3$souche=='tela' & lmia3$frac!='fam'),]
tela$frac = factor(tela$frac)
tela$var = factor(tela$var)
tela$var = factor(tela$var)
tela$plq = factor(tela$plq)

#-- data struc
table(tela$plq) #-- 2 plates
table(tela$var,tela$frac)

#-- Reformat for analysis
tela$mol = paste(tela$var,tela$frac,sep='.')
tela$solution='bla'
tela$solution[tela$mol=='e063.fa'] ='E063.Alkaloids'
tela$solution[tela$mol=='ener.fa'] ='ENERGY.Alkaloids'
tela$solution[tela$mol=='e063.fna'] ='E063.Non-alkaloids'
tela$solution[tela$mol=='ener.fna'] ='ENERGY.Non-alkaloids'
tela$solution[tela$mol=='e063.macer'] ='E063-Tot. Extract'
tela$solution[tela$mol=='ener.macer'] ='ENERGY-Tot. Extract'
tela$solution[tela$mol=='tm.pos'] ='Levamisole'
tela$solution[tela$mol=='tm.neg'] ='Water'
#tela$solution = factor(tela$solution)
tela$solution = factor(tela$solution,levels=c('Water', 'Levamisole',
                                            'E063-Tot. Extract','E063.Alkaloids','E063.Non-alkaloids',
                                            'ENERGY-Tot. Extract','ENERGY.Alkaloids','ENERGY.Non-alkaloids'))

m2 = lme(1-pmig ~ solution ,
         random=~1|plq, data=tela)

summary(m2)

# Fixed effects: pmig ~ solution 
#                               Value Std.Error DF t-value p-value
# (Intercept)                   1.000     0.059 39   16.95  0.0000
# solutionLevamisole           0.737     0.055 39   13.40  0.0000
# solutionE063-Tot. Extract    0.980     0.055 39   17.82  0.0000
# solutionE063.Alkaloids       0.510     0.055 39    9.28  0.0000
# solutionE063.Non-alkaloids   0.217     0.055 39    3.94  0.0003
# solutionENERGY-Tot. Extract  0.902     0.055 39   16.40  0.0000
# solutionENERGY.Alkaloids     0.723     0.055 39   13.16  0.0000
# solutionENERGY.Non-alkaloids 0.060     0.055 39    1.09  0.2818

summary(glht(m2,mcp(solution="Tukey")))
#   Levamisole - Water == 0                           0.7367     0.0550   13.40   <0.001 ***
#   E063-Tot. Extract - Water == 0                    0.9800     0.0550   17.82   <0.001 ***
#   E063.Alkaloids - Water == 0                       0.5100     0.0550    9.28   <0.001 ***
#   E063.Non-alkaloids - Water == 0                   0.2167     0.0550    3.94   0.0021 ** 
#   ENERGY-Tot. Extract - Water == 0                  0.9017     0.0550   16.40   <0.001 ***
#   ENERGY.Alkaloids - Water == 0                     0.7233     0.0550   13.16   <0.001 ***
#   ENERGY.Non-alkaloids - Water == 0                 0.0600     0.0550    1.09   0.9588    
#   E063-Tot. Extract - Levamisole == 0               0.2433     0.0550    4.43   <0.001 ***
#   E063.Alkaloids - Levamisole == 0                 -0.2267     0.0550   -4.12   0.0010 ** 
#   E063.Non-alkaloids - Levamisole == 0             -0.5200     0.0550   -9.46   <0.001 ***
#   ENERGY-Tot. Extract - Levamisole == 0             0.1650     0.0550    3.00   0.0546 .  
#   ENERGY.Alkaloids - Levamisole == 0               -0.0133     0.0550   -0.24   1.0000    
#   ENERGY.Non-alkaloids - Levamisole == 0           -0.6767     0.0550  -12.31   <0.001 ***
#   E063.Alkaloids - E063-Tot. Extract == 0          -0.4700     0.0550   -8.55   <0.001 ***
#   E063.Non-alkaloids - E063-Tot. Extract == 0      -0.7633     0.0550  -13.88   <0.001 ***
#   ENERGY-Tot. Extract - E063-Tot. Extract == 0     -0.0783     0.0550   -1.42   0.8459    
#   ENERGY.Alkaloids - E063-Tot. Extract == 0        -0.2567     0.0550   -4.67   <0.001 ***
#   ENERGY.Non-alkaloids - E063-Tot. Extract == 0    -0.9200     0.0550  -16.73   <0.001 ***
#   E063.Non-alkaloids - E063.Alkaloids == 0         -0.2933     0.0550   -5.34   <0.001 ***
#   ENERGY-Tot. Extract - E063.Alkaloids == 0         0.3917     0.0550    7.12   <0.001 ***
#   ENERGY.Alkaloids - E063.Alkaloids == 0            0.2133     0.0550    3.88   0.0026 ** 
#   ENERGY.Non-alkaloids - E063.Alkaloids == 0       -0.4500     0.0550   -8.18   <0.001 ***
#   ENERGY-Tot. Extract - E063.Non-alkaloids == 0     0.6850     0.0550   12.46   <0.001 ***
#   ENERGY.Alkaloids - E063.Non-alkaloids == 0        0.5067     0.0550    9.22   <0.001 ***
#   ENERGY.Non-alkaloids - E063.Non-alkaloids == 0   -0.1567     0.0550   -2.85   0.0830 .  
#   ENERGY.Alkaloids - ENERGY-Tot. Extract == 0      -0.1783     0.0550   -3.24   0.0259 *  
#   ENERGY.Non-alkaloids - ENERGY-Tot. Extract == 0  -0.8417     0.0550  -15.31   <0.001 ***
#   ENERGY.Non-alkaloids - ENERGY.Alkaloids == 0     -0.6633     0.0550  -12.06   <0.001 ***

c=summary(lsmeans(m2,'solution'))[,c('lsmean','SE','lower.CL','upper.CL')]
c$solution=levels(tela$solution)
c$std = aggregate(1-pmig ~ solution, data=tela,FUN=sd)[,2]
c
# lsmean    SE lower.CL upper.CL solution                std
# 1.76e-16 0.059  -0.7497    0.750 Water                0.1349
# 7.37e-01 0.059  -0.0130    1.486 Levamisole           0.0668
# 9.80e-01 0.059   0.2303    1.730 E063-Tot. Extract    0.0268
# 5.10e-01 0.059  -0.2397    1.260 E063.Alkaloids       0.0914
# 2.17e-01 0.059  -0.5330    0.966 E063.Non-alkaloids   0.0819
# 9.02e-01 0.059   0.1520    1.651 ENERGY-Tot. Extract  0.0895
# 7.23e-01 0.059  -0.0263    1.473 ENERGY.Alkaloids     0.0979
# 6.00e-02 0.059  -0.6897    0.810 ENERGY.Non-alkaloids 0.1881

p1c = ggplot(c,aes(x=solution,y=lsmean,fill=solution))+
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('#a50026','#d73027','#f46d43',
                             '#66bd63','#1a9850','#006837',
                             '#542788','#92c5de'))+
  theme_classic() +
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(-0.31,1)) +
  geom_errorbar(aes(ymin=(c$lsmean-c$SE),ymax=(c$lsmean + c$SE),width=0.1)) +
  ylab('% inhibition relative to control') +
  xlab('') +
  ggtitle('c')+
  theme(legend.pos='none',
        plot.title = element_text(size=20),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45,hjust = 1))

###### FIGURE 2
pdf(file=paste0(draft,'Figure2.pdf'),width = 14, height=14)
multiplot(p1a,p1c,p1b,cols=2)
dev.off()
######

########--------------------------------------------
#= EC50 alcaloïdes E063 ENERGY WEY KOK
#= ggplot with drc http://stackoverflow.com/questions/38169765/drc-drc-plot-with-ggplot2
########--------------------------------------------
mes = read.csv(file="./LMIA_M_EA_FNA/EC50_ALC_COMBILEV_300615.csv",
             sep=";",
             header=TRUE)

head(mes)
#   souche  mol conc plq rep l3
# 1    wey ener 8.05   1   1 12
# 2    wey ener 8.05   1   2 14
# 3    wey ener 6.04   1   1 20
# 4    wey ener 6.04   1   2 22
# 5    wey ener 4.03   1   1 49
# 6    wey ener 4.03   1   2 42

#########----------------------	QC check 
aggregate(l3 ~ souche,data=mes[mes$mol=="neg",],FUN=mean)
# souche   l3
# 1    kok 54.2
# 2    wey 76.5

aggregate(l3 ~ plq + souche,data=mes[mes$mol=="neg",],FUN=mean)
#   plq souche   l3
# 1   1    kok 54.7
# 2   2    kok 50.0
# 3   3    kok 51.0
# 4   4    kok 61.0
# 5   1    wey 77.3
# 6   2    wey 75.0
# 7   3    wey 78.0
# 8   4    wey 75.7

#########---------------	Calcul des pourcentages de migration
mes$pmig=NA
mcont.1=mean(mes$l3[which(mes$souche=="wey" & mes$mol=="neg" & mes$plq==1)])
mcont.2=mean(mes$l3[which(mes$souche=="wey" & mes$mol=="neg" & mes$plq==2)])
mcont.3=mean(mes$l3[which(mes$souche=="wey" & mes$mol=="neg" & mes$plq==3)])
mcont.4=mean(mes$l3[which(mes$souche=="wey" & mes$mol=="neg" & mes$plq==4)])

mcont.5=mean(mes$l3[which(mes$souche=="kok" & mes$mol=="neg" & mes$plq==1)])
mcont.6=mean(mes$l3[which(mes$souche=="kok" & mes$mol=="neg" & mes$plq==2)])
mcont.7=mean(mes$l3[which(mes$souche=="kok" & mes$mol=="neg" & mes$plq==3)])
mcont.8=mean(mes$l3[which(mes$souche=="kok" & mes$mol=="neg" & mes$plq==4)])

mes$pmig[which(mes$souche=="wey" & mes$plq==1)]=mes$l3[which(mes$souche=="wey" & mes$plq==1)]/mcont.1
mes$pmig[which(mes$souche=="wey" & mes$plq==2)]=mes$l3[which(mes$souche=="wey" & mes$plq==2)]/mcont.2
mes$pmig[which(mes$souche=="wey" & mes$plq==3)]=mes$l3[which(mes$souche=="wey" & mes$plq==3)]/mcont.3
mes$pmig[which(mes$souche=="wey" & mes$plq==4)]=mes$l3[which(mes$souche=="wey" & mes$plq==4)]/mcont.4

mes$pmig[which(mes$souche=="kok" & mes$plq==1)]=mes$l3[which(mes$souche=="kok" & mes$plq==1)]/mcont.5
mes$pmig[which(mes$souche=="kok" & mes$plq==2)]=mes$l3[which(mes$souche=="kok" & mes$plq==2)]/mcont.6
mes$pmig[which(mes$souche=="kok" & mes$plq==3)]=mes$l3[which(mes$souche=="kok" & mes$plq==3)]/mcont.7
mes$pmig[which(mes$souche=="kok" & mes$plq==4)]=mes$l3[which(mes$souche=="kok" & mes$plq==4)]/mcont.8

#########-------------------------- LMIA ALK SUSCEPTIBLE ISOLATE

#---- Effet dose des alcaloides de ENERGY sur Weybridge
#- Jeux de donnees WEYBRIDGE ENER
ENERw=mes[mes$souche=="wey" & mes$mol=="ener"|
            mes$souche=="wey" & mes$mol=="neg",]
ENERw$mol="ENERGY"

#- Eliminer conc=2.01
ENERw1=ENERw[ENERw$conc!="2.01",]

#----- Effet dose des alcaloides de E063 sur Weybridge
#- Jeux de donnees
E063w=mes[mes$souche=="wey" & mes$mol=="e063"|mes$souche=="wey" & mes$mol=="neg",]
E063w$mol="E063"

#Eliminer conc=2.5
E063w1=E063w[E063w$conc!="2.5",]

#-----	Combinaison des donnees E063 et ENERGY pour WEYBRIDGE
COMBW=rbind(E063w1,ENERw1)
COMBW=COMBW[COMBW$plq==1 |COMBW$plq==3,]
COMBW$mol=factor(COMBW$mol)

#-- 1st model : different slopes, different ED betw plates 
combw = drm(pmig ~ conc, mol,data=COMBW,fct = LL.2())

# predictions and confidence intervals.
demo.fits <- expand.grid(conc=exp(seq(log(1.00e-04), log(10), length=100))) 
demo.fits$mol=rep(levels(COMBW$mol)[1])
demo.fits1 <- expand.grid(conc=exp(seq(log(1.00e-04), log(10), length=100))) 
demo.fits1$mol=rep(levels(COMBW$mol)[2])
demo.fits=rbind(demo.fits,demo.fits1)

# new data with predictions
pm <- predict(combw, newdata=demo.fits, interval="confidence") 
demo.fits$p <- pm[,1]
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]

COMBW$XX = COMBW$conc
COMBW$XX[COMBW$XX == 0] = 1.00e-09

p3a=ggplot(COMBW, aes(x = XX, y =(1-pmig),gp=mol)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~mol)+
  geom_ribbon(data=demo.fits, 
              aes(x=conc, y=1-p, ymin=(1-pmin), ymax=(1-pmax),
                  fill=mol), alpha=0.2) +
  geom_line(data=demo.fits, aes(x=conc, y=(1-p))) +
  ggtitle('c')+
  scale_fill_discrete(name = "Lupin variety") +
  ylab('% inhibition relative to control') +
  xlab('Concentration (mg/mL)') +
  theme(legend.position='bottom',
        title = element_text(size=18),
        axis.title=element_text(size=16),
        axis.text=element_text(size=14))

ED(combw,50,interva="delta")
#        Estimate Std. Error   Lower   Upper
#e063:50  9.29796    0.72644 7.81826 10.7777
#ener:50  4.53852    0.20871 4.11340  4.9636

compParm(combw,'e','-')
#             Estimate Std. Error t-value p-value
# E063-ENERGY    4.759      0.756   6.297       0

#########-------------------------- LMIA ALK RESISTANT ISOLATE

#----------------- Effet dose des alcaloides de ENERGY sur KOKSTAD
# Jeux de donnees
ENERk=mes[mes$souche=="kok" & mes$mol=="ener"|mes$souche=="kok" & mes$mol=="neg",]
ENERk$mol="ENERGY"

#eliminer conc=2.01 et plq 3 rep 1 conc 4.025
ENERk1=ENERk[ENERk$conc!="2.01",]

#----------------- Effet dose des alcaloedes de E063 sur KOKSTAD
# Jeux de donnees
E063k=mes[mes$souche=="kok" & mes$mol=="e063"|mes$souche=="kok" & mes$mol=="neg",]
E063k$mol="E063"

#eliminer conc=2.5
E063k1=E063k[E063k$conc!="2.5",]

##--------------	Combinaison des donnees E063 et ENERGY pour KOKSTAD
COMBK=rbind(E063k1,ENERk1)
COMBK=COMBK[COMBK$plq==1 |COMBK$plq==3,]
COMBK$mol=factor(COMBK$mol)

#-- Test models
combk2 = drm(pmig ~ conc,mol,data=COMBK,fct = LL.2())
combk3 = drm(pmig ~ conc,mol,data=COMBK,fct = LL.3())
combk4 = drm(pmig ~ conc,mol,data=COMBK,fct = LL.4())
anova(combk2,combk3)
anova(combk2,combk4)

#-- 1st model : different slopes, different ED betw plates 
combk = drm(pmig ~ conc,mol,data=COMBK,fct = LL.2())
summary(combk)

# predictions and confidence intervals.
demo.fits <- expand.grid(conc=exp(seq(log(1.00e-04), log(10), length=100))) 
demo.fits$mol=rep(levels(COMBK$mol)[1])
demo.fits1 <- expand.grid(conc=exp(seq(log(1.00e-04), log(10), length=100))) 
demo.fits1$mol=rep(levels(COMBK$mol)[2])
demo.fits=rbind(demo.fits,demo.fits1)

# new data with predictions
pm <- predict(combk, newdata=demo.fits, interval="confidence") 
demo.fits$p <- pm[,1]
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]

COMBK$XX = COMBK$conc
COMBK$XX[COMBK$XX == 0] = 1.00e-09

p3b=ggplot(COMBK, aes(x = XX, y = (1-pmig),gp=mol)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~mol)+
  geom_ribbon(data=demo.fits, 
              aes(x=conc, y=1-p, ymin=(1-pmin), ymax=(1-pmax),
                  fill=mol), alpha=0.2) +
  geom_line(data=demo.fits, aes(x=conc, y=(1-p))) +
  scale_fill_discrete(name = "Lupin variety") +
  ggtitle('d')+
  ylab('% inhibition relative to control') +
  xlab('Concentration (mg/mL)') +
  theme(legend.position='bottom',
        title = element_text(size=18),
        axis.title=element_text(size=16),
        axis.text=element_text(size=14))

ED(combk,50,interva="delta")
#             Estimate Std. Error Lower Upper
# e:E063:50       8.38       1.15  6.05 10.72
# e:ENERGY:50     5.59       0.47  4.64  6.55

compParm(combk,'e','-')
#             Estimate Std. Error t-value p-value
# E063-ENERGY     2.79       1.24    2.26    0.03

#### FIGURE 3
pdf(file = paste0(draft,'Figure3.pdf'),width=14,height=6)
multiplot(p3a,p3b,cols=2)
dev.off()
#### FIGURE 3

########--------------------------------------------
##--- Fraction assays (ALMA)
########--------------------------------------------
require(RColorBrewer)
require(dplyr)

alma = read.csv(file=paste0(draft,'../Fractions/alma_fractions.csv'),header=T,sep=';')

## convert to compound numbers
alma$Compound = 'unknown'
alma$Compound[alma$Mol=='D1'] = 'C1'
alma$Compound[alma$Mol=='D3'] = 'C2'
alma$Compound[alma$Mol=='A1'] = 'C3'
alma$Compound[alma$Mol=='L4'] = 'C4'
alma$Compound[alma$Mol=='D5'] = 'C5'
alma$Compound[alma$Mol=='J1'] = 'C6'
alma$Compound[alma$Mol=='Levamisole'] = 'LEV'
alma$Compound[alma$Mol=='Control'] = 'CTL'
alma$Compound[alma$Mol=='Alc'] = 'ALC'

## Convert to compound numbers
alma = alma[alma$Compound !='unknown',]

alma$Drug = 'Drug-susceptible'
alma$Drug[alma$Isolate=='Resistant'] = 'Drug-resistant'

## Work with 4h incubation time
alma4 = alma[alma$Time ==4 ,]
alma4$Mol = factor(alma4$Mol)
alma4$Compound = factor(alma4$Compound)
alma4$MI = factor(paste0(alma4$Isolate,'-',alma4$Compound))

a = aggregate(MigFraction ~ MI,data = alma4,FUN=mean)
a$std = aggregate(MigFraction ~ MI,data = alma4,FUN=sd)[,2]
a
#                 MI MigFraction       std
# 1    Resistant-ALC    78.66667 3.5725808
# 2     Resistant-C1    62.40000 2.8792360
# 3    Resistant-CTL   100.00000 2.3643181
# 4    Resistant-LEV    87.86667 7.9103308
# 5  Susceptible-ALC    68.46667 3.6115555
# 6   Susceptible-C1    53.63333 3.4239354
# 7   Susceptible-C2    56.30000 6.1538606
# 8   Susceptible-C3    61.93333 3.7527767
# 9   Susceptible-C4    76.00000 9.6814255
# 10  Susceptible-C5    92.46667 3.6170891
# 11  Susceptible-C6    94.80000 9.4015956
# 12 Susceptible-CTL   100.00000 5.8283788
# 13 Susceptible-LEV     9.80000 0.6244998

alma4 <- alma4 %>% group_by(MI) %>% mutate(med = median(MigFraction))
#lab = unique(as.character(alma4$Compound))[order(unique(as.numeric(alma4$Compound)))]
# brk = as.numeric(sort(unique(as.numeric(alma4$Compound))))

p_alma1 = ggplot(alma4[alma4$Isolate=='Susceptible',],aes(x = Compound,y = MigFraction,col = Compound,group = Compound)) +
  geom_point(size = 3) + facet_wrap(~ Drug) +
  theme_bw() + xlab('Compound') + ylab('Migration fraction (%)') +
  geom_segment(aes(x = as.numeric(Compound) - 0.2, xend = as.numeric(Compound) + 0.2, y = med, yend = med)) +
  scale_colour_manual(values = c("#081d58","#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8", "#FEE090","#FDAE61","#F46D43")) + # "#D73027" "#A50026"rev(brewer.pal(11,"RdYlBu"))) +
  scale_y_continuous(limits = c(0,110),breaks = seq(0,110,20)) +
  theme(text=element_text(size = 16),legend.position = 'none')

rm(df2)
df2 = alma4[alma4$Isolate=='Resistant',]
df2$Compound = factor(df2$Compound)
p_alma2 = ggplot(df2,aes(x = Compound,y = MigFraction,col = Compound,group = Compound)) +
  geom_point(size = 3) + facet_wrap(~ Drug) +
  theme_bw() + xlab('Compound') + ylab('Migration fraction (%)') +
  geom_segment(aes(x = as.numeric(Compound) - 0.2, xend = as.numeric(Compound) + 0.2, y = med, yend = med)) +
  scale_colour_manual(values = c("#081d58","#313695","#FDAE61","#F46D43")) + #rev(brewer.pal(11,"RdYlBu"))) +
  scale_y_continuous(limits = c(0,110),breaks = seq(0,110,20)) +
  theme(text=element_text(size = 16),legend.position = 'none')

library(gridExtra)
gA <- ggplotGrob(p_alma1)
gB <- ggplotGrob(p_alma2)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)

##########------- FIGURE 4 --------------#########

pdf(file=paste0(draft,'Figure4alma.pdf'),width = 14,height = 11)
plot(p4 <- arrangeGrob(
  gA, gB, ncol = 2, widths = c(1.20, 0.40)))
dev.off()

###--- Potency - Drug Susceptible
alma4s = data.frame(alma4[alma4$Isolate=='Susceptible',])
alma4s$Mol <- relevel(alma4s$Mol,ref="Control")
alma4s$Compound <- relevel(alma4s$Compound,ref="CTL")

arr = array(0,(length(levels(alma4s$Mol))-1))
arr2 = array(0,(length(levels(alma4s$Mol))-1))

for(i in c(match(levels(alma4s$Mol),levels(alma4s$Mol)))[-1]){
  df = alma4s[alma4s$Mol %in% c("Control",levels(alma4s$Mol)[i]),]
  df$Mol = factor(df$Mol)
  k = wilcox.test(x=df$MigFraction[df$Mol==levels(alma4s$Mol)[i]],y=df$MigFraction[df$Mol=='Control'],
                  paired=FALSE,alternative=c("less"))
  arr[i-1] = k$p.value
  arr2[i-1] = k$statistic
  rm(df)
}

res = data.frame(cbind(arr,arr2))
colnames(res) = c('p.value','statistic')
res$Mol = c(levels(alma4s$Compound)[-1])
res

#   p.value statistic Mol
# 1    0.05         0 ALC
# 2    0.05         0  C1
# 3    0.05         0  C2
# 4    0.05         0  C3
# 5    0.20         2  C4
# 6    0.50         4  C5
# 7    0.05         0  C6
# 8    0.05         0 LEV

###--- Potency - Drug Resistant
alma4r = data.frame(alma4[alma4$Isolate=='Resistant',])
alma4r$Mol <- relevel(alma4r$Mol,ref="Control")
df = alma4r[alma4r$Mol %in% c("Control","D1"),]
df$Mol = factor(df$Mol)
wilcox.test(x=df$MigFraction[df$Mol=="D1"],y=df$MigFraction[df$Mol=='Control'],
                paired=FALSE,alternative=c("less"))
# data:  df$MigFraction[df$Mol == "D1"] and df$MigFraction[df$Mol == "Control"]
# W = 0, p-value = 0.05
# alternative hypothesis: true location shift is less than 0

df = alma4r[alma4r$Mol %in% c("Control","Alc"),]
wilcox.test(x=df$MigFraction[df$Mol=="Alc"],y=df$MigFraction[df$Mol=='Control'],
            paired=FALSE,alternative=c("less"))
# data:  df$MigFraction[df$Mol == "Alc"] and df$MigFraction[df$Mol == "Control"]
# W = 0, p-value = 0.05
# alternative hypothesis: true location shift is less than 0

########--------------------------------------------
#- EC50 LUPANINNE /// SUSCEPTIBLE ISOLATE ONLY 
########--------------------------------------------
setwd("./ODubois/res/")

###############-- lupanine
lup = read.csv(file='./ALCALOIDES/test de migration-gamme lévamisole (KOK WEY)- gamme lupanine (WEY) 030715.csv',header=T,sep=';')
lup2 = lup[lup$souche=='wey' & (lup$plq==2 | lup$plq==4),]

#- QC check: keep plate 4 only
aggregate(l3 ~ plq,data=lup2[lup2$mol=='tm',],FUN=mean)
#   plq l3
# 1   2 37
# 2   4 57

lup2 = lup2[lup2$plq==4,]
lup2$souche=NULL
lup2$pmig=lup2$l3/(mean(lup2$l3[lup2$mol=='tm']))
lup2$mol[lup2$mol=='tm']='lup'

##############-- Add levamisole reference 
lev=lup[lup$plq==3 & ((lup$mol=='lev' & lup$souche=='wey') | (lup$mol=='tm' & lup$souche=='wey')) ,]
lev$souche=NULL
lev$pmig=lev$l3/(mean(lev$l3[lev$mol=='tm']))
lev$mol[lev$mol=='tm']='lev'

alc = rbind(lup2,lev)
alc$mol=factor(alc$mol)

##############-- DRC model
#-- Test models
alcLMIA = drm(pmig ~ conc,mol,data=alc,fct = LL.2())
modelFit(alcLMIA)

#--  Model with 2 parameters 
summary(alcLMIA)
# Parameter estimates:
#       Estimate Std. Error  t-value p-value
# b:lup 1.22e-01   2.44e-02 5.01e+00    0.00
# b:lev 1.01e+00   1.94e-01 5.22e+00    0.00
# e:lup 1.50e+07   2.03e+07 7.41e-01    0.46
# e:lev 1.61e+00   2.71e-01 5.95e+00    0.00

ED(alcLMIA,50,interva="delta")
#           Estimate Std. Error     Lower     Upper
# e:lev:50  1.61e+00   2.71e-01  1.07e+00  2.15e+00
# e:lup:50  1.50e+07   2.03e+07 -2.55e+07  5.55e+07

compParm(alcLMIA,'e','-')
#           Estimate Std. Error   t-value p-value
# cyt-lev  1.53e+03   1.15e+03  1.33e+00    0.19
# lup-lev  1.50e+07   2.03e+07  7.41e-01    0.46

# predictions and confidence intervals.
# demo.fits <- expand.grid(conc=exp(seq(log(1.00e-04), log(3500), length=100))) 
# demo.fits$mol=rep(levels(alc$mol)[1])
demo.fits1 <- expand.grid(conc=exp(seq(log(0.000001), log(3500), length=100))) 
demo.fits1$mol=rep(levels(alc$mol)[2])
demo.fits2 <- expand.grid(conc=exp(seq(log(0.000001), log(3500), length=100)))
demo.fits2$mol=rep(levels(alc$mol)[3])

demo.fits=rbind(demo.fits1,demo.fits2)

# new data with predictions
pm <- predict(alcLMIA, newdata=demo.fits, interval="confidence") 
demo.fits$p <- pm[,1]
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]
demo.fits$mole='Lupanine'
demo.fits$mole[demo.fits$mol=='lev']='Levamisole'
demo.fits$mole=factor(demo.fits$mole)

alc$XX = alc$conc
alc2$mol=factor(alc2$mol)
alc2$mole='Lupanine'
alc2$mole[alc2$mol=='lev']='Levamisole'
alc2$mole=factor(alc2$mole)

####### Supplementary Figure 3

library(scales)
alc2$conc[alc2$conc==0]=0.000001
alc2$facet <- ifelse(alc2$conc == min(alc2$conc), 1, 2)
demo.fits$facet <- ifelse(demo.fits$conc == min(demo.fits$conc), 1, 2)

# Discard low conc values from demo.fits for aesthetics
demo.fits = subset(demo.fits, conc>.001)
pdf(file=paste0(draft,'supplementary_figure3.pdf'))
ggplot(alc2, aes(x = conc, y = (1-pmig),gp=mol)) +
  geom_point() + 
  theme_bw() +
  geom_ribbon(data = demo.fits,
              aes(x=conc, y=1-p, ymin=(1-pmin), ymax=(1-pmax),
                  fill=mole), alpha=0.2) +
  geom_line(data=demo.fits, aes(x=conc, y=(1-p))) +
  annotation_logticks(scaled = TRUE,sides="b") +
  scale_x_log10(breaks = c(0.000001, 10^(-2:4)), 
                labels = c(0, math_format()(-2:4))) +
  scale_fill_discrete(name = "Molecule") +
  facet_grid(~facet, scales = 'free', space = 'free') +
  ylab('% inhibition relative to control') +
  xlab('Concentration (µM)') +
  theme(legend.position='bottom',
        axis.title=element_text(size=16),
        axis.text=element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_blank())
dev.off()



