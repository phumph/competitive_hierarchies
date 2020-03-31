#### HEADER ####
### loading comp2_final and performing correlation analyses and preparing figures!
### PTH 10-APRIL-2015

# load comp2_final:

comp2 <- read.table("comp2_final.txt",T,"\t",row.names=1)

# first order of business.. perform PCA on all life history traits and look at how the components load up:
# look at ggplot version of biplot.. grab code from old r file:

library(lme4)
library(devtools)
library(ggbiplot)
library(ggplot2)
library(gplots)
library(gridExtra)


#### FIGURE 1A: PHYLOGENY WITH TRAITS AS HEATMAP COLS ####
# trying first with ggplot2 after collecting which trait data to display:

comp3 <- with(comp2,cbind(mu.mm,lambda.mm,A.mm,c.i,c.h,c.w,c.r,c.t))
# put all P. syr rows to NA for c.t:
comp3[1:26,8] <- NA
heatmap.2(comp3,dendrogram=c("none"),trace=c("none"),
          density.info=c("none"),Rowv=FALSE,Colv=FALSE,
          na.color="gray",col=colorRampPalette(c("midnightblue","white","darkorange2")),
          notecex=0.6,notecol="black",symbreaks=TRUE,cexRow=0.6,cexCol=0.6,scale="column")


#### FIGURE 1B: density plots (stacked) for both clades: ####
# use ggplot2..
psd1 <- ggplot(comp2, aes(mu.mm)) + geom_density(aes(linetype=clade)) +
            theme_bw() + 
            scale_x_continuous(limits=c(0,0.0015)) + #scale_y_continuous(limits=c(0,3000)) + 
            theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())

psd2 <- ggplot(comp2, aes(lambda.mm)) + geom_density(aes(linetype=clade)) +
  theme_bw() + 
  scale_x_continuous(limits=c(0,3000)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

psd3 <- ggplot(comp2, aes(A.mm)) + geom_density(aes(linetype=clade)) +
  theme_bw() + 
  scale_x_continuous(limits=c(0,1.0)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

psd4 <- ggplot(comp2, aes(c.i)) + geom_density(aes(linetype=clade)) +
  theme_bw() + 
  scale_x_continuous(limits=c(0,1)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

psd5 <- ggplot(comp2, aes(c.h)) + geom_density(aes(linetype=clade)) +
  theme_bw() + 
  scale_x_continuous(limits=c(0,1)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

psd6 <- ggplot(comp2, aes(c.w)) + geom_density(aes(linetype=clade)) +
  theme_bw() + 
  scale_x_continuous(limits=c(-1,1)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

psd7 <- ggplot(comp2, aes(c.r)) + geom_density(aes(linetype=clade)) +
  theme_bw() + 
  scale_x_continuous(limits=c(0,1)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

psd8 <- ggplot(subset(comp2,clade=="Pflu"), aes(c.t)) + geom_density() +
  theme_bw() + 
  scale_x_continuous(limits=c(0,1)) + #scale_y_continuous(limits=c(0,3000)) + 
  theme(legend.position = "none",plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

grid.arrange(psd1,psd2,psd3,psd4,psd5,psd6,psd7,psd8,nrow=8)

# now spitting out summary statistics for Fig 1B right margins..
comp3b <- as.data.frame(cbind(comp3,comp2[,9]))
dimnames(comp3b)[[2]][9] <- c("clade")
dimnames(comp3b)[[1]] <- rownames(comp2)
comp3b$clade <- factor(comp3b$clade)

# now export table of means as well as 95% CI

sem <- function(x,ci = FALSE)
{ 
  sem_out <- sd(x,na.rm=T)/sqrt(length(na.omit(x))) 
  
  if (ci == TRUE) {
    return(1.96*sem_out) 
  }
  else {
      return(sem_out)
  }
}
  
semd <- function(x,y,ci = FALSE)
{ 
  semd_out <- sqrt( ( (var(x, na.rm=T)/length(na.omit(x))) + (var(y, na.rm=T)/length(na.omit(y))) ) )
  
  if (ci == TRUE) {
    return(1.96*semd_out) 
  }
  else {
    return(semd_out)
  }
}

# capture all summary statistics
tmp1 <- data.frame(matrix(ncol=6,nrow=8))
names(tmp1) <- c("mean.Pf","mean.Ps","ci.Pf","ci.Ps","m.d","sem.d")
for (i in 1:(length(comp3[1,])))
{
  tmp1[i,1:2] <- tapply(comp3b[,i],comp3b$clade,mean,na.rm=T) # for trait i, calculate means of both clades
  tmp1[i,3:4] <- tapply(comp3b[,i],comp3b$clade,sem)*1.96     # for trait i, calculate sem of both clades
  tmp1[i,5] <- tmp1[i,1] - tmp1[i,2]
  tmp1[i,6] <- 1.96*semd(comp3b[,i][comp3b$clade=="1"],comp3b[,i][comp3b$clade=="2"]) # calculate sem of the difference between the means
}

write.table(tmp1,file="trait_summary_stats.txt",row.names=F,quote=F,sep="\t")


#### FIGURE 2A: CORRELATION MATRIXES FOR WITHIN-GROUP TRAIT CORRELATIONS ####
# first for P. syringae: minus row 20, minus column 8:
comp4 <- comp3[-c(20,27:40),] # for P. syr
comp5 <- comp3[-c(1:26),] # for P. flu

library(psych)

cPsyr <- corr.test(comp4,method="pearson",adjust="fdr")
cPflu <- corr.test(comp5,method="pearson",adjust="fdr")

# write tables for coloring squares by FDR-corrected p-value..
write.table(cPsyr$p,file="cPsyr_p.txt",sep="\t")
write.table(cPflu$p,file="cPflu_p.txt",sep="\t")

# now make heatmap similarly to the one in MolEcol 2014..
palette.breaks <- seq(-1,1,0.125)
#color.palette  <- colorRampPalette(c("midnightblue","white","darkorange2"))(length(palette.breaks) - 1)
color.palette  <- colorRampPalette(c("midnightblue","white","darkorange2"))(length(palette.breaks) - 1)
h1 <- heatmap.2(cPsyr$r,dendrogram = "none",scale = "none", trace = "none", key = TRUE, keysize = 1.5, density.info=c("none"), labRow = NA,labCol = NA,
                Rowv       = FALSE,
                Colv       = FALSE,  
                col    = color.palette,
                #col     = cm.colors,
                breaks = palette.breaks,
                symbreaks = TRUE,
                cellnote = round(cPsyr$r,2),
                notecex  = 0.9,
                notecol  = "black",
                colsep,
                rowsep,
                sepcolor="white",
                sepwidth=c(0.1,0.1)
)

h2 <- heatmap.2(cPflu$r,dendrogram = "none",scale = "none", trace = "none", key = TRUE, keysize = 1.5, density.info=c("none"), labRow = NA,labCol = NA,
                Rowv       = FALSE,
                Colv       = FALSE,  
                col    = color.palette,
                #col     = cm.colors,
                breaks = palette.breaks,
                symbreaks = TRUE,
                cellnote = round(cPflu$r,2),
                notecex  = 0.9,
                notecol  = "black",
                colsep,
                rowsep,
                sepcolor="white",
                sepwidth=c(0.1,0.1)
)


#### FIGURE 2B: PCA PLOTS ####
# perform pca on sub-matrix of comp2, including the following:
# mu.mm
# lambda.mm
# A.mm
# int.mm
# c.i
# c.h
# c.w
# c.r
# c.t

pca_traits1 <- with(comp3[-20,], cbind(mu.mm,lambda.mm,A.mm,c.w))
pca_traits1s <- with(comp3[-c(20,27:40),], cbind(mu.mm,lambda.mm,A.mm,c.w,c.r))
pca_traits1f <- with(comp3[-c(1:25),], cbind(mu.mm,lambda.mm,A.mm,c.w,c.t))
pca1 <- prcomp(pca_traits1,center=TRUE,scale.=TRUE)
pca1s <- prcomp(pca_traits1s,center=TRUE,scale.=TRUE)
pca1f <- prcomp(pca_traits1f,center=TRUE,scale.=TRUE)
plot(pca1)
summary(pca1)

clades <- comp3[-20,]$clade
strain.id <- row.names(comp3[-20,])
strain.id <- row.names(comp3[-c(20,27:40),])
strain.id <- row.names(comp3[-20,])
# primary plot without dot labels
g12 <- ggbiplot(pca1s, obs.scale = 1, var.scale = 1, choices=1:2,labels=strain.id,
                ellipse = TRUE, ellipse.prob=0.95,
                circle = FALSE)

g12 <- ggbiplot(pca1f, obs.scale = 1, var.scale = 1, choices=1:2,
                ellipse = TRUE, ellipse.prob=0.95,
                circle = FALSE)


g12 <- ggbiplot(pca1, obs.scale = 1, var.scale = 1, choices=1:2,
                groups = clades, ellipse = TRUE, ellipse.prob=0.95,
                circle = FALSE)
g12 <- g12 + theme_bw() + theme(legend.position = c(1, 1),legend.justification = c(1,1))
g12 <- g12 + scale_color_manual(values = c("gray", "black"))

# spit out identical plot with plot labels instead; adjust the position to have them on same 
g12b <- ggbiplot(pca1, obs.scale = 1, var.scale = 1, choices=1:2,
                groups = clades, ellipse = TRUE, ellipse.prob=0.95, labels=strain.id,
                circle = FALSE)
g12b <- g12b + theme_bw() + theme(legend.position = c(1, 1),legend.justification = c(1,1))
g12b <- g12b + scale_color_manual(values = c("gray", "black"))


# exporting additional PCA with only core canonical variables present in both clades (i.e. exclude c.t, c.r)
pca_traits2 <- with(comp2[-20,], cbind(mu.mm,lambda.mm,A.mm,c.w))
pca2 <- prcomp(pca_traits2,center=TRUE,scale.=TRUE)
g12c <- ggbiplot(pca2, obs.scale = 1, var.scale = 1, choices=1:2,
                groups = clades, ellipse = TRUE, ellipse.prob=0.95,
                circle = FALSE)
g12c <- g12c + theme_bw() + theme(legend.position = c(1, 1),legend.justification = c(1,1))
g12c <- g12c + scale_color_manual(values = c("gray", "black"))

# spit out identical plot with plot labels instead; adjust the position to have them on same 
g12d <- ggbiplot(pca2, obs.scale = 1, var.scale = 1, choices=1:2,
                 groups = clades, ellipse = TRUE, ellipse.prob=0.95, labels=strain.id,
                 circle = FALSE)
g12d <- g12d + theme_bw() + theme(legend.position = c(1, 1),legend.justification = c(1,1))
g12d <- g12d + scale_color_manual(values = c("gray", "black"))

grid.arrange(g12,g12b,g12c,g12d,ncol=2,nrow=2)

# ultimately used the first PCA since c.r is important, and c.t. is quite negatively correlated with c.r.

#### FIGURE 3: LOGISTIC REGRESSIONG FROM GLMS OF INTERACTION TYPE BY DISTANCE ####

# performing GLMs and logistic regressions with data from interaction frequencies
# first: examining models to produce GLM table:

## examining correlations between all pairs of distance matrixes.. phylo vs. traits.

## loading all dists:
d2 <- read.table("../dists/dist_mm_new.txt",T,"\t",row.names=1)
muD <- read.table("../dists/mu_mm_new.txt",T,"\t",row.names=1)
lamD <- read.table("../dists/lambda_mm_new.txt",T,"\t",row.names=1)
aD <- read.table("../dists/A_mm_new.txt",T,"\t",row.names=1)
cwD <- read.table("../dists/cw_dist_new.txt",T,"\t",row.names=1)
coD <- read.table("../dists/co_dist_new.txt",T,"\t",row.names=1)
cdD <- read.table("../dists/cd_dist_new.txt",T,"\t",row.names=1)
crD <- read.table("../dists/cr_dist_new.txt",T,"\t",row.names=1)
ctD <- read.table("../dists/ct_dist_new.txt",T,"\t",row.names=1)

library(vegan)
mantel(as.dist(d2),as.dist(muD),na.rm=TRUE) # **p=0.001
mantel(as.dist(d2),as.dist(lamD),na.rm=TRUE) # **p=0.001
mantel(as.dist(d2),as.dist(aD),na.rm=TRUE) # *p=0.018

# now using competitive traits:
mantel(as.dist(d2),as.dist(cwD),na.rm=TRUE) # ***p<0.001
mantel(as.dist(d2),as.dist(coD),na.rm=TRUE) # ***p<0.001
mantel(as.dist(d2),as.dist(cdD),na.rm=TRUE) # ***p<0.001
mantel(as.dist(d2),as.dist(crD),na.rm=TRUE) # p=0.3


# trying only within P. syringae!
mantel(as.dist(d2[1:26,1:26]),as.dist(muD[1:26,1:26]),na.rm=TRUE) 
mantel(as.dist(d2[1:26,1:26]),as.dist(lamD[1:26,1:26]),na.rm=TRUE) 
mantel(as.dist(d2[1:26,1:26]),as.dist(aD[1:26,1:26]),na.rm=TRUE) 

mantel(as.dist(d2[1:26,1:26]),as.dist(cwD[1:26,1:26]),na.rm=TRUE) 
mantel(as.dist(d2[1:26,1:26]),as.dist(coD[1:26,1:26]),na.rm=TRUE) 
mantel(as.dist(d2[1:26,1:26]),as.dist(cdD[1:26,1:26]),na.rm=TRUE) 

# trying only within P. fluorescens!
mantel(as.dist(d2[27:40,27:40]),as.dist(muD[27:40,27:40]),na.rm=TRUE)
mantel(as.dist(d2[27:40,27:40]),as.dist(lamD[27:40,27:40]),na.rm=TRUE)
mantel(as.dist(d2[27:40,27:40]),as.dist(aD[27:40,27:40]),na.rm=TRUE)

mantel(as.dist(d2[27:40,27:40]),as.dist(cwD[27:40,27:40]),na.rm=TRUE)
mantel(as.dist(d2[27:40,27:40]),as.dist(coD[27:40,27:40]),na.rm=TRUE)
mantel(as.dist(d2[27:40,27:40]),as.dist(cdD[27:40,27:40]),na.rm=TRUE)

mantel(as.dist(d2[27:40,27:40]),as.dist(ctD[27:40,27:40]),na.rm=TRUE) # ns


mantel(as.dist(muD),as.dist(lamD),na.rm=TRUE) # ns
mantel(as.dist(muD),as.dist(aD),na.rm=TRUE) # **p=0.007
mantel(as.dist(lamD),as.dist(aD),na.rm=TRUE) # ns

int_glm <- read.table("ints_for_glms.txt","\t",header=T) # OK all set for GLMs:
glm_res <- list()

glm_res[[1]] <- summary(glm1 <- glm(RI~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="W" & j.clade=="Pflu")),family="binomial"))$coef
glm_res[[2]] <- summary(glm1 <- glm(RI~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="W" & j.clade=="Psyr")),family="binomial"))$coef
glm_res[[3]] <- summary(glm1 <- glm(RI~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="B")),family="binomial"))$coef
glm_res[[4]] <- summary(glm1 <- glm(RI~pdist+tdist1+tdist2+tdist3,data=int_glm,family="binomial"))$coef

glm_res[[5]] <- summary(glm1 <- glm(RNI~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="W" & j.clade=="Pflu")),family="binomial"))$coef
glm_res[[6]] <- summary(glm1 <- glm(RNI~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="W" & j.clade=="Psyr")),family="binomial"))$coef
glm_res[[7]] <- summary(glm1 <- glm(RNI~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="B")),family="binomial"))$coef
glm_res[[8]] <- summary(glm1 <- glm(RNI~pdist+tdist1+tdist2+tdist3,data=int_glm,family="binomial"))$coef

glm_res[[9]] <- summary(glm1 <- glm(ASYM~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="W" & j.clade=="Pflu")),family="binomial"))$coef
glm_res[[10]] <- summary(glm1 <- glm(ASYM~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="W" & j.clade=="Psyr")),family="binomial"))$coef
glm_res[[11]] <- summary(glm1 <- glm(ASYM~pdist+tdist1+tdist2+tdist3,data=subset(int_glm,(WB=="B")),family="binomial"))$coef
glm_res[[12]] <- summary(glm1 <- glm(ASYM~pdist+tdist1+tdist2+tdist3,data=int_glm,family="binomial"))$coef

write.table(glm_res,file="glm_res_full.txt",row.names=T,sep="\t")

## FIGURE 3A ## small barplots for proportion of interactions of each time within each bin:
# using only non-self interactions
  # Within - Psyr; Pflu
  # Between
  # Total

wPs <- with(subset(int_glm,(int_glm$j.clade=="Psyr" & int_glm$WB=="W")), cbind(sum(RI),sum(RNI),sum(ASYM)))
wPf <- with(subset(int_glm,(int_glm$j.clade=="Pflu" & int_glm$WB=="W")), cbind(sum(RI),sum(RNI),sum(ASYM)))
Btw <- with(subset(int_glm,(int_glm$WB=="B")), cbind(sum(RI),sum(RNI),sum(ASYM)))
iAll <- with(int_glm, cbind(sum(int_glm$RI),sum(RNI),sum(ASYM)))

int_tab <- t(rbind(wPs,wPf,Btw,iAll))
rownames(int_tab) <- c("RI","RNI","ASYM")
colnames(int_tab) <- c("wPsur","wPflu","Btwn","Tots")
int_prop <- prop.table(int_tab,margin=2)

# Plotting code for barplot:
par(mar=c(5, 4, 4, 2) + 0.1)
barplot(int_prop, col = c("black","gray60","white"), width = 1,las=2,cex.axis=0.75)
legend("topright", inset = c(-0.25, 0), fill = c("black","gray60","white"),legend = rownames(int_tab),bty="n",cex=0.75)

## FIGURE 3B ## logistic regressions on phylogenetic distance within, between, and total:
# import int_glm with self interactions in there; export frequencies also as barplots.

# load int_glm_s including self interactions for p-dist calculation
int_glm_s <- read.table("ints_for_glms_w_self.txt","\t",header=T) # OK all set for GLMs:


gx1 <- ggplot(int_glm_s, aes(x = pdist, y = ASYM)) +
  geom_jitter(position = position_jitter(height = .025), alpha=0.5) + scale_y_continuous(limits = c(-0.025, 1.025)) +
  #geom_point() + 
  stat_smooth( aes(y = ASYM),  method="glm", family="binomial", se=T) +
  theme_bw()

gx2 <- ggplot(int_glm_s, aes(x = pdist, y = RI)) +
  geom_jitter(position = position_jitter(height = .025), alpha=0.5) + scale_y_continuous(limits = c(-0.025, 1.025)) +
  #geom_point() + 
  stat_smooth( aes(y = RI),  method="glm", family="binomial", se=T) +
  theme_bw()

gx3 <- ggplot(int_glm_s, aes(x = pdist, y = RNI)) +
  geom_jitter(position = position_jitter(height = .025), alpha=0.5) + scale_y_continuous(limits = c(-0.025, 1.025)) +
  #geom_point() + 
  stat_smooth( aes(y = RNI),  method="glm", family="binomial", se=T) +
  theme_bw()

grid.arrange(gx1,gx2,gx3,nrow=3)

# now spit out GLM results from these plots, both with and without self:
glm_res_self <- list()

glm_res_self[[1]] <- summary(glm1 <- glm(RI~pdist,data=int_glm_s,family="binomial"))$coef
glm_res_self[[2]] <- summary(glm1 <- glm(RNI~pdist,data=int_glm_s,family="binomial"))$coef
glm_res_self[[3]] <- summary(glm1 <- glm(ASYM~pdist,data=int_glm_s,family="binomial"))$coef

glm_res_self[[4]] <- summary(glm1 <- glm(RI~pdist,data=int_glm,family="binomial"))$coef
glm_res_self[[5]] <- summary(glm1 <- glm(RNI~pdist,data=int_glm,family="binomial"))$coef
glm_res_self[[6]] <- summary(glm1 <- glm(ASYM~pdist,data=int_glm,family="binomial"))$coef

write.table(glm_res_self,file="glm_res_self.txt",sep="\t",quote=F)

#
### SUPPLEMENTAL FIGURE 1: GROWTH CURVES ####
# adopt code from previous try at plotting:
# load growth rate data for each clade separately
g1 <- read.table("MM_gcurvedata_Pflu.txt",T,"\t")
g2 <- read.table("MM_gcurvedata_Psyr.txt",T,"\t")

par(mfrow=c(2,3),mai=c(0.3,0.4,0.3,0.3))
plot(g1[,1],g1[,15],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.6),las=2,cex.axis=0.6)
points(g1[,1],g1[,14],type="l",lty=2)
points(g1[,1],g1[,13],type="l",lty=3)
legend("topleft",c(dimnames(g1)[[2]][15],dimnames(g1)[[2]][14],dimnames(g1)[[2]][13]),cex=0.66,lty=c(1,2,3),bty="n")

# row 2
plot(g1[,1],g1[,12],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.6),las=2,cex.axis=0.66)
points(g1[,1],g1[,11],type="l",lty=2)
legend("topleft",c(dimnames(g1)[[2]][12],dimnames(g1)[[2]][11]),cex=0.66,lty=c(1,2),bty="n")

# row 3
plot(g1[,1],g1[,10],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.6),las=2,cex.axis=0.66)
points(g1[,1],g1[,9],type="l",lty=2)
legend("topleft",c(dimnames(g1)[[2]][10],dimnames(g1)[[2]][9]),cex=0.66,lty=c(1,2),bty="n")

# row 4
plot(g1[,1],g1[,7],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.6),las=2,cex.axis=0.66)
points(g1[,1],g1[,6],type="l",lty=2)
points(g1[,1],g1[,5],type="l",lty=3)
points(g1[,1],g1[,4],type="l",lty=4)
legend("topleft",c(dimnames(g1)[[2]][7],dimnames(g1)[[2]][6],dimnames(g1)[[2]][5],dimnames(g1)[[2]][4]),cex=0.66,lty=c(1,2,3,4),bty="n")

# row 5
plot(g1[,1],g1[,3],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.6),las=2,cex.axis=0.66)
points(g1[,1],g1[,2],type="l",lty=2)
legend("topleft",c(dimnames(g1)[[2]][3],dimnames(g1)[[2]][2]),cex=0.66,lty=c(1,2),bty="n")

# row 6
plot(g1[,1],g1[,8],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.6),las=2,cex.axis=0.66)
legend("topleft",c(dimnames(g1)[[2]][8]),cex=0.66,lty=c(1),bty="n")

# now do P. syringae plots.. save separately
par(mfrow=c(3,3),mai=c(0.3,0.4,0.3,0.3)) # 9 plots

# first plot: 10A, 21B, 20A, 17A
plot(g2[,1],g2[,9],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,11],type="l",lty=2)
points(g2[,1],g2[,12],type="l",lty=3)
points(g2[,1],g2[,10],type="l",lty=4)
legend("topleft",c(dimnames(g2)[[2]][9],dimnames(g2)[[2]][11],dimnames(g2)[[2]][12],dimnames(g2)[[2]][10]),cex=0.66,lty=c(1,2,3,4),bty="n")

# second plot: 24A, 39C, 14B, 16A
plot(g2[,1],g2[,13],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,14],type="l",lty=2)
points(g2[,1],g2[,15],type="l",lty=3)
points(g2[,1],g2[,16],type="l",lty=4)
legend("topleft",c(dimnames(g2)[[2]][13],dimnames(g2)[[2]][14],dimnames(g2)[[2]][15],dimnames(g2)[[2]][16]),cex=0.66,lty=c(1,2,3,4),bty="n")

# thid plot: 46B 19A
plot(g2[,1],g2[,25],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,26],type="l",lty=2)
legend("topleft",c(dimnames(g2)[[2]][25],dimnames(g2)[[2]][26]),cex=0.66,lty=c(1,2),bty="n")

# fourth plot: 06C, 02A, 08A, Psm4326
plot(g2[,1],g2[,24],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,22],type="l",lty=2)
points(g2[,1],g2[,23],type="l",lty=3)
points(g2[,1],g2[,21],type="l",lty=4)
legend("topleft",c(dimnames(g2)[[2]][24],dimnames(g2)[[2]][22],dimnames(g2)[[2]][23],dimnames(g2)[[2]][21]),cex=0.66,lty=c(1,2,3,4),bty="n")

# fifth plot: 22D, 07A, 22B
plot(g2[,1],g2[,19],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,20],type="l",lty=2)
points(g2[,1],g2[,18],type="l",lty=3)
legend("topleft",c(dimnames(g2)[[2]][19],dimnames(g2)[[2]][20],dimnames(g2)[[2]][18]),cex=0.66,lty=c(1,2,3),bty="n")

# sixth plot: 26B, 39F
plot(g2[,1],g2[,17],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,2],type="l",lty=2)
legend("topleft",c(dimnames(g2)[[2]][17],dimnames(g2)[[2]][2]),cex=0.66,lty=c(1,2),bty="n")

# seventh plot: 08C, 23A
plot(g2[,1],g2[,7],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,8],type="l",lty=2)
legend("topleft",c(dimnames(g2)[[2]][7],dimnames(g2)[[2]][8]),cex=0.66,lty=c(1,2),bty="n")

# eigth plot: 21A, 20B
plot(g2[,1],g2[,4],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
points(g2[,1],g2[,5],type="l",lty=2)
legend("topleft",c(dimnames(g2)[[2]][4],dimnames(g2)[[2]][5]),cex=0.66,lty=c(1,2),bty="n")

# ninth and last plot: 08B
plot(g2[,1],g2[,6],type="l",lty=1,xlim=c(0,3600),ylim=c(0,0.8),las=2,cex.axis=0.6)
legend("topleft",c(dimnames(g2)[[2]][6]),cex=0.66,lty=c(1),bty="n")



#### SUPPLEMENTAL FIGURE 2: MM V KB COMPARISONS/TRADE-OFFS AND WITHIN-KB CORRELATIONS ####

# first part: special line graphs; second part: correlation heatmap, as in Fig 2A (main text):
par(mfrow=c(1,3),mai=c(0.4,0.4,0.26,0.26))

# first for r_MM v r_KB
plot(comp2[,1],comp2[,5],pch="",ylim=c(0,0.005),xlim=c(0,0.0012), cex.axis=0.75,ylab="",xlab="",las=2,bty="n")
rug(comp2[comp2$clade=="Psyr",1],side=1,col="darkslateblue",lwd=1)
rug(comp2[comp2$clade=="Pflu",1],side=1,col="darkorange2",lwd=1)
rug(comp2[comp2$clade=="Psyr",5],side=2,col="darkslateblue",lwd=1)
rug(comp2[comp2$clade=="Pflu",5],side=2,col="darkorange2",lwd=1)
for (i in 1:length(comp2[,1]))
{
  if (comp2$clade[i]=="Pflu") { lines(c(comp2[i,1],0),c(0,comp2[i,5]),col="darkorange2") }
  else { lines(c(comp2[i,1],0),c(0,comp2[i,5]),col="darkslateblue",lty="dashed") }
}

# second for L_MM v L_KB
plot(comp2[,2],comp2[,6],pch="",ylim=c(0,500),xlim=c(0,3000), cex.axis=0.75,ylab="",xlab="",las=2,bty="n")
rug(comp2[comp2$clade=="Psyr",2],side=1,col="darkslateblue",lwd=1)
rug(comp2[comp2$clade=="Pflu",2],side=1,col="darkorange2",lwd=1)
rug(comp2[comp2$clade=="Psyr",6],side=2,col="darkslateblue",lwd=1)
rug(comp2[comp2$clade=="Pflu",6],side=2,col="darkorange2",lwd=1)
for (i in 1:length(comp2[,2]))
{
  if (comp2$clade[i]=="Pflu") { lines(c(comp2[i,2],0),c(0,comp2[i,6]),col="darkorange2") }
  else { lines(c(comp2[i,2],0),c(0,comp2[i,6]),col="darkslateblue",lty="dashed") }
}

# third for K_MM v K_KB
plot(comp2[,3],comp2[,7],pch="",ylim=c(0,1.2),xlim=c(0,1.0), cex.axis=0.75,ylab="",xlab="",las=2,bty="n")
rug(comp2[comp2$clade=="Psyr",3],side=1,col="darkslateblue",lwd=1)
rug(comp2[comp2$clade=="Pflu",3],side=1,col="darkorange2",lwd=1)
rug(comp2[comp2$clade=="Psyr",7],side=2,col="darkslateblue",lwd=1)
rug(comp2[comp2$clade=="Pflu",7],side=2,col="darkorange2",lwd=1)
for (i in 1:length(comp2[,3]))
{
  if (comp2$clade[i]=="Pflu") { lines(c(comp2[i,3],0),c(0,comp2[i,7]),col="darkorange2") }
  else { lines(c(comp2[i,3],0),c(0,comp2[i,7]),col="darkslateblue") }
}

# now spit out correlations between traits for each of the clades, to add to right upper corner of plots:

with(subset(comp2,comp2$clade=="Psyr"), cor.test(mu.mm,mu.kb))
with(subset(comp2,comp2$clade=="Pflu"), cor.test(mu.mm,mu.kb))

with(subset(comp2,comp2$clade=="Psyr"), cor.test(lambda.mm,lambda.kb))
with(subset(comp2,comp2$clade=="Pflu"), cor.test(lambda.mm,lambda.kb))

with(subset(comp2,comp2$clade=="Psyr"), cor.test(A.mm,A.kb))
with(subset(comp2,comp2$clade=="Pflu"), cor.test(A.mm,A.kb))

# second part: correlations:
KBcomp <- with(comp2,cbind(mu.kb,lambda.kb,A.kb))
KBcomp1 <- KBcomp[-c(27:40),] # for P. syr
KBcomp2 <- KBcomp[-c(1:26),] # for P. flu

library(psych)

KBcPsyr <- corr.test(KBcomp1,method="pearson",adjust="fdr")
KBcPflu <- corr.test(KBcomp2,method="pearson",adjust="fdr")

# write tables for deleting non-sig p-values..
write.table(KBcPsyr$p,file="KBcPsyr_p.txt",sep="\t")
write.table(KBcPflu$p,file="KBcPflu_p.txt",sep="\t")

# now spit out heatmaps..
palette.breaks <- seq(-1,1,0.125)
#color.palette  <- colorRampPalette(c("midnightblue","white","darkorange2"))(length(palette.breaks) - 1)
color.palette  <- colorRampPalette(c("midnightblue","white","darkorange2"))(length(palette.breaks) - 1)
h1 <- heatmap.2(KBcPsyr$r,dendrogram = "none",scale = "none", trace = "none", key = TRUE, keysize = 1.5, density.info=c("none"), labRow = NA,labCol = NA,
                Rowv       = FALSE,
                Colv       = FALSE,  
                col    = color.palette,
                #col     = cm.colors,
                breaks = palette.breaks,
                symbreaks = TRUE,
                cellnote = round(KBcPsyr$r,2),
                notecex  = 0.9,
                notecol  = "black",
                colsep,
                rowsep,
                sepcolor="white",
                sepwidth=c(0.1,0.1)
)

h2 <- heatmap.2(KBcPflu$r,dendrogram = "none",scale = "none", trace = "none", key = TRUE, keysize = 1.5, density.info=c("none"), labRow = NA,labCol = NA,
                Rowv       = FALSE,
                Colv       = FALSE,  
                col    = color.palette,
                #col     = cm.colors,
                breaks = palette.breaks,
                symbreaks = TRUE,
                cellnote = round(KBcPflu$r,2),
                notecex  = 0.9,
                notecol  = "black",
                colsep,
                rowsep,
                sepcolor="white",
                sepwidth=c(0.1,0.1)
)


#### SUPPLEMENTAL FIGURE 3: COMPETITION PHYLOGENY PAIRWISE BONANZA ####
# this was all done in Illusrator; no scripts required.

#### SUPPLEMENTAL FIGURES 4 & 5: ALL TRAIT CORRELATIONS AS SCATTERPLOTS BY CLADE ####
# do without marginal density plots since these are already in the main text as Fig 1B.

clade_cols <- as.character(c("white","black"))

r_l <- ggplot(comp2, aes(x=mu.mm, y=lambda.mm)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(0,3000)) +
  theme_bw() + theme(legend.position = "none",
                     legend.justification = "none",
                     panel.background = element_blank(),
                     axis.title.x = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank())

r_k <- ggplot(comp2, aes(x=mu.mm, y=A.mm)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(0,0.8)) +
  theme_bw() + theme(legend.position = "none",
                     legend.justification = "none",
                     panel.background = element_blank(),
                     axis.title.x = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank())

l_k <- ggplot(comp2, aes(x=lambda.mm, y=A.mm)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,3000)) + scale_y_continuous(limits=c(0,0.8)) +
  theme_bw() + theme(legend.position = "none",
                     legend.justification = "none",
                     panel.background = element_blank(),
                     axis.title.x = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank())

#grid.arrange(r_l,r_k,l_k,ncol=3) # make each plot 250 x 250 pixels

## ROW 2
r_ci <- ggplot(comp2, aes(x=mu.mm, y=c.i)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

l_ci <- ggplot(comp2, aes(x=lambda.mm, y=c.i)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,3000)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

k_ci <- ggplot(comp2, aes(x=A.mm, y=c.i)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.8)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

#grid.arrange(r_ci,l_ci,k_ci,ncol=3)
## ROW 3
r_ch <- ggplot(comp2, aes(x=mu.mm, y=c.h)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

l_ch <- ggplot(comp2, aes(x=lambda.mm, y=c.h)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,3000)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

k_ch <- ggplot(comp2, aes(x=A.mm, y=c.h)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.8)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

#grid.arrange(r_ch,l_ch,k_ch,ncol=3)

## ROW 4
r_cw <- ggplot(comp2, aes(x=mu.mm, y=c.w)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(-1,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

l_cw <- ggplot(comp2, aes(x=lambda.mm, y=c.w)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,3000)) + scale_y_continuous(limits=c(-1,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

k_cw <- ggplot(comp2, aes(x=A.mm, y=c.w)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.8)) + scale_y_continuous(limits=c(-1,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

#grid.arrange(r_cw,l_cw,k_cw,ncol=3)

## ROW 5: correlations within defense.. 
ci_ch <- ggplot(comp2, aes(x=c.i, y=c.h)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

cw_ci <- ggplot(comp2, aes(x=c.w, y=c.i)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

cw_ch <- ggplot(comp2, aes(x=c.w, y=c.h)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

#grid.arrange(ci_ch,cw_ci,cw_ch,ncol=3)

## ROW 6
cw_cr <- ggplot(comp2, aes(x=c.w, y=c.r)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

cw_ct <- ggplot(subset(comp2,comp2$clade=="Pflu"), aes(x=c.w, y=c.t)) + geom_point(size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = c("white")) + scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

ct_cr <- ggplot(subset(comp2,comp2$clade=="Pflu"), aes(x=c.t, y=c.r)) + geom_point(size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = c("white")) + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) +
  theme_bw() + theme(legend.position = "none",legend.justification = "none",panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

#grid.arrange(cw_cr,cw_ct,ct_cr,ncol=3)

grid.arrange(r_l,r_k,l_k,r_ci,l_ci,k_ci,r_ch,l_ch,k_ch,r_cw,l_cw,k_cw,ci_ch,cw_ci,cw_ch,cw_cr,cw_ct,ct_cr,ncol=3,nrow=6)



# OK now I need to spit out the axes.. use same plots as above, but chuck the blank elements:
r_l2 <- ggplot(comp2, aes(x=mu.mm, y=lambda.mm)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(0,3000)) +
  theme_bw() + theme(legend.position = "none",
                     legend.justification = "none",
                     panel.background = element_blank())

r_k2 <- ggplot(comp2, aes(x=mu.mm, y=A.mm)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,0.001)) + scale_y_continuous(limits=c(0,0.8)) +
  theme_bw() + theme(legend.position = "none",
                     legend.justification = "none",
                     panel.background = element_blank())

l_k2 <- ggplot(comp2, aes(x=lambda.mm, y=A.mm)) + geom_point(aes(fill=comp2$clade),size=2,pch=21,colour="black") +
  stat_smooth(method=lm,se=TRUE, aes(lty=paste(as.numeric(comp2$clade)))) + theme(legend.position = "none") +
  scale_shape_manual(values=c(21,21)) +
  scale_fill_manual(values = clade_cols) + scale_x_continuous(limits=c(0,3000)) + scale_y_continuous(limits=c(0,0.8)) +
  theme_bw() + theme(legend.position = "none",
                     legend.justification = "none",
                     panel.background = element_blank())
grid.arrange(r_l2,r_k2,l_k2,ncol=3)


#### polynomial regression for P. syringae relationship between offense and defense.. ####

fit1 <- lm(sample1$Population ~ sample1$Year)
fit2 <- lm(sample1$Population ~ sample1$Year + I(sample1$Year^2))


pn1 <- lm(c.h~c.i,data=subset(comp2,comp2$clade=="Psyr"))
pn2 <- lm(c.h~c.i + I(c.i^2),data=subset(comp2,comp2$clade=="Psyr"))
pn3 <- lm(c.h~poly(c.i,3,raw=TRUE),data=subset(comp2,comp2$clade=="Psyr"))

anova(pn3,pn2)
anova(pn2,pn1)

with(subset(comp2,comp2$clade=="Psyr"),plot(c.h~c.i,xlim=c(0,1),ylim=c(0,1)))
with(subset(comp2,comp2$clade=="Psyr"),points(c.h~predict(pn2),type="l"))     



ggplot(data=comp2,aes(x=c.i,y=c.h,col=clade)) + geom_point(aes(group=clade,size=2)) + theme_bw() + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1))
