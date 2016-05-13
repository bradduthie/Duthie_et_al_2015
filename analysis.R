# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# The subsequent code replicates the analyses used in Duthie et al. (2015) American Naturalist. All data were collected in 2010 in the field from populations of Sonoran Desert rock fig (Ficus petiolaris) trees. For this analysis, three files are needed:

# The first file includes counts of all female non-pollinating fig wasps from the syconia of the Sonoran Desert rock fig, Ficus petiolaris (`wasp_data.csv'). The data file includes the site, tree and syconia (`Fruit') labels from each syconia sampled (rows); the data also includes columns indicating tree lat-lon coordinates, pollinating foundresses arriving to the syconia (corpses remain within syconia), fruit volume, counts of all wasp species, and the number of neighbours a tree has within a 1 km radius. The pollinator (`Poll') is an unnamed (as of 2015) species of Pegoscapus. Wasps `LO1', `SO1', and `SO2' are all unnamed species of the genus Idarnes, and wasps `Het1' and `Het2' are both unnamed species of the genus Heterandrium. The wasps `Physo' and `Aepoc' are wasps of unnamed species of the genera Physothorax and Aepocerus, respectively; `NA' denotes missing values.

# The second file (`wing_loadings.csv') includes measurements to estimate wing loading for 83 wasps (rows) of the species of interest in Duthie et al. (2015): LO1, SO1, SO2, Het1, and Het2. Each row is a single wasps, and columns show the site, tree, and fruit from which the wasp was sampled. Lengths and widths of wasp heads, thoraxes, and abdomens are included as measured at their widest points (e.g., for the wasp abdomen, the whole length of the segment is reported, along with an estimate of its width at the widest point). Wing areas were calculated using wing images and ImageJ software.

# The third file (`egg_loads.csv') includes estimates of wasp egg loads from 54 wasps (rows) of the following species:  LO1, SO1, SO2, Het1, and Het2. Both mature and immature egg counts were estimated, and columns include the site, tree, and fruit from which the wasp was sampled.

# The code below will replicate the statistical analysis of Duthie et al. (2015). This analysis includes: 1) The estimation of a colonization index for each species, and the correlation of this index with species mean egg load. 2) The correlation between mean species wing loading and mean species fecundity. 3) The correlation betwen species wing loading and the effect (regression slope) that the number of neighbouring fig trees within a 1 km radius of the tree from which wasps were sampled had on the species abundance; the same correlation is estimated with egg load instead of wing loading. Finally, 4) the differences between species fecundities/wing loadings as correlated with species among-tree density correlations (i.e., individual points are absolute differences in fecundities or wing loadings, and are compared with individual points that are correlations between the densities of species on a fig tree).

# After the analyses are replicated, the remaining code re-builds the figures presented in Duthie et al., 2015. 

# Any enquiries about these data or the analysis and code that follows can be made to Brad Duthie (aduthie@abdn.ac.uk; brad.duthie@gmail.com).
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #


setwd(path.expand("~/"));

# Read files of wasp abundances, wing loadings, and egg counts
wasps <- read.table(file="wasp_data.csv",header=TRUE,sep="\t");
wings <- read.table(file="wing_loadings.csv",header=TRUE,sep="\t");
eggs  <- read.table(file="egg_loads.csv",header=TRUE,sep="\t");

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
#  XXX XXX ANALYSES ARE SHOWN BELOW, FOLLOWED BY FIGURES XXX XXX  #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# Calculate mean, standard deviation, and sample sizes for eggs
eggs.means <- tapply(X=eggs[,6],INDEX=eggs[,1],FUN=mean);
eggs.stdev <- tapply(X=eggs[,6],INDEX=eggs[,1],FUN=sd);
eggs.sampl <- tapply(X=eggs[,6],INDEX=eggs[,1],FUN=length);

# ---------------------------------------------------------------
# XXX CORRELATION BETWEEN COLONIZATION INDEX AND FECUNDITY  XXX #
# ---------------------------------------------------------------
# Calculate the colonization index for each species:
# missing identifies fruits with needed values that are missing
missing  <- apply(X=wasps[,7:15],MARGIN=1,FUN=sum);
dat      <- wasps[-which(is.na(missing)),]; 
Allwasps <- as.vector(apply(X=dat[,8:15],MARGIN=1,FUN=sum));
Allnpoll <- as.vector(apply(X=dat[,9:15],MARGIN=1,FUN=sum));
# Mean densities of each non-pollinator calculated below
SO1.den  <- tapply(X=dat[,9], INDEX=dat[,2], FUN=mean);
SO2.den  <- tapply(X=dat[,10],INDEX=dat[,2], FUN=mean);
LO1.den  <- tapply(X=dat[,11],INDEX=dat[,2], FUN=mean);
Het1.den <- tapply(X=dat[,12],INDEX=dat[,2], FUN=mean);
Het2.den <- tapply(X=dat[,13],INDEX=dat[,2], FUN=mean);
# Get the mean wasp count and syconia volume per tree below
Tree.wsp <- tapply(X=Allwasps, INDEX=dat[,2], FUN=mean);
Tree.vol <- tapply(X=dat[,7] , INDEX=dat[,2], FUN=mean);
# Below gets the colonization difficult (wasps/volume on tree)
Col.diff <- Tree.wsp / Tree.vol;
# Below gets the colonization index for each species
SO1.sl   <- -1*as.numeric(lm(SO1.den ~Col.diff)$coefficients[2]);
SO2.sl   <- -1*as.numeric(lm(SO2.den ~Col.diff)$coefficients[2]);
LO1.sl   <- -1*as.numeric(lm(LO1.den ~Col.diff)$coefficients[2]);
Het1.sl  <- -1*as.numeric(lm(Het1.den~Col.diff)$coefficients[2]);
Het2.sl  <- -1*as.numeric(lm(Het2.den~Col.diff)$coefficients[2]);
# Combines this into a vector equivalent to egg.means;
CI.means <- c(Het1.sl,Het2.sl,LO1.sl,SO1.sl,SO2.sl);

# Test the correlation between colonization index and fecundity;
cor.test(CI.means, eggs.means); # P = 0.0022; R^2 = 0.9702
# ---------------------------------------------------------------

# ---------------------------------------------------------------
#  XXX XXX  CORRELATION BETWEEN WING LOADING AND FECUNDITY  XXX #
# ---------------------------------------------------------------
# Calculate mean & sd for wing wing length, area, & wing loading:
# ---------------------------------------------------------------
# (Length+Width)/4 gives the radius for the sphere
sphere.vol <- function(Length,Width) 4/3*3.14*(((Length+Width)/4)^3);
# Length refers to segment length, width as segment width
ellips.vol <- function(Length,Width) 4/3*3.14*(Length/2)*(Width/2)^2;
# Length refers to ovipositor length, width as ovipositor width
cylind.vol <- function(Length,Width) 3.14*((0.5*Width)^2)*Length;
# Using each to get body segment volumes:
Head.vols  <- sphere.vol(Length=wings[,5], Width=wings[,6]);
Thor.vols  <- ellips.vol(Length=wings[,7], Width=wings[,8]);
Abdo.vols  <- ellips.vol(Length=wings[,9], Width=wings[,10]);
Ovip.vols  <- cylind.vol(Length=wings[,11],Width=wings[,12]);
wasp.vols  <- Head.vols + Thor.vols + Abdo.vols +  Ovip.vols;

# Now get total wing area (Forwing + Hindwing times two);
wasp.wing  <- (wings[,13] + wings[,14]) * 2;
sp.wasp.wi <- tapply(X=wasp.wing, INDEX=wings[,1],FUN=mean);
sp.wasp.sd <- tapply(X=wasp.wing, INDEX=wings[,1],FUN=sd);
sp.wasp.le <- tapply(X=wasp.wing, INDEX=wings[,1],FUN=length);
sp.wasp.se <- sp.wasp.sd / sqrt(sp.wasp.le);
wasp.WL    <- wasp.vols / wasp.wing;

# Mean, sd, and sample size of species-specific Volume, WL, length
sp.vol.mu  <- tapply(X=wasp.vols, INDEX=wings[,1],FUN=mean);
sp.vol.sd  <- tapply(X=wasp.vols, INDEX=wings[,1],FUN=sd);
sp.vol.sa  <- tapply(X=wasp.vols, INDEX=wings[,1],FUN=length);
sp.WL.mu   <- tapply(X=wasp.WL,   INDEX=wings[,1],FUN=mean);
sp.WL.sd   <- tapply(X=wasp.WL,   INDEX=wings[,1],FUN=sd);
sp.WL.sa   <- tapply(X=wasp.WL,   INDEX=wings[,1],FUN=length);
sp.wlen.mu <- tapply(X=wings[,15],INDEX=wings[,1],FUN=mean);
sp.wlen.sd <- tapply(X=wings[,15],INDEX=wings[,1],FUN=sd);
sp.wlen.sa <- tapply(X=wings[,15],INDEX=wings[,1],FUN=length);

# Test the correlation between wasp wing loading and fecundity:
cor.test(eggs.means,sp.WL.mu); # P = 0.0387; R^2 = 0.8058
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# XXX WING LOADING, EGG LOAD & EFFECT OF NEIGHBOURING TREES XXX #
# ---------------------------------------------------------------

# Below gets vector for number of 1 km neighbours a tree has
tX1km <- as.vector(tapply(X=dat[,16],INDEX=dat[,2],FUN=mean));
# Species-specific regression of densities against tree neighbours
rSO1  <- lm((dat[,9]/dat[,7]) ~dat[,16])$coefficients[2];
rSO2  <- lm((dat[,10]/dat[,7])~dat[,16])$coefficients[2];
rLO1  <- lm((dat[,11]/dat[,7])~dat[,16])$coefficients[2];
rHet1 <- lm((dat[,12]/dat[,7])~dat[,16])$coefficients[2];
rHet2 <- lm((dat[,13]/dat[,7])~dat[,16])$coefficients[2];
# Orders the species in the same way as eggs.means & sp.WL.mu
respo <- c(rHet1,rHet2,rLO1,rSO1,rSO2);
# Test the correlation between species-specific regression of
# densities against tree neighbours versus mean egg loads:
cor.test(respo,eggs.means); # P = 0.0371; R^2 = 0.8111
# Test the correlation between species-specific regression of
# densities against tree neighbours versus mean wing loadings:
cor.test(respo,sp.WL.mu); # P = 0.0379; R^2 = 0.8085
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# XXX DIFFERENCE IN FECUNDITY VERSUS WASP DENSITY CORRELATIONS  #
# ---------------------------------------------------------------
# This gets the differences between all species fecundities:
DiffFecu <- c( abs(eggs.means[1]-eggs.means[2]), 
               abs(eggs.means[1]-eggs.means[3]), 
               abs(eggs.means[1]-eggs.means[4]),
               abs(eggs.means[1]-eggs.means[5]), 
               abs(eggs.means[2]-eggs.means[3]), 
               abs(eggs.means[2]-eggs.means[4]),
               abs(eggs.means[2]-eggs.means[5]), 
               abs(eggs.means[3]-eggs.means[4]), 
               abs(eggs.means[3]-eggs.means[5]),
               abs(eggs.means[4]-eggs.means[5]));

# This gets the differences between all species wing loadings:
DiffWL  <-  c(  abs(sp.WL.mu[1]-sp.WL.mu[2]), 
                abs(sp.WL.mu[1]-sp.WL.mu[3]),
                abs(sp.WL.mu[1]-sp.WL.mu[4]),
                abs(sp.WL.mu[1]-sp.WL.mu[5]), 
                abs(sp.WL.mu[2]-sp.WL.mu[3]),
                abs(sp.WL.mu[2]-sp.WL.mu[4]),
                abs(sp.WL.mu[2]-sp.WL.mu[5]),
                abs(sp.WL.mu[3]-sp.WL.mu[4]),
                abs(sp.WL.mu[3]-sp.WL.mu[5]),
                abs(sp.WL.mu[4]-sp.WL.mu[5]));

# This gets the differences between all species densities:
# Note: ordered in the same way as the egg means
dHet1    <- tapply(X=(dat[,12]/dat[,7]),INDEX=dat[,2],FUN=mean);
dHet2    <- tapply(X=(dat[,13]/dat[,7]),INDEX=dat[,2],FUN=mean);
dLO1     <- tapply(X=(dat[,11]/dat[,7]),INDEX=dat[,2],FUN=mean);
dSO1     <- tapply(X=(dat[,9] /dat[,7]),INDEX=dat[,2],FUN=mean);
dSO2     <- tapply(X=(dat[,10]/dat[,7]),INDEX=dat[,2],FUN=mean);

# This gets the correlation between species densities
DensCor  <- c(  cor(dHet1, dHet2), cor(dHet1, dLO1),
                cor(dHet1, dSO1),  cor(dHet1, dSO2), 
                cor(dHet2, dLO1),  cor(dHet2, dSO1), 
                cor(dHet2, dSO2),  cor(dLO1,  dSO1), 
                cor(dLO1,  dSO2),  cor(dSO1,  dSO2)); 

# Density correlations against differences in egg load:
cor.test(DensCor,DiffFecu); # P =  0.3773; R^2 = 0.09845
# Density correlations against differences in wing loading:
cor.test(DensCor,DiffWL); # P =  0.0462; R^2 = 0.4097


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
#    XXX XXX FIGURES ARE MADE AS IN DUTHIE ET AL. 2015   XXX XXX  #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
library("extrafont");
loadfonts(device="postscript");

eggs.se  <- eggs.stdev / sqrt(eggs.sampl); # SE for egg means
# Below gets the standard errors for the colonization index
SO1.se   <- as.numeric(summary(lm(SO1.den ~Col.diff))$coefficients[2,2]);
SO2.se   <- as.numeric(summary(lm(SO2.den ~Col.diff))$coefficients[2,2]);
LO1.se   <- as.numeric(summary(lm(LO1.den ~Col.diff))$coefficients[2,2]);
Het1.se  <- as.numeric(summary(lm(Het1.den~Col.diff))$coefficients[2,2]);
Het2.se  <- as.numeric(summary(lm(Het2.den~Col.diff))$coefficients[2,2]);
CI.se    <- c(Het1.se, Het2.se, LO1.se, SO1.se, SO2.se);
# Below gets the standard errors for the effects of neighbouring density
rSO1.se  <- as.numeric(summary(lm((dat[,9]/dat[,7]) ~dat[,16]))$coefficients[2,2]);
rSO2.se  <- as.numeric(summary(lm((dat[,10]/dat[,7])~dat[,16]))$coefficients[2,2]);
rLO1.se  <- as.numeric(summary(lm((dat[,11]/dat[,7])~dat[,16]))$coefficients[2,2]);
rHet1.se <- as.numeric(summary(lm((dat[,12]/dat[,7])~dat[,16]))$coefficients[2,2]);
rHet2.se <- as.numeric(summary(lm((dat[,13]/dat[,7])~dat[,16]))$coefficients[2,2]);
r.se     <- c(rHet1.se, rHet2.se, rLO1.se, rSO1.se, rSO2.se);

# ---------------------------------------------------------------
# XXX               CODE TO REPLICATE FIGURE 1              XXX #
# ---------------------------------------------------------------
setEPS();
postscript("Figure1.eps",width=12,height=7,family="Arial");
par(mar=c(6,5,1,1),mfrow=c(1,2))
# ------------ Subfigure `a'
plot(x=CI.means,y=eggs.means,pch=20,xlim=c(-200,125),xaxt="n",cex.axis=1.75,
	ylab="Egg load",xlab="Colonization index",cex.lab=2,cex=1.5);
axis(side=1,at=c(-200,-100,0,100),cex.axis=1.75);
for(i in 1:length(eggs.means)){
	arrows( x0=CI.means[i],x1=CI.means[i],
            y0=eggs.means[i]-eggs.se[i],
            y1=eggs.means[i]+eggs.se[i],
	        angle=90,code=3,length=0.05,lwd=2);
	arrows( x0=CI.means[i]-CI.se[i],
            x1=CI.means[i]+CI.se[i],
            y0=eggs.means[i],y1=eggs.means[i],
            angle=90,code=3,length=0.05,lwd=2);
}
B0 <- lm(eggs.means~CI.means)$coefficients[1];
B1 <- lm(eggs.means~CI.means)$coefficients[2];
x  <- seq(from=-123,to=65,by=1);
y  <- B0 + B1*x;
points(x=x,y=y,type="l",lwd=1.5,lty="dashed");
text(y=eggs.means[1],x=CI.means[1]-85,labels=expression(Het1),cex=1.5);
text(y=eggs.means[2],x=CI.means[2]+35,labels=expression(Het2),cex=1.5);
text(y=eggs.means[3],x=CI.means[3]+95,labels=expression(LO1),cex=1.5);
text(y=eggs.means[4],x=CI.means[4]-50,labels=expression(SO1),cex=1.5);
text(y=eggs.means[5],x=CI.means[5]-90,labels=expression(SO2),cex=1.5);
text(y=eggs.means[3]-1,x=117,labels="a",cex=3);
# ------------ Subfigure `b'
sp.WL.se <- sp.WL.sd / sqrt(sp.WL.sa);
par(mar=c(6,1,1,2.5));
plot(x=sp.WL.mu,y=eggs.means,pch=20,xlim=c(0.125,0.175),xaxt="n",
	ylab="",xlab=expression(paste("Wing loading (",mm^3/mm^2,")")),
	yaxt="n",cex.lab=2,cex.axis=1.75,cex=1.5);
axis(side=1,at=c(0.13,0.15,0.17),cex.axis=1.75);

for(i in 1:length(eggs.means)){
	arrows(x0=sp.WL.mu[i],x1=sp.WL.mu[i],
           y0=eggs.means[i]-eggs.se[i],
           y1=eggs.means[i]+eggs.se[i],
           angle=90,code=3,length=0.05,lwd=2);
	arrows(x0=sp.WL.mu[i]-sp.WL.se[i],
           x1=sp.WL.mu[i]+sp.WL.se[i],
           y0=eggs.means[i],y1=eggs.means[i],
           angle=90,code=3,length=0.05,lwd=2)
}
B0 <- lm(eggs.means~sp.WL.mu)$coefficients[1];
B1 <- lm(eggs.means~sp.WL.mu)$coefficients[2];
x  <- seq(from=0.125,to=0.175,by=0.0001);
y  <- B0 + B1*x;
points(x=x,y=y,type="l",lwd=1.5,lty="dashed");
text(y=eggs.means[1],x=sp.WL.mu[1]+0.011,labels=expression(Het1),cex=1.5);
text(y=eggs.means[2],x=sp.WL.mu[2]+0.0085,labels=expression(Het2),cex=1.5);
text(y=eggs.means[3],x=sp.WL.mu[3]-0.012,labels=expression(LO1),cex=1.5);
text(y=eggs.means[4],x=sp.WL.mu[4]+0.011,labels=expression(SO1),cex=1.5);
text(y=eggs.means[5],x=sp.WL.mu[5]+0.011,labels=expression(SO2),cex=1.5);
text(y=eggs.means[3]-1,x=0.126,labels="b",cex=3);
dev.off();
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# XXX               CODE TO REPLICATE FIGURE 2              XXX #
# ---------------------------------------------------------------
sp.WL.se <- sp.WL.sd / sqrt(sp.WL.sa);
setEPS();
postscript("Figure2.eps",width=12,height=7,family="Arial");
par(mar=c(6,5,1,1),mfrow=c(1,2))

plot(x=sp.WL.mu,y=respo,pch=20,xlim=c(0.125,0.175),
    ylab=expression(paste("Effect of neighboring trees on wasp density")),
    xlab=expression(paste("Wing loading (",mm^3/mm^2,")")),
    cex.lab=1.75,cex.axis=1.5,
    ylim=c((respo[1]-r.se[1]),(respo[3]+r.se[3])),cex=1.5);

for(i in 1:length(sp.WL.mu)){
	arrows(x0=sp.WL.mu[i],x1=sp.WL.mu[i],
           y0=respo[i]-r.se[i],
           y1=respo[i]+r.se[i],
	angle=90,code=3,length=0.05,lwd=2);
	arrows(x0=sp.WL.mu[i]-sp.WL.se[i],
           x1=sp.WL.mu[i]+sp.WL.se[i],
           y0=respo[i],y1=respo[i],
	angle=90,code=3,length=0.05,lwd=2);
}
B0 <- lm(respo~sp.WL.mu)$coefficients[1];
B1 <- lm(respo~sp.WL.mu)$coefficients[2];
x  <- seq(from=0.125,to=0.175,by=0.0001);
y  <- B0 + B1*x;
points(x=x,y=y,type="l",lwd=1.5,lty="dashed");
text(y=respo[1],x=sp.WL.mu[1]+sp.WL.se[1]+0.0025,labels=expression(Het1),cex=1.125);
text(y=respo[2],x=sp.WL.mu[2]+sp.WL.se[2]+0.0026,labels=expression(Het2),cex=1.125);
text(y=respo[3],x=sp.WL.mu[3]-sp.WL.se[3]-0.0021,labels=expression(LO1),cex=1.125);
text(y=respo[4],x=sp.WL.mu[4]-sp.WL.se[4]-0.0025,labels=expression(SO1),cex=1.125);
text(y=respo[5],x=sp.WL.mu[5]+sp.WL.se[5]+0.0025,labels=expression(SO2),cex=1.125);
text(y=respo[1]-1,x=0.80,labels="b",cex=2);
text(y=3.9e-05,x=0.126,cex=2.5,labels="a");

par(mar=c(6,1,1,3));
plot(x=eggs.means,y=respo,pch=20,xlim=c(30,80),
    ylab="",xlab=expression(paste("Egg load")),yaxt="n",
    cex.lab=1.5,cex.axis=1.5,
    ylim=c((respo[1]-r.se[1]),(respo[3]+r.se[3])),cex=1.5);

for(i in 1:length(eggs.means)){
    arrows(x0=eggs.means[i],x1=eggs.means[i],
           y0=respo[i]-r.se[i],
           y1=respo[i]+r.se[i],
           angle=90,code=3,length=0.05,lwd=2);
	arrows(x0=eggs.means[i]-eggs.se[i],
           x1=eggs.means[i]+eggs.se[i],
           y0=respo[i],y1=respo[i],
           angle=90,code=3,length=0.05,lwd=2);
}
B0 <- lm(respo~eggs.means)$coefficients[1];
B1 <- lm(respo~eggs.means)$coefficients[2];
x  <- seq(from=30,to=78,by=0.001);
y  <- B0 + B1*x;
points(x=x,y=y,type="l",lwd=1.5,lty="dashed");
text(y=respo[1],x=eggs.means[1]+eggs.se[1]+3,labels=expression(Het1),cex=1.125);
text(y=respo[2],x=eggs.means[2]+eggs.se[2]+2.5,labels=expression(Het2),cex=1.125);
text(y=respo[3],x=eggs.means[3]-eggs.se[3]-2.25,labels=expression(LO1),cex=1.125);
text(y=respo[4],x=eggs.means[4]-eggs.se[4]-2.5,labels=expression(SO1),cex=1.125);
text(y=respo[5],x=eggs.means[5]+eggs.se[5]+2.5,labels=expression(SO2),cex=1.125);
text(y=respo[1]-1,x=0.80,labels="b",cex=2);
text(y=3.9e-05,x=30.5,cex=2.5,labels="b");

dev.off();
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# XXX               CODE TO REPLICATE A2                    XXX #
# ---------------------------------------------------------------

setEPS();
postscript("FigureA2.eps",width=12,height=8,family="Arial");
par(mar=c(6,6,1,0.5),mfrow=c(1,2));
plot(x=DiffFecu,y=DensCor,pch=20,cex=1.75,
	xlab="Species difference in egg load",
	ylab="Correlation of species density",
	cex.lab=1.5,cex.axis=1.75);
text(x=41,y=0.79,labels="a",cex=3.25);
par(mar=c(6,1,1,1));
plot(x=DiffWL,y=DensCor,pch=20,cex=1.75,
	xlab=expression(paste("Species difference in wing loading (",mm^3/mm^2,")")),
	ylab="",yaxt="n",
	cex.lab=1.5,cex.axis=1.75);
text(x=0.0295,y=0.79,labels="b",cex=3.25);

dev.off();
# ---------------------------------------------------------------







