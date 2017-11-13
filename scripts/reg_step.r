# @author:    joao lopes
# @workplace: University of Reading
# @date:      13th May 2009

# Function to perform the regression step on a file with rejection step results.
# It also plots histograms of the obtained posterior distributions.
# Creates .eps image files in the directory from where it is run

# @arg data_file - file with summary of 'real' data (.trg)
# @arg rej_file  - file with rejection step results (.rej)
# @arg pri_file  - file with sample of priors (.pri)
reg_step <- function(data_file , rej_file , pri_file){

	#import locfit library
	library(locfit)

	#import Mark Beaumont's scripts
	source("make_pd2.r")
	source("loc2plot_d.r")

	#demographic parameters' priors (minimum and maximum values)
	mintev <- 0
	maxtev <- 10000
	minNe1 <- 0
	maxNe1 <- 1000
	minNe2 <- 0
	maxNe2 <- 2000
	minNeA <- 0
	maxNeA <- 1500
	minmig1 <- 0
	maxmig1 <- 0.0001
	minmig2 <- 0
	maxmig2 <- 0.0001

	#import the .rej file
	abc.rej <- data.matrix(read.table(rej_file))

	#import the .pri files
	priors <- data.matrix(read.table(pri_file))

	#import the .trg files
	target <- data.matrix(read.table(data_file))

	# Plot line for sequence DNA mutation rate and save it in a .eps file
	abc.reg.mut <- makepd4(target,abc.rej[,5],abc.rej[,17:34],tol=1,rej=F,transf="none")
	plot(locfit(~abc.reg.mut$x,weight=abc.reg.mut$wt),main="prior (black) and posterior (blue) distributions",xlab="mut rate",col="blue")
	plot(locfit(~priors[,5]),col="black",add=T)
	dev.copy2eps(file="mutrate_reg.eps", horizontal=F)
	print("mut rate done.")

	# Plot line for sequence DNA recombination rate and save it in a .eps file
	abc.reg.rec <- makepd4(target,abc.rej[,9],abc.rej[,17:34],tol=1,rej=F,transf="none")
	plot(locfit(~abc.reg.rec$x,weight=abc.reg.rec$wt),main="prior (black) and posterior (blue) distributions",xlab="rec rate",col="blue")
	plot(locfit(~priors[,9]),col="black",add=T)
	dev.copy2eps(file="recrate_reg.eps", horizontal=F)
	print("rec rate done.")

	# Plot line for splitting time and save it in a .eps file
	abc.reg.tev <- makepd4(target,abc.rej[,11],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(mintev,maxtev))
	plot(locfit(~abc.reg.tev$x,weight=abc.reg.tev$wt,xlim=c(mintev,maxtev)),main="prior (black) and posterior (blue) distributions",xlab="tev",col="blue")
	plot(locfit(~priors[,11],xlim=c(mintev,maxtev)),col="black",add=T)
	dev.copy2eps(file="tev_reg.eps", horizontal=F)
	print("tev done.")

	# Plot line for effective size of population 1 and save it in a .eps file
	abc.reg.ne1 <- makepd4(target,abc.rej[,12],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minNe1,maxNe1))
	plot(locfit(~abc.reg.ne1$x,weight=abc.reg.ne1$wt,xlim=c(minNe1,maxNe1)),main="prior (black) and posterior (blue) distributions",xlab="Ne1",col="blue")
	plot(locfit(~priors[,12],xlim=c(minNe1,maxNe1)),col="black",add=T)
	dev.copy2eps(file="Ne1_reg.eps", horizontal=F)
	print("Ne1 done.")

	# Plot line for effective size of population 2 and save it in a .eps file
	abc.reg.ne2 <- makepd4(target,abc.rej[,13],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minNe2,maxNe2))
	plot(locfit(~abc.reg.ne2$x,weight=abc.reg.ne2$wt,xlim=c(minNe2,maxNe2)),main="prior (black) and posterior (blue) distributions",xlab="Ne2",col="blue")
	plot(locfit(~priors[,13],xlim=c(minNe2,maxNe2)),col="black",add=T)
	dev.copy2eps(file="Ne2_reg.eps", horizontal=F)
	print("Ne2 done.")

	# Plot line for effective size of ancestor population and save it in a .eps file
	abc.reg.nea <- makepd4(target,abc.rej[,14],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minNeA,maxNeA))
	plot(locfit(~abc.reg.nea$x,weight=abc.reg.nea$wt,xlim=c(minNeA,maxNeA)),main="prior (black) and posterior (blue) distributions",xlab="NeA",col="blue")
	plot(locfit(~priors[,14],xlim=c(minNeA,maxNeA)),col="black",add=T)
	dev.copy2eps(file="NeA_reg.eps", horizontal=F)
	print("NeA done.")

	# Plot line for migration rate of population 1 and save it in a .eps file
	abc.reg.mig1 <- makepd4(target,abc.rej[,15],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minmig1,maxmig1))
	plot(locfit(~abc.reg.mig1$x,weight=abc.reg.mig1$wt,xlim=c(minmig1,maxmig1)),main="prior (black) and posterior (blue) distributions",xlab="mig1",col="blue")
	plot(locfit(~priors[,15],xlim=c(minmig1,maxmig1)),col="black",add=T)
	dev.copy2eps(file="mig1_reg.eps", horizontal=F)
	print("mig1 done.")

	# Plot line for migration rate of population 2 and save it in a .eps file
	abc.reg.mig2 <- makepd4(target,abc.rej[,16],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minmig2,maxmig2))
	plot(locfit(~abc.reg.mig2$x,weight=abc.reg.mig2$wt,xlim=c(minmig2,maxmig2)),main="prior (black) and posterior (blue) distributions",xlab="mig2",col="blue")
	plot(locfit(~priors[,16],xlim=c(minmig2,maxmig2)),col="black",add=T)
	dev.copy2eps(file="mig2_reg.eps", horizontal=F)
	print("mig2 done.")

}
