# @author:    joao lopes
# @workplace: University of Reading
# @date:      13th May 2009

# Function to plot density lines of the posterior distributions
# of an ABC approach using the package 'locfit'
# Creates .eps image files in the directory from where it is run.

# @arg rej_file - file with rejection step results (.rej)
# @arg pri_file - file with sample of priors (.pri)
plot_line <- function(rej_file , pri_file){

	#import locfit library
	library(locfit)

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

	# Plot line for sequence DNA mutation rate and save it in a .eps file
	plot(locfit(~abc.rej[,5]),main="prior (black) and posterior (blue) distributions",xlab="mut rate",col="blue")
	plot(locfit(~priors[,5]),col="black",add=T)
	dev.copy2eps(file="mutrate_line.eps", horizontal=F)
	print("mut rate done.")

	# Plot line for sequence DNA recombination rate and save it in a .eps file
	plot(locfit(~abc.rej[,9]),main="prior (black) and posterior (blue) distributions",xlab="rec rate",col="blue")
	plot(locfit(~priors[,9]),col="black",add=T)
	dev.copy2eps(file="recrate_line.eps", horizontal=F)
	print("rec rate done.")

	# Plot line for splitting time and save it in a .eps file
	plot(locfit(~abc.rej[,11],xlim=c(mintev,maxtev)),main="prior (black) and posterior (blue) distributions",xlab="tev",col="blue")
	plot(locfit(~priors[,11],xlim=c(mintev,maxtev)),col="black",add=T)
	dev.copy2eps(file="tev_line.eps", horizontal=F)
	print("tev done.")

	# Plot line for effective size of population 1 and save it in a .eps file
	plot(locfit(~abc.rej[,12],xlim=c(minNe1,maxNe1)),main="prior (black) and posterior (blue) distributions",xlab="Ne1",col="blue")
	plot(locfit(~priors[,12],xlim=c(minNe1,maxNe1)),col="black",add=T)
	dev.copy2eps(file="Ne1_line.eps", horizontal=F)
	print("Ne1 done.")

	# Plot line for effective size of population 2 and save it in a .eps file
	plot(locfit(~abc.rej[,13],xlim=c(minNe2,maxNe2)),main="prior (black) and posterior (blue) distributions",xlab="Ne2",col="blue")
	plot(locfit(~priors[,13],xlim=c(minNe2,maxNe2)),col="black",add=T)
	dev.copy2eps(file="Ne2_line.eps", horizontal=F)
	print("Ne2 done.")

	# Plot line for effective size of ancestor population and save it in a .eps file
	plot(locfit(~abc.rej[,14],xlim=c(minNeA,maxNeA)),main="prior (black) and posterior (blue) distributions",xlab="NeA",col="blue")
	plot(locfit(~priors[,14],xlim=c(minNeA,maxNeA)),col="black",add=T)
	dev.copy2eps(file="NeA_line.eps", horizontal=F)
	print("NeA done.")

	# Plot line for migration rate of population 1 and save it in a .eps file
	plot(locfit(~abc.rej[,15],xlim=c(minmig1,maxmig1)),main="prior (black) and posterior (blue) distributions",xlab="mig1",col="blue")
	plot(locfit(~priors[,15],xlim=c(minmig1,maxmig1)),col="black",add=T)
	dev.copy2eps(file="mig1_line.eps", horizontal=F)
	print("mig1 done.")

	# Plot line for migration rate of population 2 and save it in a .eps file
	plot(locfit(~abc.rej[,16],xlim=c(minmig2,maxmig2)),main="prior (black) and posterior (blue) distributions",xlab="mig2",col="blue")
	plot(locfit(~priors[,16],xlim=c(minmig2,maxmig2)),col="black",add=T)
	dev.copy2eps(file="mig2_line.eps", horizontal=F)
	print("mig2 done.")

}
