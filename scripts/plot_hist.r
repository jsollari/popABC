# @author:    joao lopes
# @workplace: University of Reading
# @date:      13th May 2009

# Function to plot histograms of the posterior distributions of an ABC approach.
# Creates .eps image files in the directory from where it is run.

# @arg rej_file - file with rejection step results (.rej)
# @arg pri_file - file with sample of priors (.pri)
plot_hist <- function(rej_file , pri_file){

	#import the .rej file
	abc.rej <- data.matrix(read.table(rej_file))

	#import the .pri files
	priors <- data.matrix(read.table(pri_file))

	# Plot histograms for sequence DNA mutation rate and save it in a .eps file
	hist(abc.rej[,5],freq=F,main="prior (black) and posterior (blue) distributions",xlab="mut rate",col="blue",density=9)
	hist(priors[,5],freq=F,col="white",add=T)
	dev.copy2eps(file="mutrate.eps", horizontal=F)
	print("mut rate done.")

	# Plot histograms for sequence DNA recombination rate and save it in a .eps file
	hist(abc.rej[,9],freq=F,main="prior (black) and posterior (blue) distributions",xlab="rec rate",col="blue",density=9)
	hist(priors[,9],freq=F,col="white",add=T)
	dev.copy2eps(file="recrate.eps", horizontal=F)
	print("rec rate done.")

	# Plot histograms for splitting time and save it in a .eps file
	hist(abc.rej[,11],freq=F,main="prior (black) and posterior (blue) distributions",xlab="tev",col="blue",density=9)
	hist(priors[,11],freq=F,col="white",add=T)
	dev.copy2eps(file="tev.eps", horizontal=F)
	print("tev done.")

	# Plot histograms for effective size of population 1 and save it in a .eps file
	hist(abc.rej[,12],freq=F,main="prior (black) and posterior (blue) distributions",xlab="Ne1",col="blue",density=9)
	hist(priors[,12],freq=F,col="white",add=T)
	dev.copy2eps(file="Ne1.eps", horizontal=F)
	print("Ne1 done.")

	# Plot histograms for effective size of population 2 and save it in a .eps file
	hist(abc.rej[,13],freq=F,main="prior (black) and posterior (blue) distributions",xlab="Ne2",col="blue",density=9)
	hist(priors[,13],freq=F,col="white",add=T)
	dev.copy2eps(file="Ne2.eps", horizontal=F)
	print("Ne2 done.")

	# Plot histograms for effective size of ancestor population and save it in a .eps file
	hist(abc.rej[,14],freq=F,main="prior (black) and posterior (blue) distributions",xlab="NeA",col="blue",density=9)
	hist(priors[,14],freq=F,col="white",add=T)
	dev.copy2eps(file="NeA.eps", horizontal=F)
	print("NeA done.")

	# Plot histograms for migration rate of population 1 and save it in a .eps file
	hist(abc.rej[,15],freq=F,main="prior (black) and posterior (blue) distributions",xlab="mig1",col="blue",density=9)
	hist(priors[,15],freq=F,col="white",add=T)
	dev.copy2eps(file="mig1.eps", horizontal=F)
	print("mig1 done.")

	# Plot histograms for migration rate of population 2 and save it in a .eps file
	hist(abc.rej[,16],freq=F,main="prior (black) and posterior (blue) distributions",xlab="mig2",col="blue",density=9)
	hist(priors[,16],freq=F,col="white",add=T)
	dev.copy2eps(file="mig2.eps", horizontal=F)
	print("mig2 done.")

}
