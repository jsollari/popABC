# @author:    joao lopes
# @workplace: Instituto Gulbenkian de Ciencia
# @date:      12th May 2011

# Function to perform the regression step on a file with rejection step results.
# It also estimate the mode and the 0.95 credible interval of the obtained posterior distributions.
# Creates .txt file with results in the directory from where it is run

# @arg data_file - file with summary of 'real' data (.trg)
# @arg rej_file  - file with rejection step results (.rej)
get_modes <- function(data_file , rej_file){

	#import locfit library
	library(locfit)

	#import Mark Beaumont's scripts
	source("loc2plot_d.r")
	source("make_pd2.r")

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

    #credible interval to use
    credible <- 0.05
    
	#import the .rej file
	abc.rej <- data.matrix(read.table(rej_file))

	#import the .trg files
	target <- data.matrix(read.table(data_file))

    #create .txt file to store estimated results
    write("mode 0.025 0.975\n", "estimates.txt", append=F)

	# Plot line for sequence DNA mutation rate and save it in a .eps file
	abc.reg.mut <- makepd4(target,abc.rej[,5],abc.rej[,17:34],tol=1,rej=F,transf="none")
    write("mutation rate:", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.mut$x,wt=abc.reg.mut$wt, prob=credible, alpha=1.0, xlim=c(min(abc.reg.mut$x),max(abc.reg.mut$x))), "estimates.txt", append=T)
	print("mut rate done.")

	# Plot line for sequence DNA recombination rate and save it in a .eps file
	abc.reg.rec <- makepd4(target,abc.rej[,9],abc.rej[,17:34],tol=1,rej=F,transf="none")
    write("recombination rate:", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.rec$x,wt=abc.reg.rec$wt, prob=credible, alpha=1.0, xlim=c(min(abc.reg.rec$x),max(abc.reg.rec$x))), "estimates.txt", append=T)
	print("rec rate done.")

	# Plot line for splitting time and save it in a .eps file
	abc.reg.tev <- makepd4(target,abc.rej[,11],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(mintev,maxtev))
    write("splitting time (tev):", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.tev$x,wt=abc.reg.tev$wt, prob=credible, alpha=1.0, xlim=c(mintev,maxtev)), "estimates.txt", append=T)
	print("tev done.")

	# Plot line for effective size of population 1 and save it in a .eps file
	abc.reg.ne1 <- makepd4(target,abc.rej[,12],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minNe1,maxNe1))
    write("effective population size (Ne1):", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.ne1$x,wt=abc.reg.ne1$wt, prob=credible, alpha=1.0, xlim=c(minNe1,maxNe1)), "estimates.txt", append=T)
	print("Ne1 done.")

	# Plot line for effective size of population 2 and save it in a .eps file
	abc.reg.ne2 <- makepd4(target,abc.rej[,13],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minNe2,maxNe2))
    write("effective population size (Ne2):", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.ne2$x,wt=abc.reg.ne2$wt, prob=credible, alpha=1.0, xlim=c(minNe2,maxNe2)), "estimates.txt", append=T)
	print("Ne2 done.")

	# Plot line for effective size of ancestor population and save it in a .eps file
	abc.reg.nea <- makepd4(target,abc.rej[,14],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minNeA,maxNeA))
    write("effective population size (NeA):", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.nea$x,wt=abc.reg.nea$wt, prob=credible, alpha=1.0, xlim=c(minNeA,maxNeA)), "estimates.txt", append=T)
	print("NeA done.")

	# Plot line for migration rate of population 1 and save it in a .eps file
	abc.reg.mig1 <- makepd4(target,abc.rej[,15],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minmig1,maxmig1))
    write("migration rate (mig1):", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.mig1$x,wt=abc.reg.mig1$wt, prob=credible, alpha=1.0, xlim=c(minmig1,maxmig1)), "estimates.txt", append=T)
	print("mig1 done.")

	# Plot line for migration rate of population 2 and save it in a .eps file
	abc.reg.mig2 <- makepd4(target,abc.rej[,16],abc.rej[,17:34],tol=1,rej=F,transf="none",bb=c(minmig2,maxmig2))
    write("migration rate (mig2):", "estimates.txt", append=T)
	write(loc1statsx(x=abc.reg.mig2$x,wt=abc.reg.mig2$wt, prob=credible, alpha=1.0, xlim=c(minmig2,maxmig2)), "estimates.txt", append=T)
	print("mig2 done.")

}
