# @author:    joao lopes
# @workplace: University of Reading
# @date:      20th April 2010

# Function to perform model-choice on a file with rejection step results.
# It also plots a bar plor of the obtained posterior distribution.
# Creates .eps image files in the directory from where it is run

# @arg data_file - file with summary of 'real' data (.trg)
# @arg rej_file  - file with rejection step results (.rej)
# @arg pri_file  - file with sample of priors (.pri)
model_choice <- function(data_file , rej_file ){

	#import VGAM library
	library(VGAM)

	#import Mark Beaumont's scripts
	source("calmod.r")

	#import the .rej file
	abc.rej <- data.matrix(read.table(rej_file))

	#import the .trg files
	target <- data.matrix(read.table(data_file))

	#Bar plot for model-choice (rejection)
	top1<-length(which(abc.rej[,2]==1))/length(abc.rej[,2])
	top2<-length(which(abc.rej[,2]==2))/length(abc.rej[,2])
	print(c(top1,top2))
	write(c(top1,top2),"modelprob_rej.txt")
	barplot(c(top1,top2),names=c("Model1","Model2"))
	abline(h=1.0/2,col="red")

	#Bar plot for model-choice (regression)
	res.topol<-calmod(target,abc.rej[,2],abc.rej[,17:34],1,rej=F)
	print(res.topol$x2)
	write(res.topol$x2,"modelprob_reg.txt")
	barplot(res.topol$x2,names=c("Model1","Model2"))
	abline(h=1.0/2,col="red")
	print("model-choice done.")

}
