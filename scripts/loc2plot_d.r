# This plots countours for bivariate density 
#x = random variable 1
#y = random variable 2
#cprob = vector of hpd values for contours; depending on how many values in cprob 
#        determines the number of contour lines
#alpha = the smoothing parameter for locfit; small numbers, more jaggedy; large numbers
#          more smooth.
#xlim = the support limits of the density given as (xmin,ymin,xmax,ymax)
#gxlim = the limits of the density that are actually plotted.
#maxk = some number, the bigger it is the less likely locfit is to complain 
#       (but makes it slower)

loc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100)
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}


#like above, but allows for the points to be weighted - i.e. for output of ABC regression
loc2plotw <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,wt,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100,weight=wt)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100,weight=wt)
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

#gives the HPD value for an observation px,py in a density constructed from x, y.
gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100)
#	d1 <- (x-px)^2+(y-py)^2
#	best <- d1 == min(d1)
#	lev <- mean(fitted(fit)[best])
	lev <- predict(fit,list(px,py))
	slev <- sort(fitted(fit))
	indic <- slev <= lev
	sum(indic)/length(x)
}
	
	
# finds the mode for a bivariate density

loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	maxk=maxk,mint=100,cut=0.8,maxit=100)
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
	xlim=xlim,maxk=maxk,mint=100,cut=0.8,maxit=100)
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}


# this works for univariate data; gives a vector with 
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data. 
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,prob,alpha=0.5,xlim,...) 
{
	if(missing(xlim)){
		fit <- locfit(~x,alpha=alpha)
	}
	else {
		fit <- locfit(~x,alpha=alpha,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}
		
	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}

#this does the hpd calculation in a different way. Should give similar answers
#to loc1stats. It is best to put xlim in directly. wt is a vector of weights and is optional.
#numpoint is the number of points to do interpolation - the more the better
loc1statsx <- function(x,prob,alpha=0.5,xlim,wt,numpoint=10000,...) 
{
	if(missing(xlim)){
		if(min(x) < 0)x.min <- 1.1*min(x)
		else x.min <- min(x)*0.9
		if(max(x) < 0)x.max <- 0.9*max(x)
		else x.max <- 1.1*max(x)
		print(paste("putting in these xlimits from the data:",x.min,x.max))
		xlim <- c(x.min,x.max)
	}

	if(missing(wt))fit <- locfit(~x,alpha=alpha,xlim=xlim)
	else fit <- locfit(~x,alpha=alpha,xlim=xlim,weight=wt)
	xx <- seq(xlim[1],xlim[2],len=numpoint)
	yy <- predict(fit,xx)
	sum1 <- sum(yy)
	x.modef <- max(yy)
	x.mode <- xx[yy == x.modef]
	if(length(x.mode)>1)x.mode <- x.mode[1]
	
	yy2 <- sort(yy)
	pval <- 0
	for(j in 1:numpoint){
		pval <- pval+yy2[j]/sum1
		if(pval > prob)break
	}
	lev <- yy2[j]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	ii <- 2
	flip <- TRUE
	for(j in 2:length(xx)){
		if(flip && yy[j] > lev){
			l1[[ii]] <- xx[j-1]
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && yy[j] < lev){
			l1[[ii]] <- xx[j]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip &&  j == length(xx)){
			l1[[ii]] <- xx[j]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)

}

#gives the HPD value for an observation px,in a univariate density constructed from x.
#uses same method as loc1statsx, so need to specify xlim. wt is optional weight.
#numpoint is the number of points to do interpolation - the more the better
gethpdprob1 <- function(x,px,alpha=0.5,xlim,wt,numpoint=10000,...)
{

	if(missing(wt))fit <- locfit(~x,alpha=alpha,xlim=xlim)
	else fit <- locfit(~x,alpha=alpha,xlim=xlim,weight=wt)
	xx <- seq(xlim[1],xlim[2],len=numpoint)
	yy <- predict(fit,xx)
	lev <- predict(fit,px)
	sum1 <- sum(yy)

	indic <- yy <= lev
	sum(yy[indic])/sum1
}


#given a function described by x,y points, this works like loc1stats.
#uses spline interpolation between points.
#numpoint is the number of points to do interpolation - the more the better
linestats <- function(x,y,prob,numpoint=10000)
{
	n <- length(x)
	minval <- x[1]
	maxval <- x[n]
	# This bit just guarantees that irrespective of what is 
	# in x and y, we have numpoint *evenly* spaced (interpolated) points between
	# max and min
	res1 <- spline(x,y,numpoint,xmin=minval,xmax=maxval)
	xx <- res1$x
	yy <- res1$y
	sum1 <- sum(yy)
	x.modef <- max(yy)
	x.mode <- xx[yy == x.modef]
	if(length(x.mode)>1)x.mode <- x.mode[1]
	
	yy2 <- sort(yy)
	pval <- 0
	for(j in 1:numpoint){
		pval <- pval+yy2[j]/sum1
		if(pval > prob)break
	}
	lev <- yy2[j]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	ii <- 2
	flip <- TRUE
	for(j in 2:length(xx)){
		if(flip && yy[j] > lev){
			l1[[ii]] <- xx[j-1]
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && yy[j] < lev){
			l1[[ii]] <- xx[j]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip &&  j == length(xx)){
			l1[[ii]] <- xx[j]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}

