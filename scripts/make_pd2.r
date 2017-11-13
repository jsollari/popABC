makepd4 <- function(target,x,sumstat,tol,gwt,rejmethod=T,transf="none",bb=c(0,0))
{
# target is the set of target summary stats
# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
# sumstat is an array of simulated summary stats (i.e. independent variables).
# NBB this function originally used lm() and assumed 4 summary stats, and I edited by hand for other numbers. 
# NBB I've now modified it using lsfit() (following Shola Ajayi) so that it will take an arbitrary number of summary stats.
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x regression adjusted values
# $vals - unadjusted values in rejection region (i.e. normal rejection)
# $wt - the regression weight (i.e. the Epanechnikov weight)
# $ss - the sumstats corresponding to these points
# $predmean - estimate of the posterior mean
# $fv - the fitted value from the regression

if(sum(transf == c("none","log","logit")) == 0){
	stop("transf must be none, log, or logit")
	
}
if(transf=="logit"){
	if(bb[1] >= bb[2]){
		stop("bounds wrong for logit")
		
	}
}

if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))

nss <- length(sumstat[1,])


# scale everything 

    scaled.sumstat <- sumstat
    
    for(j in 1:nss){
    
    	scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    }
    target.s <- target
    
    for(j in 1:nss){
    
    	target.s[j] <- normalise(target[j],sumstat[,j][gwt])
    }
    
# calc euclidean distance

    sum1 <- 0
    for(j in 1:nss){
    	sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
   }
   dst <- sqrt(sum1)
# includes the effect of gwt in the tolerance
    dst[!gwt] <- floor(max(dst[gwt])+10)
    

# wt1 defines the region we're interested in 
    abstol <- quantile(dst,tol)
    wt1 <- dst < abstol
    
    if(transf == "log"){
    	if(min(x) <= 0){
    		print("log transform: val out of bounds - correcting")
    		x.tmp <- ifelse(x <= 0,max(x),x)
    		x.tmp.min <- min(x.tmp)
    		x <- ifelse(x <= 0, x.tmp.min,x)
    	}
    	x <- log(x)
    }
    else if(transf == "logit"){
    	if(min(x) <= bb[1]){
    		x.tmp <- ifelse(x <= bb[1],max(x),x)
    		x.tmp.min <- min(x.tmp)
    		x <- ifelse(x <= bb[1], x.tmp.min,x)
    		
    	}
    	if(max(x) >= bb[2]){
    		x.tmp <- ifelse(x >= bb[2],min(x),x)
    		x.tmp.max <- max(x.tmp)
    		x <- ifelse(x >= bb[2], x.tmp.max,x)
    		
    	}
    	x <- (x-bb[1])/(bb[2]-bb[1])
    	x <- log(x/(1-x))
    }
    
    if(rejmethod){
        l1 <- list(x=x[wt1],wt=0)
    }
    else{
        regwt <- 1-dst[wt1]^2/abstol^2
        
        fit1 <- lsfit(scaled.sumstat[wt1,],x[wt1],wt=regwt)
        
        predmean <- fit1$coeff %*% c(1,target.s)
        
        l1 <- list(x=fit1$residuals+predmean,vals=x[wt1],wt=regwt,ss=sumstat[wt1,],predmean=predmean,fv = x[wt1]-fit1$residuals,transf=transf)
        
    }
    if(transf == "log"){
    	l1$x <- exp(l1$x)
    	l1$vals <- exp(l1$vals)
    }
    else if(transf == "logit"){
    	l1$x <- exp(l1$x)/(1+exp(l1$x))
    	l1$x <- l1$x*(bb[2]-bb[1])+bb[1]
    	l1$vals <- exp(l1$vals)/(1+exp(l1$vals))
    	l1$vals <- l1$vals*(bb[2]-bb[1])+bb[1]
    }
    
    l1
}


normalise <- function(x,y){

if(var(y) == 0)return (x - mean(y))
(x-(mean(y)))/sqrt(var(y))
}


