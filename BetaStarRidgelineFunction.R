################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script defines the function "tbs.ridgeline" that calculates beta-star
# based on the ridgeline manifold for two multivariate normal distributions
# (equation 4 in in page 2045 of Ray and Lindsay. 2005. The topography of
# multivariate normal mixtures. The Annals of Statistics 33: 2042-2065).
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# This script may run on its own, without any data files.
#
################################################################################
################################################################################

#The function "tbs.ridgeline" estimates true.beta.star between two species, say
# species A and B, with the following arguments:
#A.M, the vector of means of species A
#B.M, the vector of means of species B
#A.VCV, the variance-covariance matrix of species A 
#B.VCV, the variance-covariance matrix of species B

tbs.ridgeline <- function(A.M, B.M, A.VCV, B.VCV)
{
	#Check validity of argument values.
	if(!identical(length(A.M), dim(A.VCV)[1], dim(A.VCV)[2])) stop("A.M and A.VCV are non-conformable arguments")
	if(!identical(length(B.M), dim(B.VCV)[1], dim(B.VCV)[2])) stop("B.M and B.VCV are non-conformable arguments")
	if(!identical(length(A.M), length(B.M))) stop("A.M and B.M have different lengths")

	#Calculate ridgeline manifold.
	#The code for the ridgeline manifold is based on equation 4 in page 2045 of
	#Ray and Lindsay, 2005. "The topography of multivariate normal mixtures" The
	#Annals of Statistics 33: 2042-2065. the ridgeline manifold is evaluated at
	#various values of alpha and the resulting coordinates (in phenotypic space)
	#are captured in matrix ridgeline.coor.
	alpha <- seq(0,1,0.001)
	p.dimensions <- length(A.M)
	ridgeline.coor <- matrix(NA, nrow=length(alpha), ncol=p.dimensions)
	for (i in 1:length(alpha))
	{
		a <- solve((1-alpha[i])*solve(A.VCV)  +  alpha[i]*solve(B.VCV))
		b <- (1-alpha[i])*solve(A.VCV)%*%A.M  +  alpha[i]*solve(B.VCV)%*%B.M
		d <- a %*% b
		ridgeline.coor[i,] <- d
	}

	#Calculate true.beta.star
	md.A <- mahalanobis(ridgeline.coor, A.M, A.VCV)
	md.B <- mahalanobis(ridgeline.coor, B.M, B.VCV)
	true.beta.vector.A <- pchisq(md.A, p, lower.tail = TRUE, log.p = FALSE)
	true.beta.vector.B <- pchisq(md.B, p, lower.tail = TRUE, log.p = FALSE)
	beta.dif <- abs(true.beta.vector.A - true.beta.vector.B)
	o.beta.dif <- order(beta.dif)
	true.beta.star <- mean(c(true.beta.vector.A[o.beta.dif[1]], true.beta.vector.B[o.beta.dif[1]]))

	#Return a list with the estimate of true.beta.star and a matrix with alpha
	#values, the respective phenotypic coordinates of the ridgeline manifold,
	#and the respective proportion of the distribution of each species contained
	#within (hyper)elliptical regions.
	ridgeline <- cbind(alpha, ridgeline.coor, true.beta.vector.A, true.beta.vector.B)
	colnames(ridgeline) <- c("alpha", paste("ridgeline.coor", 1:p.dimensions, sep="."), "true.beta.vector.A", "true.beta.vector.B")	
	return(list(true.beta.star = true.beta.star, ridgeline = ridgeline))
}
