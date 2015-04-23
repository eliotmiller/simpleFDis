#this function will take community-specific ordination results (i.e. a subset of total
#niche space, aka ordination.results) and a community road map, optionally a cov.matrix,
#and will return a new ordination space (i.e. a new community-level niche space) with the
#same number of species and observations per species as the input. First a multivariate
#normal distribution with species' centroids characteristics of comm.space is simulated.
#These are used as new species' centroids. Then the original species' points are used to
#simulate new multivariate normal distributions around these new centroids.

synthComm <- function(comm.space, comm.map, cov.matrix=NULL)
{
	if(is.data.frame(comm.map) != TRUE)
	{
		stop("comm.map must be a data frame")
	}

	#calculate each species' centroid, and the means of these centroids
	
	calcCenters <- centers(comm.space, comm.map)
	meanCenters <- apply(calcCenters, 2, mean)
	
	#if cov.matrix is missing, calculate it based on the species' centroids. otherwise
	#use the cov.matrix passed down (e.g. from a larger ordination space)

	if(is.null(cov.matrix))
	{
		cov.matrix <- cov(calcCenters)
	}
	
	#now use the centroids and the cov.matrix to generate a multivariate normal
	#distribution of the requisite number of species' centroids. will use these as new
	#species' centroids

	newCenters <- mvrnorm(n=dim(calcCenters)[1], mu=meanCenters, Sigma=cov.matrix)
	
	#split community ordination into single species and calculate covariance matrix for
	#each, saving each matrix into an element of a list

	covList <- list()
	
	for(i in 1:dim(comm.map)[1])
	{
		indices <- names(comm.map)[comm.map[i,]!=0]
		tempPoints <- comm.space[indices,]
		covList[[i]] <- cov(tempPoints)
	}
	
	#generate a multivariate normal distribution with these means and covariances for 
	#each species. VERY IMPORTANT. note that because we created a multivariate normal
	#distribution of new species' centroids, we are in effect randomizing the identity of
	#which new species gets which covariance structure, no need to randomize identities
	
	distributions <- list()
	
	for(i in 1:length(covList))
	{
		distributions[[i]] <- mvrnorm(n=sum(comm.map[i,] != 0), mu=newCenters[i,],
			Sigma=covList[[i]])
	}

	#add species names to the new community observations, since that is what is required
	#by nicheROVER downstream. first give each element of list a name, then rep it by the
	#length of that element. unlist into a character vector. then Reduce distributions
	#and cbind the species names in
	
	names(distributions) <- row.names(comm.map)
	species <- unlist(lapply(seq_along(distributions), function(x) 
		rep(x=names(distributions[x]), times=dim(distributions[[x]])[1])))
	
	distributions <- Reduce(rbind, distributions)
	distributions <- as.data.frame(distributions)
	distributions <- cbind(species, distributions)
	distributions
}
