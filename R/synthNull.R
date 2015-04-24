synthNull <- function(cdm, total.space, total.map, cov.matrix="comm", nsamples=1000, 
	nprob=1000, alpha=0.95, abundance.weighted=FALSE)
{
	#this has been tested with a data frame, not a matrix, so ensure the road map is of
	#class data frame
	if(is.data.frame(total.map) != TRUE)
	{
		stop("total.map must be a data frame")
	}

	#if cov.matrix is set to overall, use the covariance matrix of all species points in
	#multivariate space, else set to null and allow subsequent functions to calculate it
	#per community based on the species present in that community
	if(cov.matrix=="overall")
	{
		calcCenters <- centers(total.space, total.map)
		cov.matrix <- cov(calcCenters)
	}
	
	else if(cov.matrix=="comm")
	{
		cov.matrix <- NULL
	}

	else
	{
		stop("cov.matrix must be set to either overall or comm")
	}
	
	#set up a blank vector to save each quadrat results into
	results <- c()
	
	#for each row in the CDM
	for(i in 1:dim(cdm)[1])
	{
		#cut the cdm just down to that row
		comm.cdm <- cdm[i,]
		#exclude species that do not occur in that row
		comm.cdm <- comm.cdm[,comm.cdm != 0]
		#make code easier by just specifically identifying the species that occur
		species <- names(comm.cdm)
		#cut the total road map down to only species that occur in that row of the cdm
		comm.map <- total.map[row.names(total.map) %in% species,]
		#now cut out points that do not belong to any of the species that are included
		comm.map <- comm.map[,apply(comm.map, 2, sum) != 0]
		#cut the total ordination space down to only those points that belong to species
		#that occur in that row of the cdm
		comm.space <- total.space[row.names(total.space) %in% names(comm.map),]
		#now use these pared down data frames to calculate the pairwise niche overlap
		#of a simulated community with the same characteristics as the input row
		results[i] <- synthPairwise(comm.cdm=comm.cdm, comm.space=comm.space,
			comm.map=comm.map, cov.matrix=cov.matrix, nsamples=nsamples, nprob=nprob,
			alpha=alpha, abundance.weighted=abundance.weighted)
	}

	#convert results into a data frame, where quadrat is the first column
	output <- data.frame(quadrat=row.names(cdm), overlap=results)
	
	output
}
