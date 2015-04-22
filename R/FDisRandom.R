#' @export

#this function explores how FDis changes with increasing observations. no road map is
#needed, as it is generated within the function. samples is a vector over which the
#exploration is desired

FDisRandom <- function(ordination.results, samples, iterations, cores=FALSE)
{	
	if(cores==FALSE)
	{
		#set up a blank results matrix to save into.
		results <- matrix(nrow=iterations, ncol=length(samples))
	
		for(i in 1:iterations)
		{
			#set up a blank road map to fill with 1s according to desired parameters
			roadMap <- matrix(ncol=dim(ordination.results)[1], nrow=length(samples), 0)
			row.names(roadMap) <- samples
			colnames(roadMap) <- 1:dim(ordination.results)[1]
			
			for(j in 1:length(samples))
			{
				#sample some of the 0s and convert to 1s according to desired number of sample
				roadMap[j,][sample(colnames(roadMap),samples[j])] <- 1
			}
		
			temp <- FDis(ordination.results, roadMap)
		
			results[i,] <- temp
		}
	}

	if(cores!=FALSE)
	{
		if(!is.numeric(cores))
		{
			stop("cores must be set either to FALSE or to a numeric value")
		}
		require(doMC)
		require(foreach)
		registerDoMC(cores)

		results <-
		foreach(i = 1:iterations, .combine='rbind') %dopar%
		{
			#set up a blank road map to fill with 1s according to desired parameters
			roadMap <- matrix(ncol=dim(ordination.results)[1], nrow=length(samples), 0)
			row.names(roadMap) <- samples
			colnames(roadMap) <- 1:dim(ordination.results)[1]
			
			foreach(j = 1:length(samples), .combine='c') %do%
			{
				#sample some of the 0s and convert to 1s according to desired number of sample
				roadMap[j,][sample(colnames(roadMap),samples[j])] <- 1
			}
		
			results <- FDis(ordination.results, roadMap)
		}
	}

	results <- as.data.frame(results)
	names(results) <- samples
	
	results
}
