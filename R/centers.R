#this function returns the weighted centroids per cloud of points in multivariate space
#given some ordination results and a road map telling which point belongs to which cluster
#e.g. a file you have occasionally called "secondTable.csv"

centers <- function(ordination.results, road.map)
{
	results <- matrix(ncol=dim(ordination.results)[2], nrow=dim(road.map)[1])
	for(i in 1:dim(road.map)[1])
	{
		results[i,] <- apply(ordination.results, 2, weighted.mean, w=road.map[i,])
	}
	
	row.names(results) <- row.names(road.map)
	
	return(results)
}
