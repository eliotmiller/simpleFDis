#this function determines the weighted centroid of each cloud of points and then
#determines the mean absolute deviation (weighted) from each centroid per cloud

FDis <- function(ordination.results, road.map)
{	
	require(ecodist)

	results <- c()

	centerPoints <- centers(ordination.results, road.map)

	temp <- rbind(ordination.results, centerPoints)
	allDistances <- dist(temp)
	distMatrix <- full(allDistances)
	
	for(i in 1:dim(road.map)[1])
	{
		#figure out which column corresponds to distances from the weighted centroids
		#pull it out and get rid of the final elements (distances among centroids)
		centroidDists <- distMatrix[,dim(road.map)[2]+i]
		centroidDists <- centroidDists[1:(length(centroidDists)-(dim(road.map)[1]))]
		#now calculate the weighted mean distance from this centroid
		results[i] <- weighted.mean(centroidDists, road.map[i,])
	}
	
	names(results) <- row.names(road.map)

	return(results)
}
