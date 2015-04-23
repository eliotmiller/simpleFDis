#this function will take an input cdm, which is just a single row from the overall CDM,
#a subset of the total ordination space specific to that quadrat, and a community road 
#map, explaining which point belongs to which species. it will then use the
#results of synthComm to calculate pairwise niche overlap according to various parameters

synthPairwise <- function(comm.cdm, comm.space, comm.map, cov.matrix=NULL, nsamples=1000, 
	nprob=1000, alpha=0.95, abundance.weighted=FALSE)
{
	#simulate a new community
	newComm <- synthComm(comm.space, comm.map, cov.matrix)

	#generate relevant parameters from draws of default posteriors
	parameters <- tapply(1:nrow(newComm), newComm$species, 
		function(ii) niw.post(nsamples = nsamples, X = newComm[ii, 2:dim(newComm)[2]]))

	#calculate pairwise overlaps, i.e. A with B and B with A
	allOverlaps <- overlap(parameters, nreps=nsamples, nprob=nprob, alpha=alpha)

	#calculate the mean overlap across iterations and return in an array
	meanOverlaps <- apply(allOverlaps, c(1:2), mean) * 100

	#this gave you a matrix where the lower triangle is the overlap of species A within B,
	#and the upper is the overlap of species B within A. if you transpose the matrix then
	#combine them into two slices of an array, you can take the means and it will be the 
	#average overlap between each pairwise comparison

	transposed <- t(meanOverlaps)

	theArray <- array(c(meanOverlaps, transposed), 
		dim=c(dim(meanOverlaps)[1], dim(meanOverlaps)[1], 2))

	realOverlaps <- apply(theArray, c(1:2), mean, na.rm=TRUE)

	#give the realOverlaps species names so that you can use the comm.cdm to calculate
	#abundance-weighted mean pairwise niche overlaps. IMPORTANT! note that if you set the
	#diagonal element of the overlaps equal to 100, which is biologically correct, if you
	#use either intraspecific or complete MPD then your niche overlaps are going to be 
	#correlated with species richness. however, interspecific sets the weight of the
	#diagonal element of the overlap matrix equal to 0 for the mean calculation. because
	#the diagonal element is not used in the mean if abundance.weighted = FALSE, then 
	#setting the diagonal equal to 100 will not affect this either
	
	diag(realOverlaps) <- 100
	
	colnames(realOverlaps) <- names(comm.cdm)
	row.names(realOverlaps) <- names(comm.cdm)
	
	results <- modifiedMPD(comm.cdm, realOverlaps, abundance.weighted=abundance.weighted)
	
	results
}
