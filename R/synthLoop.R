synthLoop <- function(cdm, total.space, total.map, cov.matrix="comm", nsamples=1000, 
	nprob=1000, alpha=0.95, abundance.weighted=FALSE, iterations, cores)
{
	registerDoMC(cores)
	
	results <- foreach(i=1:iterations, .combine='rbind') %dopar%
	{
		synthNull(cdm=cdm, total.space=total.space, total.map=total.map,
			cov.matrix=cov.matrix, nsamples=nsamples, nprob=nprob, alpha=alpha,
			abundance.weighted=abundance.weighted)
	}
	
	results
}
