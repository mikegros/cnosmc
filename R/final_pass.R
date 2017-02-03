
final_pass <- function(smc_samples, n_samples,p_link,nreps = 5){
  # Feb 3 2017, 
  # use this as a final pass to correct for inconsistencies in the posterior
  # caused by mystery memory leakage issues.
  #
  #Need to either keep old_post as an item in smc_samples or run the
  # new_post line twice to define it with independent calculations as
  # new and old post (as it is done here).
  # Need to 
  
  
  paramsList          <- defaultParametersFuzzy(data, model)
  indexList           <- indexFinder(CNOlist = data, model = model)
  
  # get a kernel density estimate of the last stage output
  # consider [g,n,k,Gstring] blocks to be independent
  # consider sigma squared to also be independent
  #
  old_logpost_from_kde = 0
    for(lp in 1:dim(smc_samples$nCube)[2]){
      x = cbind(smc_samples$gCube[,lp],smc_samples$nCube[,lp],
                smc_samples$kCube[,lp],smc_samples$Gstring[,lp])
      out = kde(x,eval.points=x)
      old_logpost_from_kde = old_logpost_from_kde + log(out$estimate)
    }
  for(lp in 1:dim(smc_samples$sigsq)[2]){
    x = smc_samples$sigsq[,lp]
    out = kde(x,eval.points=x)
    old_logpost_from_kde = old_logpost_from_kde + log(out$estimate)
  }
  
    
  use_post = NULL
  # Uncomment if building with excess cluster.  Leave commented to go sequential with nreps repeats
  # Note that excess_cluster_call is not currently an input
  # Also cl1 is not defined in here
  # if( n_cores > 1 & !excess_cluster_call ){
  #   for(lp in 1:nreps){
  #     new_post <- (parSapply(cl,1:n_samples, function(samp){
  #       posterior(cl         = cl1,
  #             Bstring    = rep(1,length(smc_samples$Gstring[samp,])),
  #             Gstring    = smc_samples$Gstring[samp,],
  #             gCube      = smc_samples$gCube[samp,],
  #             nCube      = smc_samples$nCube[samp,],
  #             kCube      = smc_samples$kCube[samp,],
  #             sigsq      = smc_samples$sigsq[samp,],
  #             p_link     = p_link,
  #             model      = model,
  #             paramsList = paramsList,
  #             indexList  = indexList)
  #   }))
  #     use_post = apply(cbind(new_post,use_post),1,max)  
  #   }
  # } else {
  for(lp in 1:nreps){
    new_post <- (sapply(1:n_samples, function(samp){
        posterior(cl         = NULL,
              Bstring    = rep(1,length(smc_samples$Gstring[1,])),
              Gstring    = smc_samples$Gstring[samp,],
              gCube      = smc_samples$gCube[samp,],
              nCube      = smc_samples$nCube[samp,],
              kCube      = smc_samples$kCube[samp,],
              sigsq      = smc_samples$sigsq[samp,],
              p_link     = p_link,
              model      = model,
              paramsList = paramsList,
              indexList  = indexList)
      }))
      use_post = apply(cbind(new_post,use_post),1,max)
   }
# }

w  <- log(smc_samples$w) + use_post - old_logpost_from_kde # since w is reset to 1/n for all values after jittering
w  <- w - max(w)
w  <- exp(w)
w  <- w/sum(w)

# Resample parameters using weights

  resample_inds <- sample(1:n_samples,n_samples,replace = TRUE,prob=w)
  
  smc_samples$gCube   <- smc_samples$gCube[resample_inds,]
  smc_samples$nCube   <- smc_samples$nCube[resample_inds,]
  smc_samples$kCube   <- smc_samples$kCube[resample_inds,]
  smc_samples$Gstring <- smc_samples$Gstring[resample_inds,]
  smc_samples$sigsq   <- smc_samples$sigsq[resample_inds,]
  smc_samples$w       <- rep(1,n_samples)/n_samples
  smc_samples$use_post             =             use_post[resample_inds]
  smc_samples$old_logpost_from_kde = old_logpost_from_kde[resample_inds]
  w <- rep(1,n_samples)/n_samples

  
}
