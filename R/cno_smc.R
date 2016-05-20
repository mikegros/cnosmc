cno_smc <- function(n_samples, data, model,
                    init_links  = 1,
                    p_link      = 0.9,
                    n_mh        = 5,
                    jump_size   = rep(0.15,3),
                    sigma       = 0.1,
                    split_inhib = FALSE,
                    n_cores     = 1,
                    diagnostics = F){
  #
  # Last update: April 26, 2016, Mike Grosskopf
  #
  # [Apr 26,2016] Mike Grosskopf
  #               Took SMC code from individual scripts and made it a standalone file.
  #                 Did tons of code cleaning to get it completely working and clean as well
  #

  inhib_inds <- NULL
  if (split_inhib){
    inhib_settings <- unique(data$valueInhibitors)

    if (ncol(data$valueInhibitors) == 1){
      inhib_inds <- lapply(inhib_settings, function(x) which(data$valueInhibitors == x))
    } else {
      inhib_inds <- apply(inhib_settings,1,function(y) which(apply(data$valueInhibitors,1,function(x) all(x==y))))
    }
  }

  paramsList          <- defaultParametersFuzzy(data, model)
  indexList           <- indexFinder(CNOlist = data, model = model)

  if( n_cores > 1 ){
    # outfile arguement tells R where to write print console output (print statements, etc) when running in parallel
    # cl <- makeCluster(2,outfile="~/output_R.txt")
    cl <- makeCluster(4)
    clusterCall(cl,function(x) {library(CNORfuzzy)})
  }else{
    cl <- NULL
  }


  findMainEffects <- unlist(lapply(strsplit(model$reacID,'\\+'),function(x) {length(x)<2} ))
  n_params        <- sum(findMainEffects)

  test_bString    <- rep(0,n_params)

  init_gCube <- 1*runif(n_params*n_samples)
  init_nCube <- rexp(n_params*n_samples,1/2)
  init_kCube <- 1*runif(n_params*n_samples)

  # Turn on some links for initial subgraph
  test_bString[c(init_links)] <- 1

  if (split_inhib){
    n_models     <- length(inhib_inds)
    init_Gstring <- rbinom(n_models*n_params*n_samples,1,p_link)
  }else{
    n_models     <- 1
    init_Gstring <- rbinom(n_params*n_samples,1,p_link)
  }

  smc_samples <- list(gCube = matrix(init_gCube, ncol = n_params, nrow = n_samples, byrow = TRUE))

  smc_samples$nCube   <- matrix(init_nCube,   ncol = n_params,          nrow = n_samples, byrow = TRUE)
  smc_samples$kCube   <- matrix(init_kCube,   ncol = n_params,          nrow = n_samples, byrow = TRUE)
  smc_samples$Gstring <- matrix(init_Gstring, ncol = n_models*n_params, nrow = n_samples, byrow = TRUE)

  if (split_inhib){
    colnames(smc_samples$Gstring) <- rep(colnames(model$interMat)[1:n_params],n_models)
  }else{
    colnames(smc_samples$Gstring) <- colnames(model$interMat)[1:n_params]
  }

  old_post <- sapply(1:n_samples, function(samp){
    LogpriorGstring(smc_samples$Gstring[samp,],p_link) +
      Logpriorg(smc_samples$gCube[samp,]) +
      Logpriorn(smc_samples$nCube[samp,]) +
      Logpriork(smc_samples$kCube[samp,]) })

  for (stage in 1:n_params){

    print(paste("Stage: ",stage,sep=""))

    # Get likelihood weights
    if(n_cores>1) clusterExport(cl,varlist=ls(),envir = environment())

    new_post <- (sapply(1:n_samples, function(samp){
      posterior(cl,test_bString,
                smc_samples$Gstring[samp,],
                smc_samples$gCube[samp,],
                smc_samples$nCube[samp,],
                smc_samples$kCube[samp,],
                p_link,inhib_inds,model,
                paramsList,indexList,sigma)
    }))

    w  <- new_post - old_post
    w  <- w-max(w)
    w  <- exp(w)

    # Resample parameters using likelihood weights
    resample_inds <- sample(1:n_samples,n_samples,replace = TRUE,prob=w)

    smc_samples$gCube   <- smc_samples$gCube[resample_inds,]
    smc_samples$nCube   <- smc_samples$nCube[resample_inds,]
    smc_samples$kCube   <- smc_samples$kCube[resample_inds,]
    smc_samples$Gstring <- smc_samples$Gstring[resample_inds,]

    # Perturb resampled values with an MH step

    if(n_cores>1) clusterExport(cl,varlist=ls(),envir = environment())
    tmp <- sapply(1:n_samples,function(samp){wrapper_to_sample_all_links(cl   = cl,
                                                                         n_mh = n_mh,
                                                                         Bstring    = test_bString,
                                                                         Gstring    = smc_samples$Gstring[samp,],
                                                                         p_link     = p_link,
                                                                         gCube      = smc_samples$gCube[samp,],
                                                                         nCube      = smc_samples$nCube[samp,],
                                                                         kCube      = smc_samples$kCube[samp,],
                                                                         inhib_inds = inhib_inds,
                                                                         model      = model,
                                                                         paramsList = paramsList,
                                                                         indexList  = indexList,
                                                                         sigma      = sigma,
                                                                         jump_size  = jump_size)})

    for (samp in 1:n_samples){
      smc_samples$gCube[samp,]   <- tmp[,samp]$gCube
      smc_samples$nCube[samp,]   <- tmp[,samp]$nCube
      smc_samples$kCube[samp,]   <- tmp[,samp]$kCube
      smc_samples$Gstring[samp,] <- tmp[,samp]$Gstring
    }

    old_post <- (sapply(1:n_samples, function(samp){
      posterior(cl,test_bString,
                    smc_samples$Gstring[samp,],
                    smc_samples$gCube[samp,],
                    smc_samples$nCube[samp,],
                    smc_samples$kCube[samp,],
                    p_link,inhib_inds,model,
                    paramsList,indexList,sigma)
      }))

    new_bString  <- add_link(init_bit_string=test_bString,links_mat=model$interMat)

    # If the graph is complete, exit,
    #    otherwise, add a new link to the graph
    if (all(new_bString == test_bString)){
      break
    } else{
      new_link <- which(new_bString - test_bString == 1)
      smc_samples$gCube[,new_link] <- runif(n_samples)
      smc_samples$nCube[,new_link] <- rexp(n_samples,1/2)
      smc_samples$kCube[,new_link] <- runif(n_samples)
      if (split_inhib){
        smc_samples$Gstring[,new_link+(0:(n_models-1))*n_params] <- rbinom(n_models*n_samples,1,p_link)
      }else{
        smc_samples$Gstring[,new_link] <- rbinom(n_samples,1,p_link)
      }

      test_bString <- new_bString
    }
  }

  if(n_cores>1) stopCluster(cl)
  if(diagnostics) print(w/sum(w))
  if(diagnostics) save(smc_samples,file='smc_samples.RData')

  smc_samples$version <- "v1.0"
  smc_samples
}

