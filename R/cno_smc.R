cno_smc <- function(n_samples, data, model,
                    init_links  = 1,
                    p_link      = 0.9,
                    datalst     = NULL,
                    n_cores     = 1,
                    diagnostics = F){
  #
  # Last update: April 26, 2016, Mike Grosskopf
  #
  # [Apr 26,2016] Mike Grosskopf
  #               Took SMC code from individual scripts and made it a standalone file.
  #                 Did tons of code cleaning to get it completely working and clean as well
  #

  if (!is.null(datalst)){
    dataList=splitDataByInhibitor(data)}

  if (!is.null(datalst)){
    #changed this o accomodate data list
    paramsList          <- lapply(1:length(dataList), function(x) {defaultParametersFuzzy(dataList[[x]], model)})
    indexList           <- lapply(1:length(dataList), function(x) {indexFinder(CNOlist = dataList[[x]], model = model)})
    paramsList$dataList = SplitDataByInhibitor(data)
  }else{
    paramsList          <- defaultParametersFuzzy(data, model)
    indexList           <- indexFinder(CNOlist = data, model = model)
  }

  if( n_cores > 1 ){
    # outfile arguement tells R where to write print console output (print statements, etc) when running in parallel
    # cl <- makeCluster(2,outfile="~/output_R.txt")
    cl <- makeCluster(4)
    clusterCall(cl,function(x) {library(CNORfuzzy)})
  }else{
    cl <- NULL
  }


  findMainEffects <- unlist(lapply(strsplit(model$reacID,'\\+'),function(x) {length(x)<2} ))
  n_params        <- length(which(findMainEffects==T))

  test_bString    <- rep(0,n_params)

  init_gCube <- 1*runif(n_params*n_samples)
  init_nCube <- rexp(n_params*n_samples,1/0.2)
  init_kCube <- 1*runif(n_params*n_samples)

  # Turn on some links for initial subgraph
  test_bString[c(init_links)] <- 1

  if (!is.null(datalst)){
    nmodels=length(paramsList$dataList)
    init_Gstring <- rbinom(nmodels*n_params*n_samples,1,p_link)
  }else{
    init_Gstring <- rbinom(n_params*n_samples,1,p_link)
  }
  smc_samples <- list(gCube = matrix(init_gCube, ncol = n_params, nrow = n_samples, byrow = TRUE))

  smc_samples$nCube   <- matrix(init_nCube, ncol = n_params, nrow = n_samples, byrow = TRUE)
  smc_samples$kCube   <- matrix(init_kCube, ncol = n_params, nrow = n_samples, byrow = TRUE)
  smc_samples$Gstring <- matrix(init_Gstring, ncol =3*n_params, nrow = n_samples, byrow = TRUE)

  if (!is.null(datalst)){
    colnames(smc_samples$Gstring)=rep(colnames(model$interMat)[1:n_params],nmodels)
  }else{
    colnames(smc_samples$Gstring)=colnames(model$interMat)[1:n_params]
  }

  for (stage in 1:n_params){

    print(stage)

    # Get likelihood weights
    if(n_cores>1) clusterExport(cl,varlist=ls(),envir = environment())

    w <- (sapply(1:n_samples, function(samp){
      if (!is.null(datalst)){
        (-1/2)*getMSEfuzzyDataList(cl=cl,
                                   Bstring    = test_bString,
                                   Gstring    = smc_samples$Gstring[samp,],
                                   gCube      = smc_samples$gCube[samp,],
                                   nCube      = smc_samples$nCube[samp,],
                                   kCube      = smc_samples$kCube[samp,],
                                   model      = model,
                                   paramsList = paramsList,
                                   indexList  = indexList,
                                   sizeFac    = 0,NAFac=0,verbose = FALSE)$SSE/0.1^2
      }else{
        (-1/2)*getMSEfuzzy(cl=cl,
                           Bstring    = test_bString,
                           Gstring    = smc_samples$Gstring[samp,],
                           gCube      = smc_samples$gCube[samp,],
                           nCube      = smc_samples$nCube[samp,],
                           kCube      = smc_samples$kCube[samp,],
                           model      = model,
                           paramsList = paramsList,
                           indexList  = indexList,
                           sizeFac    = 0,NAFac=0,verbose = FALSE)$SSE/0.1^2

      }

    }))

    w  <- w-max(w)
    w  <- exp(w)
    we <- w/sum(w)

    if(diagnostics) save(we,file=paste('weights',stage,'2.RData',sep=''))
    if(diagnostics) save(smc_samples,file=paste('smc_samples_before_resample',stage,'2.RData',sep=''))

    # Resample parameters using likelihood weights
    resample_inds <- sample(1:n_samples,n_samples,replace = TRUE,prob=w)

    smc_samples$gCube   <- smc_samples$gCube[resample_inds,]
    smc_samples$nCube   <- smc_samples$nCube[resample_inds,]
    smc_samples$kCube   <- smc_samples$kCube[resample_inds,]
    smc_samples$Gstring <- smc_samples$Gstring[resample_inds,]

    if(diagnostics) save(smc_samples,file=paste('smc_samples_resampled_',stage,'2.RData'))

    # Perturb resampled values with an MH step

    if(n_cores>1) clusterExport(cl,varlist=ls(),envir = environment())
    tmp <- sapply(1:n_samples,function(samp){wrapper_to_sample_all_links(cl=cl,
                                                                         Bstring    = test_bString,
                                                                         Gstring    = smc_samples$Gstring[samp,],
                                                                         p_link     = p_link,
                                                                         gCube      = smc_samples$gCube[samp,],
                                                                         nCube      = smc_samples$nCube[samp,],
                                                                         kCube      = smc_samples$kCube[samp,],
                                                                         model      = model,
                                                                         paramsList = paramsList,
                                                                         indexList  = indexList,
                                                                         jump_size  = c(0.15,0.15,0.15),
                                                                         datalst    = datalst)})

    for (samp in 1:n_samples){
      smc_samples$gCube[samp,]   <- tmp[,samp]$gCube
      smc_samples$nCube[samp,]   <- tmp[,samp]$nCube
      smc_samples$kCube[samp,]   <- tmp[,samp]$kCube
      smc_samples$Gstring[samp,] <- tmp[,samp]$Gstring
    }

    if(diagnostics) save(smc_samples,file=paste('smc_samples_jittered_',stage,'2.RData',sep=''))

    new_bString  <- add_link(init_bit_string=test_bString,links_mat=model$interMat)

    # If the graph is complete, exit,
    #    otherwise, add a new link to the graph
    if (all(new_bString == test_bString)){
      break
    } else{
      test_bString <- new_bString
    }
  }

  if(n_cores>1) stopCluster(cl)
  if(diagnostics) print(w/sum(w))
  if(diagnostics) save(smc_samples,file='smc_samples.RData')

  smc_samples
}

