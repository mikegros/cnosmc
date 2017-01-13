cno_smc <- function(n_samples, data, model,
                    init_links  = 1,
                    p_link      = 0.9,
                    n_mh        = 5,
                    jump_size   = rep(0.15,3),
                    split_inhib = FALSE,
                    n_cores     = 1,
                    diagnostics = FALSE,
                    checkpoint  = FALSE,
                    time_diagnostics    = FALSE,
                    excess_cluster_call = FALSE){
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

    inhib_inds <- apply(inhib_settings,1,function(y) which(apply(data$valueInhibitors,1,function(x) all(x==y))))

    if(ncol(inhib_inds) == 1) {
      tmp             <- c(inhib_inds)
      inhib_inds      <- list()
      inhib_inds[[1]] <- tmp
    }
  }

  paramsList          <- defaultParametersFuzzy(data, model)
  indexList           <- indexFinder(CNOlist = data, model = model)

  if( n_cores > 1 ){
    # outfile arguement tells R where to write print console output (print statements, etc) when running in parallel
    # cl <- makeCluster(2,outfile="~/output_R.txt")

    cl <- makeCluster(n_cores)
    cl1 <- NULL
    if(excess_cluster_call) cl1 <- cl

    clusterCall(cl,function(x) {library(CNORfuzzy)})
    if(excess_cluster_call) clusterCall(cl1,function(x) {library(CNORfuzzy)})

  }else{
    cl1 <- NULL
  }


  findMainEffects <- unlist(lapply(strsplit(model$reacID,'\\+'),function(x) {length(x) < 2} ))
  n_params        <- sum(findMainEffects)

  test_bString    <- rep(0,n_params)

  init_gCube <- 1*runif(n_params*n_samples)
  init_nCube <- rexp(n_params*n_samples,1/2)
  init_kCube <- 1*runif(n_params*n_samples)
  init_sigsq <- 1/rgamma(n_samples,1.25,10^-5)

  determine_add_list <- list(c(31L, 39L, 47L),
                             c(33L, 46L),
                             c(19L, 32L, 44L),
                             c(6L, 28L, 36L, 43L),
                             c(15L, 27L, 51L),
                             c(2L, 38L, 40L, 50L),
                             c(25L, 37L),
                             c(3L, 26L, 45L),
                             34L,
                             c(18L, 23L),
                             c(10L, 12L, 21L),
                             24L,
                             4L,
                             c(5L, 13L),
                             41L, 14L, 29L,  8L, 30L, 42L, 22L, 52L, 48L,
                             11L,  9L, 17L, 49L, 53L, 20L, 35L,  1L,  7L, 16L)

  # Turn on some links for initial subgraph
  init_links               <- determine_add_list[[1]]
  test_bString[init_links] <- 1
  # Make sure a stimulus node is in the initial graph

  top_nodes <- which(apply(model$interMat[,1:n_params],1,function(x){all(x != 1)}))
  if (!any(model$interMat[,init_links,drop=FALSE][top_nodes,] == -1)) stop("Initial Graph Must Have A Stimulus Node!")

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
  smc_samples$sigsq   <- matrix(init_sigsq,   ncol = 1,                 nrow = n_samples, byrow = TRUE)
  smc_samples$w       <- 1/n_samples


  if (split_inhib){
    colnames(smc_samples$Gstring) <- rep(colnames(model$interMat)[1:n_params],n_models)
  }else{
    colnames(smc_samples$Gstring) <- colnames(model$interMat)[1:n_params]
  }

  old_post <- sapply(1:n_samples, function(samp){
    LogpriorGstring(smc_samples$Gstring[samp,],p_link) +
      Logpriorsigsq(smc_samples$sigsq[samp,],1.25,10^-5) +
      Logpriorg(smc_samples$gCube[samp,]) +
      Logpriorn(smc_samples$nCube[samp,]) +
      Logpriork(smc_samples$kCube[samp,]) })

  w <- rep(1,n_samples)/n_samples

  for (stage in 1:33){
    print(paste("Stage: ",stage,sep=""))
    if(diagnostics) save(smc_samples,file='smc_samples.RData')

    # Get likelihood weights
    if( n_cores > 1 ) clusterExport(cl,varlist=ls(),envir = environment())
    if(excess_cluster_call) clusterExport(cl1,varlist=ls(),envir = environment())

    if(time_diagnostics) t1 <- proc.time()
    if( n_cores > 1 & !excess_cluster_call ){
      new_post <- (parSapply(cl,1:n_samples, function(samp){
                                            posterior(cl         = cl1,
                                                      Bstring    = test_bString,
                                                      Gstring    = smc_samples$Gstring[samp,],
                                                      gCube      = smc_samples$gCube[samp,],
                                                      nCube      = smc_samples$nCube[samp,],
                                                      kCube      = smc_samples$kCube[samp,],
                                                      sigsq      = smc_samples$sigsq[samp,],
                                                      p_link     = p_link,
                                                      inhib_inds = inhib_inds,
                                                      model      = model,
                                                      paramsList = paramsList,
                                                      indexList  = indexList)
      }))
    } else {
      new_post <- (sapply(1:n_samples, function(samp){
                                                posterior(cl         = cl1,
                                                          Bstring    = test_bString,
                                                          Gstring    = smc_samples$Gstring[samp,],
                                                          gCube      = smc_samples$gCube[samp,],
                                                          nCube      = smc_samples$nCube[samp,],
                                                          kCube      = smc_samples$kCube[samp,],
                                                          sigsq      = smc_samples$sigsq[samp,],
                                                          p_link     = p_link,
                                                          inhib_inds = inhib_inds,
                                                          model      = model,
                                                          paramsList = paramsList,
                                                          indexList  = indexList)
      }))
    }

    w  <- log(w) + new_post - old_post
    w  <- w - max(w)
    w  <- exp(w)
    w  <- w/sum(w)

    smc_samples$w <- w

    if(time_diagnostics) t1 <- proc.time() - t1
    if(time_diagnostics) print(paste('Time elapsed for calculating weights:',round(t1[3]/60,3),'minutes'))

    ESS <- 1/(sum(w^2))
    monitorESS <- list(ESS,n_samples/2)
    if(diagnostics) print(paste('ESS:',monitorESS))
    if(diagnostics) save(monitorESS,file=paste('ESS',stage,'.RData'))

    # Resample parameters using likelihood weights
    if(time_diagnostics) t2 <- proc.time()
    if (ESS < n_samples/2){
      resample_inds <- sample(1:n_samples,n_samples,replace = TRUE,prob=w)

      smc_samples$gCube   <- smc_samples$gCube[resample_inds,]
      smc_samples$nCube   <- smc_samples$nCube[resample_inds,]
      smc_samples$kCube   <- smc_samples$kCube[resample_inds,]
      smc_samples$Gstring <- smc_samples$Gstring[resample_inds,]
      smc_samples$sigsq   <- smc_samples$sigsq[resample_inds,,drop=FALSE]
      smc_samples$w       <- rep(1,n_samples)/n_samples

      w <- rep(1,n_samples)/n_samples
    }

    if(time_diagnostics) t2 <- proc.time() - t2
    if(time_diagnostics) print(paste('Time elapsed for resampling:',round(t2[3]/60,3),'minutes'))
    # Perturb resampled values with an MH step
    if(time_diagnostics) t3 <- proc.time()
    if( n_cores > 1 ) clusterExport(cl,varlist=ls(),envir = environment())
    if(excess_cluster_call) clusterExport(cl1,varlist=ls(),envir = environment())

    if( n_cores > 1 & !excess_cluster_call ){
      tmp <- parSapply(cl,1:n_samples,function(samp){wrapper_to_sample_all_links(cl   = cl1,
                                                                                 n_mh = n_mh,
                                                                                 Bstring    = test_bString,
                                                                                 Gstring    = smc_samples$Gstring[samp,],
                                                                                 p_link     = p_link,
                                                                                 gCube      = smc_samples$gCube[samp,],
                                                                                 nCube      = smc_samples$nCube[samp,],
                                                                                 kCube      = smc_samples$kCube[samp,],
                                                                                 sigsq      = smc_samples$sigsq[samp,],
                                                                                 inhib_inds = inhib_inds,
                                                                                 model      = model,
                                                                                 paramsList = paramsList,
                                                                                 indexList  = indexList,
                                                                                 jump_size  = jump_size)})
    }else{
      tmp <- sapply(1:n_samples,function(samp){wrapper_to_sample_all_links(cl   = cl1,
                                                                           n_mh = n_mh,
                                                                           Bstring    = test_bString,
                                                                           Gstring    = smc_samples$Gstring[samp,],
                                                                           p_link     = p_link,
                                                                           gCube      = smc_samples$gCube[samp,],
                                                                           nCube      = smc_samples$nCube[samp,],
                                                                           kCube      = smc_samples$kCube[samp,],
                                                                           sigsq      = smc_samples$sigsq[samp,],
                                                                           inhib_inds = inhib_inds,
                                                                           model      = model,
                                                                           paramsList = paramsList,
                                                                           indexList  = indexList,
                                                                           jump_size  = jump_size)})
    }

    for (samp in 1:n_samples){
      smc_samples$gCube[samp,]   <- tmp[,samp]$gCube
      smc_samples$nCube[samp,]   <- tmp[,samp]$nCube
      smc_samples$kCube[samp,]   <- tmp[,samp]$kCube
      smc_samples$sigsq[samp,]   <- tmp[,samp]$sigsq
      smc_samples$Gstring[samp,] <- tmp[,samp]$Gstring

      old_post[samp]             <- tmp[,samp]$post
    }

    if(time_diagnostics) t3 <- proc.time() - t3
    if(time_diagnostics) print(paste('Time elapsed for MH step:',t3[3]/60))


    # Adding steps deterministically for now:
    new_bString  <- test_bString

    if(stage < 33) new_bString[determine_add_list[[stage+1]]] <- 1
    nParamsAdded <- length(which(new_bString==1)) - length(which(test_bString==1))

    if(diagnostics) print(paste('Number of parameters added:',nParamsAdded))

    # If the graph is complete, exit,
    #    otherwise, add a new link to the graph
    #    and populate the newly added g, n and k with samples from their priors
    if (all(new_bString == test_bString)){
      break
    } else{
      new_link <- which(new_bString - test_bString == 1)
      for (links in new_link){
        smc_samples$gCube[,new_link] <- runif(n_samples)
        smc_samples$nCube[,new_link] <- rexp(n_samples,1/2)
        smc_samples$kCube[,new_link] <- runif(n_samples)
        if (split_inhib){
          smc_samples$Gstring[,new_link+(0:(n_models-1))*n_params] <- rbinom(n_models*n_samples,1,p_link)
        }else{
          smc_samples$Gstring[,new_link] <- rbinom(n_samples,1,p_link)
        }
      }

      test_bString <- new_bString
    }
    if (checkpoint){
      print("Saving samples as checkpoint.")
      print(paste("at", Sys.time()))
      tryCatch(save(smc_samples, test_bString, old_post, tmp,
           n_samples, data, model,init_links,p_link,n_mh,jump_size,
                    split_inhib,n_cores,diagnostics,checkpoint,new_bString, paramsList,
                    time_diagnostics,excess_cluster_call,nParamsAdded,ESS, indexList,
           file = paste("smc_checkpoint_",stage,"_lots_of_stuff.RData",sep="")),
               save(smc_samples,file = paste("smc_checkpoint_",stage,".RData",sep="")) )
    }
  }

  if( n_cores > 1 ) stopCluster(cl)
  if(excess_cluster_call) {stopCluster(cl1)}

  if(diagnostics) print(w)

  smc_samples$version <- "vNew_big_example"

  smc_samples
}

