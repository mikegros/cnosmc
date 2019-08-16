cno_smc <- function(n_samples, data, model,
                    init_links  = 1,
                    p_link      = 0.9,
                    n_mh        = 5,
                    jump_size   = rep(0.15,3),
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
  n_signals       <- length(data$namesSignals)

  test_bString    <- rep(0,n_params)

  # Finding correct indices for entries of gCube, etc
  simList <- prep4simFuzzy(model = model, paramsList = paramsList, verbose = FALSE)

  if ( ncol(simList$gCube) == 2 ){
    cube_inds <- matrix(0,nrow(simList$gCube) - n_params,2)

    for (jj in (1:nrow(simList$gCube))[-(1:n_params)]) {
      tmp1 <- strsplit(colnames(model$interMat)[jj],"=")[[1]]
      tmp2 <- strsplit(tmp1[1],"\\+")[[1]]
      cube_inds[jj-n_params,1] <- which(apply(sapply(colnames(model$interMat)[1:n_params],function(x)strsplit(x,"=")[[1]]),
                                              2,function(y) all(y == c(tmp2[1],tmp1[2]))))
      cube_inds[jj-n_params,2] <- which(apply(sapply(colnames(model$interMat)[1:n_params],function(x)strsplit(x,"=")[[1]]),
                                              2,function(y) all(y == c(tmp2[2],tmp1[2]))))

    }
  } else{
    cube_inds <- NA
  }

  # Continue with initializing model
  init_gCube <- 1*runif(n_params*n_samples)
  init_nCube <- rexp(n_params*n_samples,1/2) + 1
  init_kCube <- 1*runif(n_params*n_samples)
  init_sigsq <- 1/rgamma(n_samples*n_signals,1.25,10^-5)

  # Turn on some links for initial subgraph
  test_bString[c(init_links)] <- 1

  # Make sure a stimulus node is in the initial graph

  top_nodes <- which(apply(model$interMat[,1:n_params],1,function(x){all(x != 1)}))
  if (!any(model$interMat[,init_links,drop=FALSE][top_nodes,] == -1)) stop("Initial Graph Must Have A Stimulus Node!")

  n_models     <- 1
  init_Gstring <- rbinom(n_params*n_samples,1,p_link)

  smc_samples <- list(gCube = matrix(init_gCube, ncol = n_params, nrow = n_samples, byrow = TRUE))

  smc_samples$nCube   <- matrix(init_nCube,   ncol = n_params,          nrow = n_samples, byrow = TRUE)
  smc_samples$kCube   <- matrix(init_kCube,   ncol = n_params,          nrow = n_samples, byrow = TRUE)
  smc_samples$Gstring <- matrix(init_Gstring, ncol = n_models*n_params, nrow = n_samples, byrow = TRUE)
  smc_samples$sigsq   <- matrix(init_sigsq,   ncol = n_signals,         nrow = n_samples, byrow = TRUE)
  smc_samples$w       <- 1/n_samples


  colnames(smc_samples$Gstring) <- colnames(model$interMat)[1:n_params]

  old_post <- sapply(1:n_samples, function(samp){
    LogpriorGstring(smc_samples$Gstring[samp,],p_link) +
      Logpriorsigsq(smc_samples$sigsq[samp,],1.25,10^-5) +
      Logpriorg(smc_samples$gCube[samp,]) +
      Logpriorn(smc_samples$nCube[samp,]) +
      Logpriork(smc_samples$kCube[samp,]) })

  w <- rep(1,n_samples)/n_samples

  mh_jump_size <- data.frame(
    g = rep(jump_size[1],n_params),
    n = rep(jump_size[2],n_params),
    k = rep(jump_size[3],n_params)
  )

  for (stage in 1:n_params){

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
                                                      model      = model,
                                                      paramsList = paramsList,
                                                      indexList  = indexList,
                                                      simList    = simList,
                                                      cube_inds  = cube_inds)
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
                                                          model      = model,
                                                          paramsList = paramsList,
                                                          indexList  = indexList,
                                                          simList    = simList,
                                                          cube_inds  = cube_inds)
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
      smc_samples$sigsq   <- smc_samples$sigsq[resample_inds,]
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
                                                                                 model      = model,
                                                                                 paramsList = paramsList,
                                                                                 indexList  = indexList,
                                                                                 simList    = simList,
                                                                                 cube_inds  = cube_inds,
                                                                                 jump_size  = mh_jump_size)})
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
                                                                           model      = model,
                                                                           paramsList = paramsList,
                                                                           indexList  = indexList,
                                                                           simList    = simList,
                                                                           cube_inds  = cube_inds,
                                                                           jump_size  = mh_jump_size)})
    }

    # These variables are used to count how often mh transitions are accepted
    #     for changing jump_size between stages

    g_jump <- rep(0,n_params)
    n_jump <- rep(0,n_params)
    k_jump <- rep(0,n_params)

    for (samp in 1:n_samples){
      smc_samples$gCube[samp,]   <- tmp[,samp]$gCube
      smc_samples$nCube[samp,]   <- tmp[,samp]$nCube
      smc_samples$kCube[samp,]   <- tmp[,samp]$kCube
      smc_samples$sigsq[samp,]   <- tmp[,samp]$sigsq
      smc_samples$Gstring[samp,] <- tmp[,samp]$Gstring

      g_jump <- g_jump + tmp[,samp]$accepted$g
      n_jump <- n_jump + tmp[,samp]$accepted$n
      k_jump <- k_jump + tmp[,samp]$accepted$k
      if (n_mh == 0){
        old_post[samp] <- new_post[samp]
      } else{
        old_post[samp] <- tmp[,samp]$post
      }
    }

    # Update jump size
    mh_jump_size$g <- (g_jump + 1)/(n_samples*n_mh + 2) / 0.44 * mh_jump_size$g
    mh_jump_size$n <- (n_jump + 1)/(n_samples*n_mh + 2) / 0.44 * mh_jump_size$n
    mh_jump_size$k <- (k_jump + 1)/(n_samples*n_mh + 2) / 0.44 * mh_jump_size$k

    if(time_diagnostics) t3 <- proc.time() - t3
    if(time_diagnostics) print(paste('Time elapsed for MH step:',t3[3]/60))

    new_bString  <- add_link(init_bit_string=test_bString,links_mat=model$interMat)
    nParamsAdded <- length(which(new_bString==1)) - length(which(test_bString==1))

    if(diagnostics) print(paste('Number of parameters added:',nParamsAdded))

    # If the graph is complete, exit,
    #    otherwise, add a new link to the graph
    #    and populate the newly added g, n and k with samples from their priors
    if (all(new_bString == test_bString)){
      break
    } else{
      new_link <- which(new_bString - test_bString == 1)
      smc_samples$gCube[,new_link]   <- runif(n_samples*length(new_link))
      smc_samples$nCube[,new_link]   <- rexp(n_samples*length(new_link),1/2) + 1
      smc_samples$kCube[,new_link]   <- runif(n_samples*length(new_link))
      smc_samples$Gstring[,new_link] <- rbinom(n_samples*length(new_link),1,p_link)

      mh_jump_size$g[new_link] <- mean(mh_jump_size$g[test_bString == 1])
      mh_jump_size$n[new_link] <- mean(mh_jump_size$n[test_bString == 1])
      mh_jump_size$k[new_link] <- mean(mh_jump_size$k[test_bString == 1])

      active_nodes     <- sapply(strsplit(colnames(model$interMat)[which(test_bString==1)],split = "="),function(x){x[2]})
      active_nodes_new <- sapply(strsplit(colnames(model$interMat)[which(new_bString-test_bString==1)],split = "="),function(x){x[2]})
      active_nodes_new <- intersect(active_nodes,active_nodes_new)
      active_nodes_new <- which(paramsList$data$namesSignals %in% active_nodes_new)

      smc_samples$sigsq[,active_nodes_new] <- 1/rgamma(n_samples*length(active_nodes_new),1.25,10^-5)

      test_bString <- new_bString
    }
    if (checkpoint){
      print("Saving samples as checkpoint.")
      print(paste("at", Sys.time()))
      save(list=ls(),file = paste("smc_checkpoint_",stage,"_lots_of_stuff.RData",sep=""))
    }
  }

  if( n_cores > 1 ) stopCluster(cl)
  if(excess_cluster_call) {stopCluster(cl1)}

  if(diagnostics) print(w)

  smc_samples$mh_jump_size <- mh_jump_size

  smc_samples$version <- "v1.04"
  smc_samples
}

