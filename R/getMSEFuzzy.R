getMSEFuzzy = function(cl1     = NULL,
                       Bstring    = rep(1, length(model$reacID)+length(grep("\\+",model$reacID))),
                       Gstring    = rep(1,  length(model$reacID[-grep("\\+",model$reacID)])),
                       gCube      = rep(4,  length(model$reacID[-grep("\\+",model$reacID)])),
                       nCube      = rep(.5, length(model$reacID[-grep("\\+",model$reacID)])),
                       kCube      = rep(.2, length(model$reacID[-grep("\\+",model$reacID)])),
                       sigsq,
                       sizeFac    = 0,
                       NAFac      = 0,
                       verbose    = FALSE,
                       model,
                       paramsList,
                       indexList,
                       simList,
                       cube_inds){


  #
  #
  #
  # CHANGELOG:
  #
  # [Jan 20, 2017] Changes
  #         Now sigsq is a vector but if it is a scalar then we make it into a vector
  #         output nDataP is now a vector of observation counts for everything, not just the active_nodes.
  #         active_nodes is a new output of indices of the observed columns of the data matrix.
  #
  # [April 26, 2016] Changes
  #         Renamed file for clarity of use with respect to the R package
  #         Added to R package
  #
  #
  # [April 12, 2016] Potential issues:
  #         getMSEfuzzy -- Since we replace NAs in the simulated data with 0 we now the same
  #                        with NA values in the raw data when calculating the SSE
  #
  # [April 12, 2016] Changes
  #           getMSEfuzzyDataList  =  evaluate the likelihood with partial data
  #           SplitDataByInhibitor = split the data file into a list with elements
  #                                  taking on different settings of the inhibitor
  #           Gstring was added as inputs to all the getMSE functions, the default
  #                                   value of Gstring is set to a vector of 1s
  #           Gathered additional functions from other sources that we have written
  #           Merged getMSE functions to use alternative versions of the code
  #           getMSE functions use clusterCall which is now turned off by leaving the cluster
  #                                   name at its default value of NULL
  #           Also removed outputs that were multiply named across different versions of the code
  #
  # [Jan 24, 2016] - generalized code that pointed to the Morris specific model rather than the model of interest
  #
  # [Nov 1, 2015] - changed the number of parameters.
  # We only need enough parameters for the main effects.
  # The interaction functions are automatically populated by main effect values
  # because the interactions are always MAX, MIN, AND, OR
  # So now length(gparameters) = length(nparameters) = length(kparameters) = length(# main effects)
  # Bstring remains of length (# main effects + # interactions)
  #
  #
  #Functions in this file:
  # getMSEfuzzy          = given a set of parameters obtain the SSE and MSE (negative log Gaussian likelihood)
  # getMSEfuzzyMoreInfo  = same as getMSEfuzzy except that it opens up some of the library(CNORfuzzy)
  #                        functions to extract the number of NA points in the simulated data and the
  #                        size of the network.  The size and # NAs are used as penalties in the library(CNORfuzzy)
  #                        optimization.  Because this file opens up functions I suggest using getMSEfuzzy instead
  #                        but wrote this code for future potential extraction of penalty information and diagnostic purposes
  # Sample1FromPriorGetNegLogLike = Sample 1 parameter vector for the network from the set of priors (BString,g,n,k)
  # ApplyFriendlyMSE              = This takes a vector of all parameters (BString,g,n,k) for one run of
  #                                 the network model and then returns only the negative log likelihood.  This is
  #                                 a wrapper for getMSEfuzzy so that it works nicely with
  #                                 the 'apply' family of functions which can then be performed over large number of
  #                                 parameter vectors
  # SplitDataByInhibitor     = Split the data into a list where each component of the
  #                            list uses a different setting of the inhibitor
  # getMSEfuzzyDataList      = This function is basically identical to getMSEfuzzy except that the input
  #                            ParamsList contains a dataList and not just a data file AND GString is
  #                            conditional on the inhibitor value meaning that a new value of Gstring is passed
  #                             for each inhibitor setting (and equivalently, each element of the dataList)
  #
  #
  # generatePlot      = plot the model with highlights for the best model
  #                      and probabilities for each of the links
  # SelectModel       = get the table with probabilities for each model
  #
  #
  #
  #
  #CALLED BY:  after SMC is run
  #
  #
  ########
  #
  # Last update: April 12, 2016, Dave Campbell
  #
  #
  # [Apr 12] folded in Gstring and clusterCall from Biljana's version
  # [Nov 1 ] adjusted the number of parameters for each of (g,n,k)
  #
  # This uses the 'fuzzy' logic model as opposed to the 'boolean' logic model.
  # The 'boolean' model just look at links as being 'on' or 'off'
  # the 'fuzzy' model refines the boolean by including the parameters (g,n,k) for all links
  #
  # Use this function to:
  # Given a set of parameters obtain the MSE = SSE/ (number of observation points)
  #
  #
  # CALLED BY:
  #
  #
  #==============================================================
  #                         REQUIRED INPUTS
  #==============================================================
  # cl1 = cluster as made by library(parallel);cl1 = makeCluster(4);
  #       clusterCall(cl1,function(x) {library(CNORfuzzy)});
  #       clusterExport(cl1,varlist=ls(),envir = environment())
  #     If the cluster is not being used leave it at its default of cl1=NULL
  #
  # Bstring = binary vector of 1(link is on) and 0(link is off)
  #
  # Gstring = multiplier for gCube to allow the parameter to be discrete (0) or continuous (0,\infty]
  #
  # The next 3 inputs are matrices of parameters with 2 columns
  # Reactions are row elements
  # When a reaction is a main effect, just the first column of parameters are used
  # When a reaction is a (2way) interaction, the 2 columns are the parameters for each input in the interaction
  #
  # gCube = matrix of g parameters for all of the links.  Used only if the corresponding main effect element of Bstring == 1
  # nCube = matrix of n parameters for all of the links.  Used only if the corresponding main effect element of Bstring == 1
  # kCube = matrix of k parameters for all of the links.  Used only if the corresponding main effect element of Bstring == 1
  #
  #
  # model = the reduced model after the prior knowledge network is trimmed down to exclude un-identifiable links
  #         This is fixed throughout.  Individual links are then turned on and off using Bstring
  #
  # paramsList = is a very long list that can be left at all of its default parameter values - they are mainly used
  #            for the genetic algorithm optimizer.
  #
  # indexList = the output from the function 'indexFinder'
  #
  # sizeFac = is effectively a penalty on network size.
  # NAFac = is a penalty on the number of NA values returned.
  #
  # verbose = logical value; 'TRUE' prints messages reminding you what the sizeFac and NAFac parameters do.
  #
  #===============================================================%
  #                           OUTPUTS                             %
  #===============================================================%
  # getMSEfuzzy: A list with elements:
  # model = its input value
  # MSE = Residual SSE / (# observations), note that the number of observations is the total number and not
  #                                             the number of included observations in the sub graph.
  # SSE = Residual SSE
  # NAFac   = its input value
  # sizeFac = its input value
  # SimResults = SimResults
  #===============================================================%
  #                          MODEL DETAILS                        %
  #===============================================================%
  # COMMENTS: The software evaluates the model in 2 stages, first it uses Bstring to
  #           build the discrete Boolean version of the model using 'interpretDiscreteGA'.
  #           Then it includes the model parameters to obtain the model fit using 'simFuzzyT1'
  #===============================================================%

  if(verbose){
    print("note that SSE = deviationPen + NAPen + sizePen")
    print(paste("Where: NAPen <- NAFac * length(which(is.na(simResults)) but NAFac = ",NAFac))
    print(paste("And: sizePen <- (nDataPts * sizeFac * nInputs)/nInTot, but sizeFac = ",sizeFac))
  }

  n_params <- length(Bstring)

  #############
  # Find the indices of observed nodes that are in the active subgraph
  active_nodes <- sapply(strsplit(colnames(model$interMat)[which(Bstring==1)],split = "="),function(x){x[2]})
  active_nodes <- which(paramsList$data$namesSignals %in% active_nodes)

  #############

  # NO LONGER CALL THIS HERE. IT GETS WRITTEN OVER AND EATS UP TIME
  # simList <- prep4simFuzzy(model = model, paramsList = paramsList, verbose = FALSE)

  ####################Put model parameters into the appropriate list <--THIS IS THE FUNCTION CONDITIONAL ON THE PRESENCE/ABSENCE OF A LINK
  # This was stated to be Something x 2 matrix, but was not actually....

  if ( ncol(simList$gCube) == 2 ) {
    simList$gCube[1:n_params,1]    <- gCube # matrix of __ rows, 2 columns
    simList$gCube[-(1:n_params),1] <- gCube[cube_inds[,1]] # matrix of __ rows, 2 columns
    simList$gCube[-(1:n_params),2] <- gCube[cube_inds[,2]] # matrix of __ rows, 2 columns

    simList$nCube[1:n_params,1]    <- nCube                # matrix of __ rows, 2 columns
    simList$nCube[-(1:n_params),1] <- nCube[cube_inds[,1]] # matrix of __ rows, 2 columns
    simList$nCube[-(1:n_params),2] <- nCube[cube_inds[,2]] # matrix of __ rows, 2 columns

    simList$kCube[1:n_params,1]    <- kCube                # matrix of __ rows, 2 columns
    simList$kCube[-(1:n_params),1] <- kCube[cube_inds[,1]] # matrix of __ rows, 2 columns
    simList$kCube[-(1:n_params),2] <- kCube[cube_inds[,2]] # matrix of __ rows, 2 columns

    tmpSimList <- cutSimList(simList, Bstring*Gstring)

    trim_inds <- rownames(simList$finalCube) %in% rownames(tmpSimList$finalCube)
    tmpSimList$gCube <- tmpSimList$gCube[trim_inds,]
    tmpSimList$nCube <- tmpSimList$nCube[trim_inds,]
    tmpSimList$kCube <- tmpSimList$kCube[trim_inds,]

  } else{
    simList$gCube <- matrix(gCube*Bstring*Gstring,ncol=1)
    simList$nCube <- matrix(nCube,ncol=1)
    simList$kCube <- matrix(kCube,ncol=1)

    tmpSimList <- cutSimList(simList, Bstring*Gstring)

    trim_inds <- rownames(simList$finalCube) %in% rownames(tmpSimList$finalCube)
    tmpSimList$gCube <- tmpSimList$gCube[trim_inds,]
    tmpSimList$nCube <- tmpSimList$nCube[trim_inds,]
    tmpSimList$kCube <- tmpSimList$kCube[trim_inds,]

  }

  # Replicate 4 times to attempt to avoid the stochastic failure issue observed previously
  if(!is.null(cl1)){
    tmpModel <- cutModel(model, Bstring*Gstring)
    SimResultsList <- clusterCall(cl1,function() {simFuzzyT1(CNOlist = paramsList$data,
                                                             model   = tmpModel,
                                                             simList = tmpSimList)})
    SimResultsList <- array(unlist(SimResultsList), dim = c(nrow(SimResultsList[[1]]), ncol(SimResultsList[[1]]), length(SimResultsList)))
    SimResultsList[is.na(SimResultsList)] = 0
    # Obtain the SSEs from the simulated results and the data
    # NOTE: Occasionally this returns NA values for some parameter settings. In these cases SimResultsList seems to give
    #       back a blank matrix of values. May be an issue with the apply statement and the structure of the output from
    #       "replicate", but I think it is just a result of being at bad places in parameter space. Will explore this more
    #       to ensure that the SSE is being handled properly
    # Since we remove the NA values from the simulated data we may as well do the same from the real data
    if(length(sigsq)==1){
        sigsq = rep(sigsq,dim(paramsList$data$valueSignals[[2]])[2])
    }
    SumSqaredResidsByRun =
                apply(SimResultsList,3,function(x) {#Get a vector of sums for each 3rd dimension element
                       apply((x[,indexList$signals[active_nodes]] -
                       paramsList$data$valueSignals[[2]][,active_nodes])^2, 2, sum,na.rm=TRUE)#/sigsq[active_nodes]# sum down columns
                     })

    SSEvectorScaled = rep(0,dim(paramsList$data$valueSignals[[2]])[2])
    SSEvectorScaled[active_nodes] = SumSqaredResidsByRun[, which.min(apply(SumSqaredResidsByRun,2,sum))]
  }else{
    tmpModel <- cutModel(model, Bstring*Gstring)
    SimResultsList <- simFuzzyT1(CNOlist = paramsList$data,
                                      model   = tmpModel,
                                      simList = tmpSimList)
    if(length(sigsq)==1){
        sigsq = rep(sigsq,dim(paramsList$data$valueSignals[[2]])[2])
    }

     SSEvectorScaled = rep(0,dim(paramsList$data$valueSignals[[2]])[2])
    # SSEvectorScaled[active_nodes] = apply((SimResultsList[,indexList$signals[active_nodes],drop=FALSE] -
    #                                          paramsList$data$valueSignals[[2]][,active_nodes,drop=FALSE])^2, 2, sum,na.rm=TRUE)# /sigsq[active_nodes]
     SSEvectorScaled[active_nodes] = apply((SimResultsList[,indexList$signals[active_nodes],drop=FALSE] -
                                              paramsList$data$valueSignals[[2]][,active_nodes,drop=FALSE])^2, 2, sum)# /sigsq[active_nodes]
     SSEvectorScaled[is.na(SSEvectorScaled)] <- 1.e10

  }

  # NOTE: Setting all NA values to 0. This is just a hack to avoid killing particles that might otherwise be at a reasonable
  #       place in the parameter space. Need to figure out why the NA's happen in their code. Could be more correct to kill these
  #       particles or add a more reasonble NA penalty.

  ################

  # Obtain the SSEs from the simulated results and the data
  # NOTE: Occasionally this returns NA values for some parameter settings. In these cases SimResultsList seems to give
  #       back a blank matrix of values. May be an issue with the apply statement and the structure of the output from
  #       "replicate", but I think it is just a result of being at bad places in parameter space. Will explore this more
  #       to ensure that the SSE is being handled properly
  # Since we remove the NA values from the simulated data we may as well do the same from the real data


  # Return Inf for the SSE and MSE if the SimResults matrix is bad
  if(is.na(sum(SSEvectorScaled))){
    return(list(model = model,MSE=Inf,SSEvectorScaled=Inf,NAFac=NAFac,sizeFac=sizeFac,SimResults=SimResults,nDataP = 1,SSE=Inf))
  }else{
    nDataP = colSums(!is.na(paramsList$data$valueSignals[[2]]))
    MSE    = SSEvectorScaled/nDataP
    return(list(model = model,MSE=MSE,NAFac=NAFac,nDataP = nDataP,
                sizeFac=sizeFac,SimResults=SimResultsList,
                SSEvectorScaled=SSEvectorScaled, active_nodes=active_nodes))
  }
}
