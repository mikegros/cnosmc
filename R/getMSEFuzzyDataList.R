########################################################
##############  getMSEfuzzyDataList
########################################################

getMSEFuzzyDataList = function(cl1=NULL,Bstring = rep(1, length(model$reacID[-grep("\\+",model$reacID)])),
                               Gstring =    rep(1,  length(model$reacID[-grep("\\+",model$reacID)])*length(paramsList$dataList)),
                               gCube   =    rep(4,  length(model$reacID[-grep("\\+",model$reacID)])),
                               nCube   =    rep(.5, length(model$reacID[-grep("\\+",model$reacID)])),
                               kCube   =    rep(.2, length(model$reacID[-grep("\\+",model$reacID)])),
                               model, paramsList,indexList,
                               sizeFac = 0,
                               NAFac   = 0,
                               verbose = TRUE){


  #Last update:  April 14 2016 Biljana (very small changes)
  # In the output list   changed  MSEpiece= MSE to  MSEpiece= MSEfuzzy$MSE
  #also in  temp = getMSEfuzzy(...) changed
  #         temp = getMSEfuzzy(..,indexList=indexList,..) into
  #         temp = getMSEfuzzy(..,indexList =indexList[[kp]],..)

  ########
  # Last update: April 12, 2016, Dave Campbell
  #
  #
  # Use this function to:
  # Get the MSE for a set of parameters and a datalist
  # The datalist is a list of valid library(CNORfuzzy) data
  # elements.  These could be from different experiments
  # with shared parameters or something like that.
  #
  # This function is basically identical to getMSEfuzzy except that the input
  # ParamsList contains a dataList and not just a data file AND GString is
  # conditional on the inhibitor value meaning that a new value of Gstring is passed
  # for each inhibitor setting (and equivalently, each element of the dataList)
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
  #            for the genetic algorithm optimizer.  The dataList is an element within this
  #
  # indexList = the output from the function 'indexFinder'
  #
  # sizeFac = is effectively a penalty on network size.
  # NAFac = is a penalty on the number of NA values returned.
  #
  # verbose = logical value; 'TRUE' prints messages reminding you what the sizeFac and NAFac parameters do.
  #

  #
  #===============================================================%
  #                           OUTPUTS                             %
  #===============================================================%
  #
  # a list of valid data files where each takes only a single inhibitor setting
  #
  #===============================================================%
  #                          MODEL DETAILS                        %
  #===============================================================%
  # COMMENTS:
  #===============================================================%

  # Possibly split the data by inihibitor
  if(is.null(dataList)){
    dataList = paramsList$dataList = SplitDataByInhibitor(paramsList$data)
  }
  # make sure the pieces are in the right places
  if(is.null(paramsList$dataList)){
    paramsList$dataList = dataList
  }
  # make sure that Gstring is present
  if(is.null(Gstring)){
    #Gstring = rep(1,length(model$reacID[-grep("\\+",model$reacID)])*length(paramsList$dataList))
    Gstring = rep(1,length(model$reacID)*length(paramsList$dataList))
  }
  # and that it is the correct length:
  if(length(Gstring)==length(model$reacID[-grep("\\+",model$reacID)])){
    Gstring = rep(Gstring,length(paramsList$dataList))
  }else if(length(Gstring)==0){
    #Gstring =    rep(1,  length(model$reacID[-grep("\\+",model$reacID)])*length(paramsList$dataList))
    Gstring =    rep(1,  length(model$reacID)*length(paramsList$dataList))
  }


  paramsListUse       =  paramsList
  #L = length(model$reacID[-grep("\\+",model$reacID)])
  L = length(model$reacID)

  MSEfuzzy            = list()
  MSEfuzzy$model      = model
  MSEfuzzy$MSE        = c()
  MSEfuzzy$NAFac      = NAFac
  MSEfuzzy$sizeFac    = sizeFac
  MSEfuzzy$SimResultspieces = list()
  MSEfuzzy$SSE        = c()
  MSEfuzzy$nDataP     = c()
  MSEfuzzy$SimResults = c()

  for(kp in (1:length(paramsList$dataList)) ){
    # since clusterCall is used from within getMSEfuzzy we can't parallelize
    # this embarassingly parallel task
    paramsListUse       =  paramsList[[kp]]
    paramsListUse$data  =  paramsList$dataList[[kp]]
    temp                =  getMSEfuzzy(cl1,Bstring = Bstring,
                                       Gstring = Gstring[L*(kp-1)+(1:L)],
                                       gCube = gCube,
                                       nCube = nCube,
                                       kCube = kCube,
                                       model = model, paramsList = paramsListUse,
                                       indexList=indexList[[kp]], sizeFac = sizeFac,
                                       NAFac=NAFac,verbose = verbose)

    MSEfuzzy$SSE[kp]          = temp$SSE
    MSEfuzzy$MSE[kp]          = temp$MSE
    MSEfuzzy$SimResultspieces[[kp]] = temp$SimResults
    MSEfuzzy$model            = temp$model
    MSEfuzzy$nDataP[kp]       = temp$nDataP
    MSEfuzzy$SimResults       = rbind(MSEfuzzy$SimResults,temp$SimResults)

  }

  return(list(model = model,
              MSEpiece        = MSEfuzzy$MSE,
              MSE             = sum(MSEfuzzy$SSE)/sum(MSEfuzzy$nDataP),
              NAFac           = NAFac,
              nDataP          = sum(MSEfuzzy$nDataP),
              nDataPpiece     = MSEfuzzy$nDataP,
              sizeFac         = sizeFac,
              SimResults      = MSEfuzzy$SimResults,
              SimResultspiece = MSEfuzzy$SimResultspieces,
              SSE             = sum(MSEfuzzy$SSE),
              SSEpieces       = MSEfuzzy$SSE))
}
