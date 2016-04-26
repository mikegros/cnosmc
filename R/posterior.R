posterior = function(cl,Bstring,Gstring,gCube,nCube,kCube,model,paramsList,indexList,datalst){

  # Last update: April 14, 2015  Biljana
  # introduced datalst parameter that swithches off or on dataList usage


  # Last update: October 30, 2015
  #================================================================#
  #   posterior distribution of the parameters g,n,k for a given sub-model
  #================================================================#

  #   Use this function to:
  #   1. Evaluate posterior distribution of g,n,k for a given sub-graph

  # CALLED BY: run_mcmc_one_link

  #==============================================================
  #                          REQUIRED INPUTS
  #==============================================================
  #Bstring    - a sequence of 0 (link-off) and 1(link -on) for the given sub-graph
  #gCube      - last sampled value g for the corresponding link
  #nCube      - last sampled value n for the corresponding link
  #kCube      - last sampled value k for the corresponding link
  #model      - subgraph - a list of matrices and vectors
  #paramsList - parameters of the sub-graph, a list of different parameters including the data
  #indexList  -
  #datalst    - swithches off or on dataList usage


  #===============================================================
  #                           OUTPUTS
  #===============================================================
  # evaluated posterior distribution

  #===============================================================
  #                          MODEL DETAILS
  #===============================================================


  #===============================================================

  prior_g       = Logpriorg(gCube)
  prior_k       = Logpriork(kCube)
  prior_n       = Logpriorn(nCube)
  prior_Gstring = LogpriorGstring(Gstring)

  if(any(c(prior_g,prior_k,prior_n) == -Inf)) return(-Inf)

  if (is.null(datalst)){
    lik = -(1/2)*getMSEfuzzy(cl,
                             Bstring    = Bstring,
                             Gstring    = Gstring,
                             gCube      = gCube,
                             nCube      = nCube,
                             kCube      = kCube,
                             model      = model,
                             paramsList = paramsList,
                             indexList  = indexList,
                             sizeFac    = 0,
                             NAFac      = 0,
                             verbose    = FALSE)$SSE/0.1^2
  }else{
    lik = -(1/2)*getMSEfuzzyDataList(cl,
                                     Bstring    = Bstring,
                                     Gstring    = Gstring,
                                     gCube      = gCube,
                                     nCube      = nCube,
                                     kCube      = kCube,
                                     model      = model,
                                     paramsList = paramsList,
                                     indexList  = indexList,
                                     sizeFac    = 0,
                                     NAFac      = 0,
                                     verbose    = FALSE)$SSE/0.1^2
  }

  return(lik + prior_g + prior_k + prior_n + prior_Gstring)
}
