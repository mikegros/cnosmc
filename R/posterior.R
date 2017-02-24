posterior = function(cl,
                     Bstring,
                     Gstring,
                     gCube,
                     nCube,
                     kCube,
                     sigsq,
                     p_link,
                     model,
                     paramsList,
                     indexList,
                     cube_inds){

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
  prior_Gstring = LogpriorGstring(Gstring,p_link)

  if(any(c(prior_g,prior_k,prior_n) == -Inf)) return(-Inf)

  get_mse <- getMSEFuzzy(cl,
                         Bstring    = Bstring,
                         Gstring    = Gstring,
                         gCube      = gCube,
                         nCube      = nCube,
                         kCube      = kCube,
                         sigsq      = sigsq,
                         model      = model,
                         paramsList = paramsList,
                         indexList  = indexList,
                         cube_inds  = cube_inds,
                         sizeFac    = 0,
                         NAFac      = 0,
                         verbose    = FALSE)

  lik = -(1/2)*sum(get_mse$SSEvectorScaled[get_mse$active_nodes]) - sum(get_mse$nDataP[get_mse$active_nodes]/2*log(sigsq[get_mse$active_nodes]))

  prior_sigsq   = Logpriorsigsq(sigsq[get_mse$active_nodes],1.25,10^-5)

  return(lik + prior_g + prior_k + prior_n + prior_sigsq + prior_Gstring)
}
