wrapper_to_sample_all_links = function(cl,
                                       n_mh,
                                       Bstring,
                                       Gstring,
                                       p_link,
                                       gCube,
                                       nCube,
                                       kCube,
                                       sigsq,
                                       inhib_inds,
                                       model,
                                       paramsList,
                                       indexList,
                                       jump_size){

  #Last update:  April 14 2016 Biljana
  #              added a datalst switch that switches on or off the dataList usage

  # Last update: October 30, 2015

  # Last update: October 30, 2015
  #================================================================#
  #   wrapper_to_sample_all_links to sample g,n,k for all active links in one sub-graph
  #================================================================#

  #   Use this function to:
  #   1. Sample g,n and k for all active links in one sub-graph

  # CALLED BY: NA

  #==============================================================
  #                          REQUIRED INPUTS
  #==============================================================
  #Bstring    - a sequence of 0 (link-off) and 1(link -on) for the given sub-graph
  #gCube      - last sampled g for the current iteration graph
  #             a matrix(nrow=length(Bstring), ncol=2)
  #nCube      - last sampled g for the current iteration graph
  #             a matrix(nrow=length(Bstring), ncol=2)
  #kCube      - last sampled g for the current iteration graph
  #             a matrix(nrow=length(Bstring), ncol=2)
  #jump_size          - proposal step, could be tuned later
  #model      - the full model - a list of matrices and vectors
  #paramsList - parameters of the full model, a list of different parameters including the data
  #indexList  -
  #accs       - number of accepted proposals, a vector (acc_g,acc_n,acc_k)
  #datalst    - swithches off or on dataList usage
  #===============================================================
  #                           OUTPUTS
  #===============================================================
  # list of sampled values:
  #gCube      - same as the gCube input matrix, updated rows corresponding to the active links from BString
  #nCube      - same as the nCube input matrix, updated rows corresponding to the active links from BString
  #kCube      - same as the kCube input matrix, updated rows corresponding to the active links from BString
  #jump_size          - tuned value, a vector  (jump_size_g,jump_size_n,jump_size_k)


  #===============================================================
  #                          MODEL DETAILS
  #===============================================================


  #===============================================================

  inds <- which(Bstring==1)

  if (n_mh == 0) return(list(gCube = gCube, nCube = nCube, kCube = kCube, Gstring = Gstring, sigsq = sigsq))

  #sample the active links only

  for (j in 1:n_mh){
    fuzzy_out <- getMSEFuzzy(cl,
                             Bstring    = Bstring,
                             Gstring    = Gstring,
                             gCube      = gCube,
                             nCube      = nCube,
                             kCube      = kCube,
                             inhib_inds = inhib_inds,
                             model      = model,
                             paramsList = paramsList,
                             indexList  = indexList,
                             sizeFac    = 0,
                             NAFac      = 0,
                             verbose    = FALSE)
    NN    <- fuzzy_out$nDataP
    SSE   <- fuzzy_out$SSE
    sigsq <- 1/rgamma(1,1.25+NN/2,10^-5+SSE/2)

    for (ind in inds){

      SampleOneLink = run_mcmc_one_link(cl=cl,
                                        Bstring    = Bstring,
                                        Gstring    = Gstring,
                                        p_link     = p_link,
                                        gCube      = gCube,
                                        nCube      = nCube,
                                        kCube      = kCube,
                                        sigsq      = sigsq,
                                        inhib_inds = inhib_inds,
                                        model      = model,
                                        paramsList = paramsList,
                                        indexList  = indexList,
                                        jump_size  = jump_size,
                                        index      = ind)

      gCube   = SampleOneLink$gCube
      nCube   = SampleOneLink$nCube
      kCube   = SampleOneLink$kCube
      Gstring = SampleOneLink$Gstring
      post    = SampleOneLink$post
    }
  }

  list(gCube = gCube, nCube = nCube, kCube = kCube, Gstring = Gstring, sigsq = sigsq, post = post)
}
