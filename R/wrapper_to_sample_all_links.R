wrapper_to_sample_all_links = function(cl,
                                       n_mh,
                                       Bstring,
                                       Gstring,
                                       p_link,
                                       gCube,
                                       nCube,
                                       kCube,
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

  #sample the active links only

  # MJG: Added variable for counting number of times that the jitter step is accepted for diagnostic
  #      purposes, though removed the print statements for now. When parallelized, be sure to set the
  #      outfile argument in the makeCluster command so that the print statements are returned.
  accepts <- rep(0,length(inds))
  for (j in 1:n_mh){
    for (ind in inds){

      SampleOneLink = run_mcmc_one_link(cl=cl,
                                        Bstring    = Bstring,
                                        Gstring    = Gstring,
                                        p_link     = p_link,
                                        gCube      = gCube,
                                        nCube      = nCube,
                                        kCube      = kCube,
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
      accepts[which(inds==ind)] <- accepts[which(inds==ind)] + SampleOneLink$accept
    }}

  list(gCube = gCube, nCube = nCube,kCube = kCube,Gstring=Gstring)
}
