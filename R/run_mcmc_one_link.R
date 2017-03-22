run_mcmc_one_link <- function(cl,
                              Bstring,
                              Gstring,
                              p_link,
                              gCube,
                              nCube,
                              kCube,
                              sigsq,
                              model,
                              paramsList,
                              indexList,
                              simList,
                              cube_inds,
                              jump_size,
                              index){

  #Last update:   April 22 2016 Biljana
  #               added a datalst switch that switches on or off the dataList usage
  #               Changed the update of Gstring, gCube, kCube, nCube separately


  # Last update: October 30, 2015
  #================================================================#
  #   MH step to sample g,n,k for only one active link
  #================================================================#

  #   Use this function to:
  #   1. Sample g,n and k for one active link

  # CALLED BY: wrapper_to_sample_all_links

  #==============================================================
  #                          REQUIRED INPUTS
  #==============================================================
  #Bstring    - a sequence of 0 (link-off) and 1(link -on) for the given sub-graph
  #gCube      - last sampled value g for the corresponding link
  #nCube      - last sampled value n for the corresponding link
  #kCube      - last sampled value k for the corresponding link
  #jump_size          - proposal step, could be tuned later
  #model      - the full model - a list of matrices and vectors
  #paramsList - parameters of the full model, a list of different parameters including the data
  #indexList  -
  #index      - current index of the link to sample
  #datalst    - swithches off or on dataList usage
  #===============================================================
  #                           OUTPUTS
  #===============================================================
  # list of sampled values:
  #gCube      - new sampled value g for the corresponding link
  #nCube      - new sampled value g for the corresponding link
  #kCube      - new sampled value k for the corresponding link
  #jump_size  - tuned value, a vector  (jump_size_g,jump_size_n,jump_size_k)


  #===============================================================
  #                          MODEL DETAILS
  #===============================================================


  #===============================================================

  #params g,k,n

  accepted <- data.frame( g=0, n=0, k=0 )

  #sample Gstring separately

  Gstring_1  <- Gstring
  Gstring_0  <- Gstring
  n_params   <- length(gCube)

  Gstring_1[index] <- 1
  Gstring_0[index] <- 0

  like_Gstring_1    <- (-1/2)*sum(getMSEFuzzy(cl=cl,
                                          Bstring    = Bstring,
                                          Gstring    = Gstring_1,
                                          gCube      = gCube,
                                          nCube      = nCube,
                                          kCube      = kCube,
                                          sigsq      = sigsq,
                                          model      = model,
                                          paramsList = paramsList,
                                          indexList  = indexList,
                                          simList    = simList,
                                          cube_inds  = cube_inds,
                                          sizeFac    = 0,NAFac=0,verbose = FALSE)$SSEvectorScaled/sigsq)
  like_Gstring_0    <- (-1/2)*sum(getMSEFuzzy(cl=cl,
                                          Bstring    = Bstring,
                                          Gstring    = Gstring_0,
                                          gCube      = gCube,
                                          nCube      = nCube,
                                          kCube      = kCube,
                                          sigsq      = sigsq,
                                          model      = model,
                                          paramsList = paramsList,
                                          indexList  = indexList,
                                          simList    = simList,
                                          cube_inds  = cube_inds,
                                          sizeFac    = 0,NAFac=0,verbose = FALSE)$SSEvectorScaled/sigsq)
  tmp <- max(c(like_Gstring_0,like_Gstring_1))
  cond_p <- p_link*exp(like_Gstring_1-tmp)/(p_link*exp(like_Gstring_1-tmp)+(1-p_link)*exp(like_Gstring_0-tmp))
  Gstring[index] <- rbinom(1,1,cond_p)

  current_post <- posterior(cl = cl,
                            Bstring    = Bstring,
                            Gstring    = Gstring,
                            gCube      = gCube,
                            nCube      = nCube,
                            kCube      = kCube,
                            sigsq      = sigsq,
                            p_link     = p_link,
                            model      = model,
                            paramsList = paramsList,
                            indexList  = indexList,
                            simList    = simList,
                            cube_inds  = cube_inds)

  # sample g
  gCube_prop        <- gCube
  gCube_prop[index] <- rnorm(1, gCube[index], jump_size$g[index])

  prop_post <- posterior(cl = cl,
                         Bstring    = Bstring,
                         Gstring    = Gstring,
                         gCube      = gCube_prop,
                         nCube      = nCube,
                         kCube      = kCube,
                         sigsq      = sigsq,
                         p_link     = p_link,
                         model      = model,
                         paramsList = paramsList,
                         indexList  = indexList,
                         simList    = simList,
                         cube_inds  = cube_inds)
  alpha_g   <- prop_post - current_post

  if (all(!is.na(alpha_g) , runif(1) < exp(alpha_g))){
    gCube        <- gCube_prop
    current_post <- prop_post
    accepted$g   <- accepted$g + 1
  }

  # sample k
  kCube_prop        <- kCube
  kCube_prop[index] <- rnorm(1, kCube[index], jump_size$k[index])

  prop_post <- posterior(cl = cl,
                         Bstring    = Bstring,
                         Gstring    = Gstring,
                         gCube      = gCube,
                         nCube      = nCube,
                         kCube      = kCube_prop,
                         sigsq      = sigsq,
                         p_link     = p_link,
                         model      = model,
                         paramsList = paramsList,
                         indexList  = indexList,
                         simList    = simList,
                         cube_inds  = cube_inds)

  alpha_k   <- prop_post - current_post

  if (all(!is.na(alpha_k) , runif(1) < exp(alpha_k))){
    kCube        <- kCube_prop
    current_post <- prop_post
    accepted$k   <- accepted$k + 1
  }

  # sample n
  nCube_prop        <- nCube
  nCube_prop[index] <- rlnorm(1, log(nCube[index]), jump_size$n[index])

  prop_post <- posterior(cl = cl,
                         Bstring    = Bstring,
                         Gstring    = Gstring,
                         gCube      = gCube,
                         nCube      = nCube_prop,
                         kCube      = kCube,
                         sigsq      = sigsq,
                         p_link     = p_link,
                         model      = model,
                         paramsList = paramsList,
                         indexList  = indexList,
                         simList    = simList,
                         cube_inds  = cube_inds)
  alpha_n   <- prop_post - current_post +
               dlnorm(nCube[index], log(nCube_prop[index]), jump_size$n[index], log=T) -
               dlnorm(nCube_prop[index], log(nCube[index]), jump_size$n[index], log=T)

  if (all(!is.na(alpha_n) , runif(1) < exp(alpha_n))){
    nCube        <- nCube_prop
    current_post <- prop_post
    accepted$n   <- accepted$n + 1
  }

  list(gCube = gCube, kCube = kCube, nCube = nCube, Gstring = Gstring, post = current_post, accepted = accepted)
}
