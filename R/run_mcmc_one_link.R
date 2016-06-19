run_mcmc_one_link <- function(cl,
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
                              sigma,
                              jump_size=rep(0.5,3),
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
  #sample Gstring separately

  if (!is.null(inhib_inds)){
    Gstring_1  <- Gstring
    Gstring_0  <- Gstring
    n_params   <- length(gCube)
    n_models   <- ncol(inhib_inds)

    for (i in (1:n_models)){

      Gstring_1[index+(i-1)*n_params] <- 1
      Gstring_0[index+(i-1)*n_params] <- 0

      like_Gstring_1    <- (-1/2)*getMSEFuzzy(cl=cl,
                                              Bstring    = Bstring,
                                              Gstring    = Gstring_1,
                                              gCube      = gCube,
                                              nCube      = nCube,
                                              kCube      = kCube,
                                              inhib_inds = inhib_inds,
                                              model      = model,
                                              paramsList = paramsList,
                                              indexList  = indexList,
                                              sizeFac    = 0,NAFac=0,verbose = FALSE)$SSE/sigma^2
      like_Gstring_0    <- (-1/2)*getMSEFuzzy(cl=cl,
                                              Bstring    = Bstring,
                                              Gstring    = Gstring_0,
                                              gCube      = gCube,
                                              nCube      = nCube,
                                              kCube      = kCube,
                                              inhib_inds = inhib_inds,
                                              model      = model,
                                              paramsList = paramsList,
                                              indexList  = indexList,
                                              sizeFac    = 0,NAFac=0,verbose = FALSE)$SSE/sigma^2

      tmp <- max(c(like_Gstring_0,like_Gstring_1))
      cond_p <- p_link*exp(like_Gstring_1-tmp)/(p_link*exp(like_Gstring_1-tmp)+(1-p_link)*exp(like_Gstring_0-tmp))
      Gstring[index+(i-1)*n_params] <- rbinom(1,1,cond_p)
    }
  }else{
    Gstring_1  <- Gstring
    Gstring_0  <- Gstring
    n_params   <- length(gCube)
    n_models   <- length(inhib_inds)


    Gstring_1[index] <- 1
    Gstring_0[index] <- 0

    like_Gstring_1    <- (-1/2)*getMSEFuzzy(cl=cl,
                                            Bstring    = Bstring,
                                            Gstring    = Gstring_1,
                                            gCube      = gCube,
                                            nCube      = nCube,
                                            kCube      = kCube,
                                            inhib_inds = inhib_inds,
                                            model      = model,
                                            paramsList = paramsList,
                                            indexList  = indexList,
                                            sizeFac    = 0,NAFac=0,verbose = FALSE)$SSE/sigma^2
    like_Gstring_0    <- (-1/2)*getMSEFuzzy(cl=cl,
                                            Bstring    = Bstring,
                                            Gstring    = Gstring_0,
                                            gCube      = gCube,
                                            nCube      = nCube,
                                            kCube      = kCube,
                                            inhib_inds = inhib_inds,
                                            model      = model,
                                            paramsList = paramsList,
                                            indexList  = indexList,
                                            sizeFac    = 0,NAFac=0,verbose = FALSE)$SSE/sigma^2
    tmp <- max(c(like_Gstring_0,like_Gstring_1))
    cond_p <- p_link*exp(like_Gstring_1-tmp)/(p_link*exp(like_Gstring_1-tmp)+(1-p_link)*exp(like_Gstring_0-tmp))
    Gstring[index] <- rbinom(1,1,cond_p)
  }
  current_post <- posterior(cl,Bstring,Gstring,gCube,nCube,kCube,p_link,inhib_inds,model,paramsList,indexList,sigma)

  # sample g
  gCube_prop        <- gCube
  gCube_prop[index] <- rnorm(1, gCube[index], jump_size[1])

  prop_post <- posterior(cl,Bstring,Gstring,gCube_prop,nCube,kCube,p_link,inhib_inds,model,paramsList,indexList,sigma)
  alpha_g   <- prop_post - current_post

  if (all(!is.na(alpha_g) , runif(1) < exp(alpha_g))){
    gCube        <- gCube_prop
    current_post <- prop_post
  }

  # sample k
  kCube_prop        <- kCube
  kCube_prop[index] <- rnorm(1, kCube[index], jump_size[2])

  prop_post <- posterior(cl,Bstring,Gstring,gCube,nCube,kCube_prop,p_link,inhib_inds,model,paramsList,indexList,sigma)
  alpha_k   <- prop_post - current_post

  if (all(!is.na(alpha_k) , runif(1) < exp(alpha_k))){
    kCube <- kCube_prop
    current_post <- prop_post
  }

  # sample n
  nCube_prop <- nCube
  delta      <- rnorm(1,0,jump_size[3])

  nCube_prop[index] <- exp(log(nCube[index]) + delta)

  prop_post <- posterior(cl,Bstring,Gstring,gCube,nCube_prop,kCube,p_link,inhib_inds,model,paramsList,indexList,sigma)
  alpha_n   <- prop_post - current_post +
               dlnorm(nCube, log(nCube_prop), jump_size[3], log=T) -
               dlnorm(nCube_prop, log(nCube), jump_size[3], log=T)

  if (all(!is.na(alpha_n) , runif(1) < exp(alpha_n))){
    nCube <- nCube_prop
    current_post <- prop_post
  }

  list(gCube=gCube,kCube=kCube,nCube=nCube,Gstring=Gstring)
}
