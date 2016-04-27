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
    acceptGstring <- 0
    Gstring_prop  <- Gstring
    n_params      <- length(gCube)
    n_models      <- length(inhib_inds)

    for (i in (1:n_models)){

      Gstring_prop[index+(i-1)*n_params] <- rbinom(1,1,Gstring[index+(i-1)*n_params]*p_link+(1-Gstring[index+(i-1)*n_params])*(1-p_link))

      alpha_Gstring = posterior(cl,Bstring,Gstring_prop,gCube, nCube, kCube,inhib_inds,model,paramsList,indexList) -
                      posterior(cl,Bstring,Gstring,     gCube, nCube, kCube,inhib_inds,model,paramsList,indexList)

      if (all(!is.na(alpha_Gstring) , runif(1) < exp(alpha_Gstring))){
        Gstring       <- Gstring_prop
        acceptGstring <- acceptGstring + 1
      }
    }
  }else{
    acceptGstring       <- 0
    Gstring_prop        <- Gstring
    Gstring_prop[index] <- rbinom(1,1,Gstring[index]*p_link+(1-Gstring[index])*(1-p_link))

    alpha_Gstring = posterior(cl,Bstring,Gstring_prop,gCube,nCube,kCube,inhib_inds,model,paramsList,indexList) -
                    posterior(cl,Bstring,Gstring,     gCube,nCube,kCube,inhib_inds,model,paramsList,indexList)

    if (all(!is.na(alpha_Gstring) , runif(1) < exp(alpha_Gstring))){
      Gstring       <- Gstring_prop
      acceptGstring <- acceptGstring + 1
    }
  }

  # sample g
  gCube_prop        <- gCube
  gCube_prop[index] <- rnorm(1, gCube[index], jump_size[1])

  alpha_g <- posterior(cl,Bstring,Gstring,gCube_prop,nCube,kCube,inhib_inds,model,paramsList,indexList) -
             posterior(cl,Bstring,Gstring,gCube,     nCube,kCube,inhib_inds,model,paramsList,indexList)


  if (all(!is.na(alpha_g) , runif(1) < exp(alpha_g))){
    gCube <- gCube_prop
  }

  # sample k
  kCube_prop        <- kCube
  kCube_prop[index] <- rnorm(1, kCube[index], jump_size[2])

  alpha_k <- posterior(cl,Bstring,Gstring,gCube,nCube,kCube_prop,inhib_inds,model,paramsList,indexList) -
             posterior(cl,Bstring,Gstring,gCube,nCube,kCube,     inhib_inds,model,paramsList,indexList)

  if (all(!is.na(alpha_k) , runif(1) < exp(alpha_k))){
    kCube <- kCube_prop
  }

  # sample n
  nCube_prop <- nCube
  delta      <- rnorm(1,0,jump_size[3])

  nCube_prop[index] <- exp(log(nCube[index]) + delta)

  alpha_n <- posterior(cl,Bstring,Gstring,gCube,nCube_prop,kCube,inhib_inds,model,paramsList,indexList) -
             posterior(cl,Bstring,Gstring,gCube,nCube,     kCube,inhib_inds,model,paramsList,indexList) +
             dlnorm(nCube, log(nCube_prop), jump_size[3], log=T) -
             dlnorm(nCube_prop, log(nCube), jump_size[3], log=T)

  if (all(!is.na(alpha_n) , runif(1) < exp(alpha_n))){
    nCube <- nCube_prop
  }

  list(gCube=gCube,kCube=kCube,nCube=nCube,Gstring=Gstring,acceptGstring=acceptGstring)
}
