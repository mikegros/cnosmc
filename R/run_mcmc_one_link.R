run_mcmc_one_link <- function(cl,
                              Bstring,
                              Gstring,
                              p_link,
                              gCube,
                              nCube,
                              kCube,
                              model,
                              paramsList,
                              indexList,
                              jump_size=rep(0.5,3),
                              index,
                              datalst){

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
  # dirak_delta =function(x,a)  {(1/(a*sqrt(pi)))*exp(-x^2/a^2)}

  if (!is.null(datalst)){
    acceptGstring                       = 0
    Gstring_prop                        = Gstring
    lgcube                              = length(gCube)
    for (i in (1:length(paramsList$dataList))){

      Gstring_prop[index+(i-1)*lgcube]    = rbinom(1,1,Gstring[index+(i-1)*lgcube]*p_link+(1-Gstring[index+(i-1)*lgcube])*(1-p_link))


      paramsListUpd                       = paramsList[[i]]
      paramsListUpd$data                  = paramsList$dataList[[i]]

      alpha_Gstring = posterior(cl,Bstring,Gstring_prop[(1:lgcube)+lgcube*(i-1)],gCube, nCube, kCube,model,paramsListUpd,indexList[[i]],datalst=NULL) -
        posterior(cl,Bstring,Gstring[(1:lgcube)+lgcube*(i-1)],     gCube, nCube, kCube,model,paramsListUpd,indexList[[i]],datalst=NULL)




      if (all(!is.na(alpha_Gstring) , runif(1) < exp(alpha_Gstring))){
        Gstring                            = Gstring_prop
        acceptGstring                      = acceptGstring+1
      }
      #
      #       #Gstring_prop=Gstring
    }
    # Gstring=Gstring_prop
  }else{
    Gstring_prop                        = Gstring
    Gstring_prop[index]                 = rbinom(1,1,Gstring[index]*p_link+(1-Gstring[index])*(1-p_link))


    alpha_Gstring = posterior(cl,Bstring,Gstring_prop[(1:lgcube)+lgcube*(i-1)],gCube, nCube, kCube,model,paramsListUpd,indexList[[i]],datalst=NULL) -
      posterior(cl,Bstring,Gstring[(1:lgcube)+lgcube*(i-1)],     gCube, nCube, kCube,model,paramsListUpd,indexList[[i]],datalst=NULL)


    if (all(!is.na(alpha_Gstring) , runif(1) < exp(alpha_Gstring))){
      Gstring                            = Gstring_prop
      acceptGstring                      = acceptGstring+1
    }


  }

  #sample g
  gCube_prop                             = gCube
  gCube_prop[index]                      = rnorm(1, gCube[index], jump_size[1])
  #   delta_g                                = rnorm(1,0,jump_size[1])
  #   gCube_prop[index]                      = exp(log(gCube[index]) + delta_g)

  alpha_g = posterior(cl,Bstring,Gstring,gCube_prop, nCube, kCube,model,paramsList,indexList,datalst) -
    posterior(cl,Bstring,Gstring,gCube,      nCube, kCube,model,paramsList,indexList,datalst)


  if (all(!is.na(alpha_g) , runif(1) < exp(alpha_g))){
    gCube   = gCube_prop
  }


  #sample k
  kCube_prop                             = kCube
  kCube_prop[index]                      = rnorm(1, kCube[index], jump_size[2])


  alpha_k = posterior(cl,Bstring,Gstring,gCube, nCube, kCube_prop,model,paramsList,indexList,datalst) -
    posterior(cl,Bstring,Gstring,gCube, nCube, kCube,model,paramsList,indexList,datalst)


  if (all(!is.na(alpha_k) , runif(1) < exp(alpha_k))){
    kCube   = kCube_prop
  }


  nCube_prop                             = nCube
  delta                                  = rnorm(1,0,jump_size[3])
  nCube_prop[index]                      = exp(log(nCube[index]) + delta)


  alpha_n = posterior(cl,Bstring,Gstring,gCube, nCube_prop, kCube,model,paramsList,indexList,datalst) -
    posterior(cl,Bstring,Gstring,gCube, nCube,      kCube,model,paramsList,indexList,datalst) +
    dlnorm(nCube, log(nCube_prop), jump_size[3], log=T) -
    dlnorm(nCube_prop, log(nCube), jump_size[3], log=T)



  if (all(!is.na(alpha_n) , runif(1) < exp(alpha_n))){
    nCube   = nCube_prop
  }


  #   if (!is.null(datalst)){
  #   alpha_n = posterior(cl,Bstring,Gstring,gCube_prop, nCube_prop, kCube_prop,model,paramsList,indexList,datalst) -
  #             posterior(cl,Bstring,Gstring,gCube,      nCube,      kCube,     model,paramsList,indexList,datalst) +
  #             dlnorm(nCube, log(nCube_prop), jump_size[3], log=T) -
  #             dlnorm(nCube_prop, log(nCube), jump_size[3], log=T)
  #
  #   accept = 0
  #
  #   if (all(!is.na(alpha_n) , runif(1) < exp(alpha_n))){
  #     gCube   = gCube_prop
  #     nCube   = nCube_prop
  #     kCube   = kCube_prop
  #     accept  = 1
  #   }
  #
  #
  #   }else{
  #     alpha_n = posterior(cl,Bstring,Gstring_prop,gCube_prop, nCube_prop, kCube_prop,model,paramsList,indexList,datalst=NULL) -
  #       posterior(cl,Bstring,Gstring,gCube,      nCube,      kCube,     model,paramsList,indexList,datalst=NULL) +
  #       dlnorm(nCube, log(nCube_prop), jump_size[3], log=T) -
  #       dlnorm(nCube_prop, log(nCube), jump_size[3], log=T)
  #
  #
  #
  #   accept = 0
  #
  #   if (all(!is.na(alpha_n) , runif(1) < exp(alpha_n))){
  #     gCube   = gCube_prop
  #     nCube   = nCube_prop
  #     kCube   = kCube_prop
  #     Gstring = Gstring_prop
  #     accept  = 1
  #   }
  #   }

  list(gCube=gCube,kCube=kCube,nCube=nCube,Gstring=Gstring,acceptGstring=acceptGstring)
}
