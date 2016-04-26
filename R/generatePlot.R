generatePlot = function(tbld,model,smc_samples,data,datalst){
  #Last update:   April 14 2016 Biljana
  #              added a datalst switch that switches on or off the dataList usage


  #==============================================================
  #Use this function to: plot the model with highlights for the best model
  #                      and probabilities for each of the links
  #CALLED BY:  after SMC is run and after the best model is selected
  #             by SelectModel
  #==============================================================
  #                         REQUIRED INPUTS
  #==============================================================
  #tbld        = output from SelectModel i.e.,
  #              table with probabilities for each model
  #smc_samples = sampled parameters from SMC
  #model       = the model
  #data        = the data
  #datalst    - swithches off or on dataList usage
  #===============================================================%
  #                           OUTPUTS                             %
  #===============================================================%
  # saved plot as  model.png file
  #===============================================================%
  #                          MODEL DETAILS                        %
  #===============================================================%
  # COMMENTS:  requires two files which are modified plotModel function
  #           which allows passig of the edge labels
  # the file names are:plotModel_helpfns.R and plotModel_fn.R
  #===============================================================%
  #

  n_params=dim(smc_samples$gCube)[2]
  if (is.null(datalst)){
    #extract the model with highest probability (the most visited model)
    bitString=tbld[which(tbld[,ncol(tbld)]==max(tbld[,ncol(tbld)])),1:(ncol(tbld)-1)]
    #bitString has to be full size including the interactions
    #so we need to add interactions in bitString
    #make color shades based on the links
    bitString=unlist(c(bitString,rep(1,length(model$reacID)-length(bitString))))
    #calculate probability of each sampled link from Gstring
    tbl=as.data.frame(smc_samples$Gstring)
    edgeLabels=apply(tbl,2,mean)
    names(edgeLabels)=colnames(tbld)[1:n_params]
    names(edgeLabels)=unlist(lapply(strsplit(names(edgeLabels),'='),function(x) {paste(x,collapse='~')} ))
    names(edgeLabels)=sub('!', '', names(edgeLabels))
    #run first the refular plotModel function to get the structure of edge labels
    pl=plotModel(model,data,bString=bitString,show=F)
    pl_label=unlist(pl$edgeAttrs$label)
    #update the structure of edge labels with names of the main effects from our table
    pl_label[which(names(pl_label) %in% names(edgeLabels))]=edgeLabels[which(names(edgeLabels) %in% names(pl_label))]

    #pass all the parameters including updated edge labels to the modified plotModel function
    plotModel_modified(model,data,bString=bitString,graphvizParams=list(edgelabels=pl_label),output='PNG')
    #find code for createEdgeAttrs
  }else{


    for (i in (1:length(data))){
      tbl=get(load(paste('ModelSelect',i,'.RData',sep='')))

      indMax=which(tbl[,ncol(tbl)]==max(tbl[,ncol(tbl)]))
      indMax=sample(indMax,1)
      bitString=tbl[indMax,1:(ncol(tbl)-1)]
      bitString=unlist(c(bitString,rep(1,length(model$reacID)-length(bitString))))

      tblEdge=as.data.frame(smc_samples$Gstring[,(1:n_params) + n_params*(i-1)])
      edgeLabels=round(apply(tblEdge,2,mean),2)
      names(edgeLabels)=colnames(tbl)[1:n_params]
      names(edgeLabels)=unlist(lapply(strsplit(names(edgeLabels),'='),function(x) {paste(x,collapse='~')} ))
      names(edgeLabels)=sub('!', '', names(edgeLabels))
      #run first the refular plotModel function to get the structure of edge labels
      pl=plotModel(model,data[[i]],bString=bitString,show=F)
      pl_label=unlist(pl$edgeAttrs$label)
      #update the structure of edge labels with names of the main effects from our table
      pl_label[which(names(pl_label) %in% names(edgeLabels))]=edgeLabels[which(names(edgeLabels) %in% names(pl_label))]

      #pass all the parameters including updated edge labels to the modified plotModel function
      plotModel_modified(model,data[[i]],bString=bitString,graphvizParams=list(edgelabels=pl_label),output='PNG',filename = paste('Plot_model_inh',i,sep=''))
      #find code for createEdgeAttrs
    }
  }
}
