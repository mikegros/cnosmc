selectModel <-function(smc_samples,model,n_params,datalst){
  #Last update:   April 14 2016 Biljana
  #              added a datalst switch that switches on or off the dataList usage


  #==============================================================
  #Use this function to: get the table with probabilities for each model
  #CALLED BY:  after SMC is run
  #==============================================================
  #                         REQUIRED INPUTS
  #==============================================================
  #smc_samples = sampled parameters from SMC
  #model       = the model
  ##datalst    - swithches off or on dataList usage
  #
  #===============================================================%
  #                           OUTPUTS                             %
  #===============================================================%
  # table with distinct visited models with their associated probabilities
  # the table object is returned and also saved as ModelSelect.RData
  #===============================================================%
  #                          MODEL DETAILS                        %
  #===============================================================%
  # COMMENTS:  require(sqldf)
  #===============================================================%
  #
  if (is.null(datalst)){
    tbl=as.data.frame(smc_samples$Gstring)
    cnames=sapply(1:ncol(tbl), function(x) paste('c',x,sep=''))
    cnames = paste(cnames, collapse=",")[1]
    colnames(tbl)=strsplit(cnames,split="," )[[1]]

    tbl1=sqldf(paste('select *,count(*) from tbl group by ', cnames,sep=''))
    tbl1[,ncol(tbl1)]=tbl1[,ncol(tbl1)]/sum(tbl1[,ncol(tbl1)])

    colnames(tbl1)=c(colnames(model$interMat)[1:n_params],'Prob.Model')
    save(tbl1,file='ModelSelect.RData')
  }else{
    tbl1=list()
    for (i in (1:(dim(smc_samples$Gstring)[2]/n_params))){
      tbl                         = as.data.frame(smc_samples$Gstring[,(1:n_params) + n_params*(i-1)])
      cnames                      = sapply(1:ncol(tbl), function(x) paste('c',x,sep=''))
      cnames                      = paste(cnames, collapse=",")[1]
      colnames(tbl)               = strsplit(cnames,split="," )[[1]]

      tbl1[[i]]                   = sqldf(paste('select *,count(*) from tbl group by ', cnames,sep=''))
      tbl1[[i]][,ncol(tbl1[[i]])] = tbl1[[i]][,ncol(tbl1[[i]])]/sum(tbl1[[i]][,ncol(tbl1[[i]])])

      colnames(tbl1[[i]])         = c(colnames(model$interMat)[1:n_params],'Prob.Model')
      tblSave                     = tbl1[[i]]
      save(tblSave,file=paste('ModelSelect',i,'.RData',sep=''))
    }
  }
  return(tbl1)
}
