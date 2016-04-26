########################################################
##############  SplitDataByInhibitor
########################################################

SplitDataByInhibitor = function(data){

  ########
  # Last update: April 11, 2016, Dave Campbell
  #
  #
  # Use this function to:
  # Split the data into a list where each component of the list
  # uses a different setting of the inhibitor
  #
  #
  # CALLED BY:
  #
  #
  #==============================================================
  #                         REQUIRED INPUTS
  #==============================================================
  # data      = the library(CNORfuzzy) data file
  #
  #
  #===============================================================%
  #                           OUTPUTS                             %
  #===============================================================%
  #
  # a list of valid data files where each takes only a single inhibitor setting
  #
  #===============================================================%
  #                          MODEL DETAILS                        %
  #===============================================================%
  # COMMENTS:
  #===============================================================%
  inhibs   = data$valueInhibitors
  siga     = data$valueSignals
  stim     = data$valueStimuli
  cues     = data$valueCues
  dataList = list()

  #Extract the data corresponding to the Inh1=0 and Inh2=0
  dataList[[1]]                      = data
  noInhibIndex                       = which(apply(data$valueInhibitors,1,sum)==0)
  dataList[[1]]$valueInhibitors      = data$valueInhibitors[noInhibIndex,]
  dataList[[1]]$valueSignals[[1]]    = data$valueSignals[[1]][noInhibIndex,]
  dataList[[1]]$valueSignals[[2]]    = data$valueSignals[[2]][noInhibIndex,]
  dataList[[1]]$valueStimuli         = data$valueStimuli[noInhibIndex,]
  dataList[[1]]$valueCues            = data$valueCues[noInhibIndex,]

  #Extract the data corresponding to the (Inh1=1,Inh2=0) and (Inh1=0,Inh2=1)
  for(InhibIndex in 1:dim(inhibs)[2]){
    dataList[[InhibIndex+1]] = data
    InhibOnIndex = which(inhibs[,InhibIndex]==1)
    dataList[[InhibIndex+1]]$valueInhibitors      = data$valueInhibitors[InhibOnIndex,]
    dataList[[InhibIndex+1]]$valueSignals[[1]]    = data$valueSignals[[1]][InhibOnIndex,]
    dataList[[InhibIndex+1]]$valueSignals[[2]]    = data$valueSignals[[2]][InhibOnIndex,]
    dataList[[InhibIndex+1]]$valueStimuli         = data$valueStimuli[InhibOnIndex,]
    dataList[[InhibIndex+1]]$valueCues            = data$valueCues[InhibOnIndex,]
  }
  return(dataList)
}
