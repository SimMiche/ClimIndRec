#!/usr/bin/env Rscript

library(stringr)
read_params = function(filepath) {
  params=list()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    spl=stringr::str_split(line,'=')
    print(spl[[1]])
    if (length(spl[[1]])>1){
        params=c(params,list(spl[[1]][length(spl[[1]])]))
    }
  }
  close(con)
  return(params)
}

params=read_params('./params.txt')

source('./ClimIndRec1.1.r')
results=apply_rec(workdir=as.character(params[[1]]),path_db=as.character(params[[3]]),path_mode=as.character(params[[4]]),y1=as.numeric(params[[5]]),y2=as.numeric(params[[6]]),method=as.character(params[[8]]),R=as.numeric(params[[7]]),freq_calib=as.numeric(params[[9]]),tests=as.logical(params[[10]]),blockstyle_cv=as.logical(params[[13]]),conf=as.numeric(params[[11]]),K_cv=as.character(params[[14]]),seed=as.numeric(params[[15]]),trace=as.logical(params[[16]]),blockstyle_holdout=as.logical(params[[12]]))

setwd(params[[1]])
dir.create(params[[2]])
setwd(paste(params[[1]],params[[2]],sep='/'))

write.csv(results[[1]],paste(params[[2]],'_renormalized_reconstruction.csv',sep=''),quote=F,row.names=F)

write.csv(results[[2]],paste(params[[2]],'_original_reconstruction.csv',sep=''),quote=F,row.names=F)                      

write.csv(results[[3]],paste(params[[2]],'_correlation_scores.csv',sep=''),quote=F,row.names=F)                           
 
write.csv(results[[4]],paste(params[[2]],'_RMSE_scores.csv',sep=''),quote=F,row.names=F) 

write.csv(results[[5]],paste(params[[2]],'_nb_records_testing.csv',sep=''),quote=F,row.names=F)

write.csv(results[[6]],paste(params[[2]],'_name_proxies_testing.csv',sep=''),quote=F,row.names=F)

write.csv(results[[7]],paste(params[[2]],'_individual_reconstructions.csv',sep=''),quote=F,row.names=F)

write.csv(results[[8]],paste(params[[2]],'_NSCE_scores.csv',sep=''),quote=F,row.names=F)

write.csv(results[[9]],paste(params[[2]],'_testing_uncertainties.csv',sep=''),quote=F,row.names=F)

write.csv(results[[10]],paste(params[[2]],'_training_samples.csv',sep=''),quote=F,row.names=F)

write.csv(results[[11]],paste(params[[2]],'_testing_samples.csv',sep=''),quote=F,row.names=F)

write.csv(results[[12]],paste(params[[2]],'_testing_ShapiroWilk_pvalues_residuals.csv',sep=''),quote=F,row.names=F)

write.csv(results[[13]],paste(params[[2]],'_namas_proxies_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[14]],paste(params[[2]],'_NSCE_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[15]],paste(params[[2]],'_correlation_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[16]],paste(params[[2]],'_RMSE_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[17]],paste(params[[2]],'_uncertainties_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[18]],paste(params[[2]],'_ShapiroWilk_pvalues_residuals_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[19]],paste(params[[2]],'_nb_records_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[20]],paste(params[[2]],'_name_proxies_final_model.csv',sep=''),quote=F,row.names=F)

write.csv(results[[21]],paste(params[[2]],'_RE_scores.csv',sep=''),quote=F,row.names=F)
