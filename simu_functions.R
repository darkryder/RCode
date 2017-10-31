require(MASS)

# data = the reference data.frame or matrix of count data
# M = the number of samples
# N = the number of genes (by default the same number than the reference matrix)
# PI0 = the proportion of genes simulated under the null hypothesis
# H1 = the H1 parameter (tau, in the manuscript)
# HCG_FILTER = have high count genes to be filtered ? 
# HCG_TH = the high cont genes threshold 
# BALANCED = are the simulated data balanced (differentially expressed genes in both under and over expression)
# EQLIBSIZE = does the samples have the same lib.size. 
# DENSFUN = the density function used to simulated. By default "Poisson". 
getSimulatedRNASeqData = function(data = data, M = 20, N = NULL, PI0 = 0.5, H1 = 0.2, HCG_FILTER = FALSE, HCG_TH = 0.0001, BALANCED = TRUE, EQLIBSIZE = TRUE, DENSFUN = "Poisson"){
  SIMU = NULL
  
  #! Data
  if (!is.null(data)){
    rc <- rowSums(data)
    data <- data[rc > 0,]
 
    if (HCG_FILTER){
      POISSON_PARAM = NULL
      for (iii in 1:nrow(data)){
        x = as.numeric(data[iii, ])
        fit = fitdistr(x, densfun = "Poisson")
        POISSON_PARAM[iii] = fit$"estimate"[1]
      }
      aa = which(POISSON_PARAM/sum(POISSON_PARAM) >= HCG_TH)
      if (length(aa) > 0)
        data = data[-aa, ]
    }
    if (is.null(N))
    	rows = 1:nrow(data)
    if (!is.null(N)){
    	rows = round(runif(N, 1, nrow(data)))
    	aa = which(rows == 0)
    	if (length(aa) > 0)
      		rows[aa] = 1
    	data = data[rows, ]
    }
  }
  
  #! Simulate H0 data
  N0 = round(nrow(data)*PI0) # number under H0
  aa0 = 1:N0 # indice of H0 genes
  for (iii in aa0)
      SIMU = rbind(SIMU, getSimulatedLine_Poisson(as.numeric(data[iii, ]), M))
  
  #! Simulate H1 data
  if (PI0 < 1){
    aa1 = (N0+1):nrow(data)
     for (iii in aa1){
      sign = 1
      if (BALANCED == TRUE)
        sign = sign(runif(1, -1, 1))
      SIMU = rbind(SIMU, getSimulatedLine_Poisson(as.numeric(data[iii, ]), M, h1Param = H1, sign = sign))
     }
  }
  
  SIMU = as.data.frame(SIMU)
  names(SIMU) = paste("M", 1:ncol(SIMU), sep = "")
  rownames(SIMU) = paste(1:length(rows), rownames(data)[rows], sep = "__")
      
  if (EQLIBSIZE == FALSE){
    for (zzz in 2:ncol(SIMU))
      SIMU[, zzz] = round(SIMU[, zzz]*(abs(rnorm(1, 1, 0.8))+1))
  }
  
  res = list(Data = SIMU, Group = factor(c(rep("0", round(M/2)), rep("1", round(M/2)))))
  return(res)
}






getSimulatedLine_Poisson = function(x, M, h1Param = 0, sign = 1){
  fit = fitdistr(x, densfun = "Poisson")
  ok = "no"
  while(ok == "no"){
    #xsim = rpois(N, lambda = fit$"estimate"[1])
    xsim = c(
            rpois(round(M/2), lambda = fit$"estimate"[1]),
            rpois(round(M/2), lambda = fit$"estimate"[1]+fit$"estimate"[1]*h1Param*sign)
            )
    if (sum(xsim) > 0)
      ok = "yes"
  }
  return(xsim)
}


