PCO_acat = function(Z.mat,Sigma,p.method="TruncPCO",method = "davies"){
  ## avoid Z.mat is a row vector which won't work later
  ## Z.mat should be a column vector object, not a 1*K matrix 
  
  ##########
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  
  if(p.method != "PCO"){
    lambdas = lambdas[lambdas>0.1]
    eigen.vec = eigen.vec[, which(lambdas>0.1)]
  }
  
  K = length(lambdas)
  w.vec = 1/lambdas
  
  z.tmp = t(Z.mat)
  PC.vec = t(eigen.vec)%*%z.tmp
  
  
  ## 1: PCMinP
  PC.std = (PC.vec)^2 * w.vec
  PC.p = pchisq(PC.std,df=1,lower.tail = FALSE)
  p.min = apply(PC.p, 2, min)
  #p.PCMinP = 1- (1-p.min)^K
  p.PCMinP = -expm1(K*log1p(-p.min))
  cat("PCMinP done.", "\n")
  
  ## 2: PCFisher
  PC.Fisher.stat = -2*colSums(log(PC.p))
  p.PCFisher = pchisq(PC.Fisher.stat,df=2*K,lower.tail = FALSE)
  cat("PCFisher done.", "\n")
  
  ## 3: PCLC
  PCLC.stat = colSums(PC.vec * w.vec)
  PCLC.stat.std = PCLC.stat^2/sum(w.vec)
  p.PCLC = pchisq(PCLC.stat.std,df=1,lower.tail = FALSE)
  cat("PCLC done.", "\n")
  
  ## 4: WI
  T.WI = colSums((PC.vec)^2)
  
  if(method=="davies"){
    p.WI = sapply(T.WI, function(x) davies(x, lambdas)$Qq)
    if(!all(p.WI > 0)){
      ind.liumod = !(p.WI > 0)
      snp.liumod = names(p.WI[ind.liumod])
      p.WI[snp.liumod] = liumod(T.WI[snp.liumod], lambdas)
    }
  }else if(method=="liu"){
    p.WI = liu(T.WI,lambdas)
  }else if(method=="liumod"){
    p.WI = liumod(T.WI,lambdas)
  }else {stop("Please specify a valid method: davies, liu, liumod.\n")}
  cat("WI done.", "\n")
  
  ## 5: Wald
  T.Wald = colSums(PC.std)
  p.Wald = pchisq(T.Wald,df=K,lower.tail=FALSE)
  cat("Wald done.", "\n")
  
  
  ## 6: VC
  T.VC = colSums((PC.vec * w.vec)^2)
  wt.VC = sort(w.vec, decreasing = TRUE)
  
  if(method=="davies"){
    p.VC = sapply(T.VC, function(x) davies(x, wt.VC)$Qq)
    if(!all(p.VC > 0)){
      ind.liumod = !(p.VC > 0)
      snp.liumod = names(p.VC[ind.liumod])
      p.VC[snp.liumod] = liumod(T.VC[snp.liumod], wt.VC)
    }
  }else if(method=="liu"){
    p.VC = liu(T.VC, wt.VC)
  }else if(method=="liumod"){
    p.VC = liumod(T.VC, wt.VC)
  }else {stop("Please specify a valid method: davies, liu, liumod.\n")}
  cat("VC done.", "\n")
  
  
  ## PCO ACAT
  #T.acat = apply(rbind(p.PCMinP,p.PCFisher,p.PCLC, p.WI,p.Wald,p.VC), 2, 
  #	      function(x){mean(tan((0.5 - x) * pi))})
  #pval = 1/2 - atan(T.acat) / pi
  pval = ACAT(rbind(p.PCMinP,p.PCFisher,p.PCLC,p.WI,p.Wald,p.VC), is.check = F)
  cat("PCO ACAT done.", "\n")
  
  return(pval)
}
