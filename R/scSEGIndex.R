#' scSEGIndex
#'
#' Calculate single-cell Stably Expressed Gene (scSEG) index
#'
#'
#' @author Shila Ghazanfar, Yingxin Lin, Pengyi Yang
#' @param exprsMat The log-transoformed single-cell data (assumed to be no batch effect and covered a wide range of cell types). A n by m matrix, where n is the number of genes and m is the number of cells
#' @param cell_type A vector indicating the cell type information for each cell in the gene expression matrix. If it is \code{NULL}, the function calculates the scSEG index without using F-statistics.
#' @param ncore Number of cores that are used in parallel
#'
#' @return A matrix of stably expressed features.
#' @importFrom foreach %dopar%
#' @export


# require(distr)
# require(scales)
# require(parallel)
# require(foreach)
# require(doSNOW)



# This is the main function for calculating stably expressed gene index
scSEGIndex <- function(exprsMat, cell_type = NULL, ncore = 1) {

  if(is.null(exprsMat)){
    stop("exprsMat is NULL.")
  }



  if(is.null(cell_type)){
    print("Calculating scSEG index without F-statistics ")
  }else{
    if(length(cell_type)!=ncol(exprsMat)){
      stop("length of cell type information is not equal to the number of column of exprsMat.")
    }
    cell_type <- as.factor(cell_type)
    print("Calculating scSEG index with F-statistics ")
  }

  ## Core feature 1: mixture modelling
  # Zero%, remove those that have more than 80% missing values before fitting the mixture model
  # rownames(exprsMat) <- toupper(rownames(exprsMat))
  z.all <- apply(exprsMat,1,function(x)sum(x==0))/dim(exprsMat)[2]
  del <- which(z.all > 0.8)
  if(length(del) > 0) {
    exprsMat_filt <- exprsMat[-del,]
    if(nrow(exprsMat_filt) < 100) {
      stop("Not enough gene pass the QC!")
    }
  }else{
    exprsMat_filt <- exprsMat
  }



  # fitting the model
  print("Fitting the mixture model...")
  paraMat <- make_para_gn_parallel(as.matrix(exprsMat_filt), ncore = ncore)
  r <- paraMat$rho
  s <- paraMat$sigma
  m <- paraMat$mu
  m.scaled <-(m-min(m))/(max(m)-min(m))
  # names(r) <- names(s) <- names(m) <- names(m.scaled) <- toupper(rownames(mat))
  names(r) <- names(s) <- names(m) <- names(m.scaled) <- rownames(exprsMat_filt)


  genes <- rownames(exprsMat_filt)
  z <- z.all[genes] * m.scaled[genes]

  ### Fitting ANOVA model for F-stats
  if(is.null(cell_type)){

    ### Combine all core features for creating scSEG index
    x1 = rank(r)/(length(r)+1)
    x2 = 1-rank(s)/(length(s)+1)
    x3 = 1-rank(z)/(length(z)+1)
    # x4 = 1-rank(f)/(length(f)+1)
    segIdx <- rowMeans(cbind(x1, x2, x3))

    resMat <- data.frame(segIdx = segIdx, rho = r,
                         sigma = s, mu = m,
                         mu.scaled = m.scaled, zero = z)

  }else{

    print("Fitting ANOVA model for F-stats...")
    aovStats <- apply(exprsMat_filt, 1, function(x){
      tryCatch(aov(as.numeric(x) ~ cell_type), error = function(e){NULL})
    })
    f <- log2(sapply(aovStats, function(x){
      tryCatch(summary(x)[[1]]$`F value`[1], error = function(e){NA})
    }))

    # f <- f.all[genes]

    ### Combine all core features for creating scSEG index
    x1 = rank(r)/(length(r)+1)
    x2 = 1-rank(s)/(length(s)+1)
    x3 = 1-rank(z)/(length(z)+1)
    x4 = 1-rank(f)/(length(f)+1)
    segIdx <- rowMeans(cbind(x1, x2, x3, x4))
    resMat <- data.frame(segIdx = segIdx, rho = r,
                         sigma = s, mu = m,
                         mu.scaled = m.scaled, zero = z, f_stats = f)

  }

  return(resMat)
}




# This function perform mixture model fitting using multiple CPUs
make_para_gn_parallel = function(exprsMat, ncore = 1) {

  cl <- parallel::makeCluster(ncore, type='SOCK')
  doSNOW::registerDoSNOW(cl)

  iter <- nrow(exprsMat)
  pb <- txtProgressBar(min = 1, max = iter,style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  # `%dopar%` <- foreach::`%dopar%`
  res <- foreach::foreach(i = iterators::icount(iter), .combine = rbind,
                 .options.snow = opts,.export = c("gammaNormMix","bic","aic","icl_bic"))%dopar%{
    setTxtProgressBar(pb,i)
    exprsMat_gene <- exprsMat[i,]
    gene_name <- rownames(exprsMat)[i]
    para <-gammaNormMix(exprsMat_gene,verbose=FALSE,plot=FALSE)

    if(sum(is.na(para))==0){
      dat_para <- c(mu=para$mu,sigma=para$sd,rho=para$rho)
    }else{
      dat_para <- c(NA,NA,NA)
    }
    dat_para <- matrix(dat_para,ncol=3)
    rownames(dat_para) <- gene_name
    return(dat_para)
  }
  close(pb)
  stopCluster(cl)
  colnames(res)<-c("mu","sigma","rho")
  res<-data.frame(res)
  return(res)
}


bic<-function(loglik,n,p){
  return(-2*loglik+p*log(n))
}


aic<-function(loglik,p){
  return(-2*loglik+2*p)
}#Tend to fit more component


icl_bic<-function(loglik,postprob,n,p){
  postprob<-postprob[postprob>0]
  EN=-sum(postprob*log(postprob))
  return(-2*loglik+2*EN+p*log(n))
}




gammaNormMix = function(data, thresh = 1e-7, maxiter = 10000,
                        removeZeroes = TRUE, plot = TRUE, hist=TRUE,hist_col="light cyan",
                        verbose = FALSE, forceExponential = FALSE,
                        calculateAreaDifference = FALSE,
                        minDataPoints = 5,onlyAddCurves=FALSE,
                        addContextData = FALSE, contextData = NULL) {


  #fitting a 2 component normal and gamma mixture model

  # add other data to fit the model with as well, but only return
  # the classification for those we're interested in

  if (addContextData) {
    nOriginal = length(data)
    data <- c(data,contextData)
  }

  # assume all values exactly zero already belong to the gamma comp
  # and remove them from the EM algorithm

  if (removeZeroes) {
    nonZeroInd = which(data>0)
    x = data[nonZeroInd]
  } else {
    x = data
  }

  if(length(x)<minDataPoints) {
    if (verbose) cat("Not enough data points to fit mixture model!")
    return(NA)
  }

  # initiate
  n = length(x)
  z = rbinom(n,1,0.5)
  z_iter = z
  mu = -100
  mu_iter = 10
  sig2 = -100
  sig2_iter = 0
  alpha = -100
  alpha_iter = 1
  beta = -100
  beta_iter = 1
  rho = -100
  rho_iter = 0.5
  niter = 0

  while(any(c(
    abs(mu - mu_iter)>thresh ,
    abs(sig2 - sig2_iter)>thresh ,
    abs(alpha - alpha_iter)>thresh ,
    abs(beta - beta_iter)>thresh ,
    abs(rho - rho_iter)>thresh
  )) & (niter<maxiter) ) {

    # save old parameters
    mu = mu_iter
    sig2 = sig2_iter
    alpha = alpha_iter
    beta = beta_iter
    rho = rho_iter
    if (forceExponential) alpha_iter = 1 # SHILA

    niter = niter + 1

    # M step
    mu_iter = sum(z_iter*x)/sum(z_iter)
    sig2_iter = sum(z_iter*(x-mu_iter)*(x-mu_iter))/sum(z_iter)
    if (sig2_iter<=0 | is.na(sig2_iter)) sig2_iter = 1e-11
    beta_iter = alpha_iter*sum(1-z_iter)/sum((1-z_iter)*x)
    if (beta_iter<=0 | is.na(beta_iter)) beta_iter = 3
    if (!forceExponential) {
      alpha_iter = distr::igamma(sum((log(beta_iter) + log(x))*(1-z_iter))/sum(1-z_iter))
    }
    if (alpha_iter>150 | is.na(alpha_iter)) alpha_iter=150
    rho_iter = sum(z_iter)/n


    # E step
    eta_iter =
      - 0.5*log(2*pi*sig2_iter) -
      ((x-mu_iter)*(x-mu_iter))/(2*sig2_iter) -
      alpha_iter*log(beta_iter) +
      log(gamma(alpha_iter)) -
      (alpha_iter-1)*log(x) +
      beta_iter*x +
      log(rho_iter/(1-rho_iter))
    z_iter = 1/(1+exp(-eta_iter))

    if (verbose) cat(niter,mu_iter,sqrt(sig2_iter),alpha_iter,beta_iter,rho_iter,"\n")
  }


  ll<-sum(log(rho_iter*dnorm(x,mu_iter,sqrt(sig2_iter)) +
                (1-rho_iter)*dgamma(x,shape=alpha_iter,rate=beta_iter)))


  xg <- seq(0,max(x)+1,length.out=300)
  c1g <- rho_iter*dnorm(xg,mu_iter,sqrt(sig2_iter))
  # c2g <- (1-rho_iter)*dexp(xg,lambda_iter)
  c2g <- (1-rho_iter)*dgamma(xg,shape=alpha_iter,rate=beta_iter)
  fg <- rho_iter*dnorm(xg,mu_iter,sqrt(sig2_iter)) +
    (1-rho_iter)*dgamma(xg,shape=alpha_iter,rate=beta_iter)
  # fg2 <- rho_true*dnorm(xg,mu_true,sqrt(sig2_true)) + (1-rho_true)*dexp(xg,lambda_true)
  if (plot) {
    if(hist){
      hist(x,probability=TRUE,col=hist_col,breaks=50,main=NA,xlab=NA,ylab="Density (zeroes removed)",ylim=c(0,0.6),xlim=c(0,20))
    }
    if (!onlyAddCurves) {
      lines(density(x,from=0),lty=2,lwd=2,col=alpha("darkgrey",0.6))
    }
    lines(xg,c1g,col=alpha("red",0.6),lwd=2)#Normal Lines
    lines(xg,c2g,col=alpha("blue",0.6),lwd=2)#Gamma lines
    lines(xg,fg,col=alpha("black",0.6),lwd=2)#Mixture model line
    # lines(xg,fg2,col="green",lwd=2)
    if (onlyAddCurves) return(list(xg=xg,c1g=c1g,c2g=c2g,fg=fg))
  }
  if (calculateAreaDifference) {
    f1 <- approxfun(xg, approxfun(density(x,from=0))(xg)-fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = integrate(f2, min(x[x!=0]), max(x))$value
  } else {
    AreaDifference = NULL
  }

  if (removeZeroes) {
    z = rep(0,length(data))
    z[nonZeroInd] <- z_iter
  } else {
    z = z_iter
  }

  # force prob expression values above the max to stay the same value
  maxdata = data[which.max(z)]
  z[which(data>maxdata)] <- max(z)



  if (addContextData) {
    z <- z[1:nOriginal]
  }
  if (plot){
    if(addContextData){
      points(data[1:nOriginal],z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }else{
      points(data,z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }
  }
  model_bic<-bic(ll,n,5)
  model_aic<-aic(ll,5)
  model_icl_bic<-icl_bic(ll,z,n,5)
  return(list(probExpressed = z,
              propExpressed = n*rho_iter/length(data),
              numExpressed=length(which(z>0.5)),
              mu = mu_iter,
              sd = sqrt(sig2_iter),
              alpha = alpha_iter,
              beta = beta_iter,
              rho = rho_iter,
              niter = niter,
              loglik=ll,
              BIC=model_bic,
              AIC=model_aic,
              ICL_BIC=model_icl_bic,
              AreaDifference=AreaDifference))
}



