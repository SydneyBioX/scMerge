#' @title Single Cell Stably Express Gene Index
#' @description This function computes the single-cell Stably Expressed Gene (scSEG) 
#' index from Lin. et al. (2019) for a given single-cell count data matrix. Each gene in the data is 
#' fitted with a gamma-normal mixture model and the final SEG index is computed as an average of 
#' key parameters that measure the expression stability of a gene. 
#' 
#' We recommend using either the pre-computed genes (see "See Also" below) or the top SEG genes 
#' from an user's own data as the control genes in the \code{scMerge} function 
#' (see the \code{ctl} argument in the \code{scMerge} function).
#'  
#' @author Shila Ghazanfar, Yingxin Lin, Pengyi Yang
#' @param exprs_mat A log-transformed single-cell data, assumed to have no batch effect and covered a wide range of cell types. 
#' A n by m matrix, where n is the number of genes and m is the number of cells.
#' @param cell_type A vector indicating the cell type information for each cell in the gene expression matrix. 
#' If it is \code{NULL}, the function calculates the scSEG index without using F-statistics.
#' @param BPPARAM A \code{BiocParallelParam} class object from the \code{BiocParallel} package is used. Default is SerialParam(progressbar = TRUE).
#' @param return_all Default to FALSE. If set to TRUE, then all genes will be returned. Otherwise, 
#' only genes with less than 80 percent zeroes will be returned. 
#' @return Returns a data frame. 
#' Each row is a gene and each column is a statistic relating to the stability of expression of each gene.
#' The main statistic is the \code{segIdx} column, which is the SEG index.
#' @importFrom BiocParallel SerialParam
#' @importFrom distr igamma
#' @export
#' @examples
#' ## Loading example data
#' data('example_sce', package = 'scMerge')
#' ## subsetting genes to illustrate usage.
#' exprs_mat = SummarizedExperiment::assay(example_sce, 'logcounts')[1:110, 1:20]
#' set.seed(1)
#' result1 = scSEGIndex(exprs_mat = exprs_mat)
#' ## If parallelisation is needed:
#' param = BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)
#' result2 = scSEGIndex(exprs_mat = exprs_mat, BPPARAM = param)
#' ## Closing the parallelisation
#' BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
#' @references Evaluating stably expressed genes in single cells (2019). doi:10.1093/gigascience/giz106. 
#' @seealso Download human SEG directly from this \href{http://www.maths.usyd.edu.au/u/pengyi/software/scHK/scHK_human.xlsx}{link}; 
#' Download mouse SEG directly from this \href{http://www.maths.usyd.edu.au/u/pengyi/software/scHK/scHK_mouse.xlsx}{link}.




# This is the main function for calculating stably expressed
# gene index
scSEGIndex <- function(exprs_mat, cell_type = NULL, BPPARAM = SerialParam(progressbar = TRUE),
                       return_all = FALSE) {
  
  if (is.null(exprs_mat)) {
    stop("exprs_mat is NULL.")
  }
  
  
  
  if (is.null(cell_type)) {
    message("Calculating scSEG index without F-statistics \n")
  } else {
    if (length(cell_type) != ncol(exprs_mat)) {
      stop("length of cell type information is not equal to the number of column of exprs_mat.")
    }
    cell_type <- as.factor(cell_type)
    message("Calculating scSEG index with F-statistics \n")
  }
  
  ## Core feature 1: mixture modelling Zero%, remove those that
  ## have more than 80% missing values before fitting the
  ## mixture model rownames(exprs_mat) <-
  ## toupper(rownames(exprs_mat))
  
  z.all <- rowMeans(exprs_mat == 0)
  del <- which(z.all > 0.8)
  if (length(del) > 0) {
    exprs_mat_filt <- exprs_mat[-del, ]
    if (nrow(exprs_mat_filt) < 100) {
      stop("Not enough gene pass the QC!")
    }
  } else {
    exprs_mat_filt <- exprs_mat
  }
  
  
  
  # fitting the model
  message("Fitting the mixture model... \n")
  paraMat <- make_para_gn_parallel(as.matrix(exprs_mat_filt), 
                                    BPPARAM = BPPARAM)
  r <- paraMat$rho
  s <- paraMat$sigma
  m <- paraMat$mu
  m.scaled <- (m - min(m))/(max(m) - min(m))
  
  names(r) <- names(s) <- names(m) <- names(m.scaled) <- rownames(exprs_mat_filt)
  
  
  genes <- rownames(exprs_mat_filt)
  z <- z.all[genes] * m.scaled[genes]
  
  ### Fitting ANOVA model for F-stats
  if (is.null(cell_type)) {
    
    ### Combine all core features for creating scSEG index
    x1 = rank(r)/(length(r) + 1)
    x2 = 1 - rank(s)/(length(s) + 1)
    x3 = 1 - rank(z)/(length(z) + 1)
    
    segIdx <- base::rowMeans(cbind(x1, x2, x3))
    
    resMat <- data.frame(segIdx = segIdx, rho = r, sigma = s, 
                         mu = m, mu.scaled = m.scaled, zero = z)
    
  } else {
    
    message("Fitting ANOVA model for F-stats... \n")
    aovStats <- apply(exprs_mat_filt, 1, function(x) {
      tryCatch(stats::aov(as.numeric(x) ~ cell_type), error = function(e) {
        NULL
      })
    })
    f <- log2(
      as.numeric(lapply(aovStats, FUN = 
                          function(x) {
                            tryCatch(summary(x)[[1]]$`F value`[1], 
                                     error = function(e) {NA})}))
    )
    
    
    
    ### Combine all core features for creating scSEG index
    x1 = rank(r)/(length(r) + 1)
    x2 = 1 - rank(s)/(length(s) + 1)
    x3 = 1 - rank(z)/(length(z) + 1)
    x4 = 1 - rank(f)/(length(f) + 1)
    segIdx <- base::rowMeans(cbind(x1, x2, x3, x4))
    resMat <- data.frame(segIdx = segIdx, rho = r, sigma = s, 
                         mu = m, mu.scaled = m.scaled, zero = z, f_stats = f)
    
  }
 if(return_all){
   resMat2 = cbind(gene = as.character(rownames(resMat)), resMat)
   all_genes = data.frame(gene = as.character(rownames(exprs_mat)),
                       zero.all = rowMeans(exprs_mat == 0))
   result = merge(x = all_genes, y = resMat2, by = "gene", all.x = TRUE)
   rownames(result) = rownames(exprs_mat)
 } else {
   result = resMat
 }
  return(result)
}


compute_data_para = function(exprs_mat_gene){
  exprs_mat_gene
  para <- gammaNormMix(exprs_mat_gene, verbose = FALSE, 
                       plot = FALSE)
  if (sum(is.na(para)) == 0) {
    data_para <- c(mu = para$mu, sigma = para$sd, rho = para$rho)
  } else {
    data_para <- c(mu = NA, sigma = NA, rho = NA)
  }
  return(data_para)
}

# This internal function perform mixture model fitting using
# multiple CPUs
make_para_gn_parallel <- function(exprs_mat, BPPARAM) {
  res_list <- bplapply(
    X = seq_len(nrow(exprs_mat)),
    FUN = function(i){compute_data_para(exprs_mat[i, ])},
    BPPARAM = BPPARAM)
  res <- do.call(rbind, res_list)
  res <- data.frame(res)
  rownames(res) = rownames(exprs_mat)
  colnames(res) <- c("mu", "sigma", "rho")
  return(res)
}


bic <- function(loglik, n, p) {
  return(-2 * loglik + p * log(n))
}


aic <- function(loglik, p) {
  return(-2 * loglik + 2 * p)
}  #Tend to fit more component


icl_bic <- function(loglik, postprob, n, p) {
  postprob <- postprob[postprob > 0]
  EN = -sum(postprob * log(postprob))
  return(-2 * loglik + 2 * EN + p * log(n))
}




gammaNormMix <- function(data, thresh = 1e-07, maxiter = 10000, 
                         removeZeroes = TRUE, plot = TRUE, hist = TRUE, hist_col = "light cyan", 
                         verbose = FALSE, forceExponential = FALSE, calculateAreaDifference = FALSE, 
                         minDataPoints = 5, onlyAddCurves = FALSE, addContextData = FALSE, 
                         contextData = NULL) {
  
  
  # fitting a 2 component normal and gamma mixture model
  
  # add other data to fit the model with as well, but only
  # return the classification for those we're interested in
  
  if (addContextData) {
    nOriginal = length(data)
    data <- c(data, contextData)
  }
  
  # assume all values exactly zero already belong to the gamma
  # comp and remove them from the EM algorithm
  
  if (removeZeroes) {
    nonZeroInd = which(data > 0)
    x = data[nonZeroInd]
  } else {
    x = data
  }
  
  if (length(x) < minDataPoints) {
    if (verbose) 
      cat("Not enough data points to fit mixture model!")
    return(NA)
  }
  
  # initiate
  n = length(x)
  z = stats::rbinom(n, 1, 0.5)
  if(sum(z) == 0){z[1] = 1} ## Break out of a sequence of zeroes error
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
  
  while (any(c(abs(mu - mu_iter) > thresh, abs(sig2 - sig2_iter) > 
               thresh, abs(alpha - alpha_iter) > thresh, abs(beta - 
                                                             beta_iter) > thresh, abs(rho - rho_iter) > thresh)) & 
         (niter < maxiter)) {
    
    # save old parameters
    mu = mu_iter
    sig2 = sig2_iter
    alpha = alpha_iter
    beta = beta_iter
    rho = rho_iter
    if (forceExponential) 
      alpha_iter = 1
    
    niter = niter + 1
    
    # M step
    mu_iter = sum(z_iter * x)/sum(z_iter)
    sig2_iter = sum(z_iter * (x - mu_iter) * (x - mu_iter))/sum(z_iter)
    if (sig2_iter <= 0 | is.na(sig2_iter)) 
      sig2_iter = 1e-11
    beta_iter = alpha_iter * sum(1 - z_iter)/sum((1 - z_iter) * 
                                                   x)
    if (beta_iter <= 0 | is.na(beta_iter)) 
      beta_iter = 3
    if (!forceExponential) {
      alpha_iter = distr::igamma(sum((log(beta_iter) + 
                                        log(x)) * (1 - z_iter))/sum(1 - z_iter))
    }
    if (alpha_iter > 150 | is.na(alpha_iter)) 
      alpha_iter = 150
    rho_iter = sum(z_iter)/n
    
    
    # E step
    eta_iter = -0.5 * log(2 * pi * sig2_iter) - ((x - mu_iter) * 
                                                   (x - mu_iter))/(2 * sig2_iter) - alpha_iter * log(beta_iter) + 
      log(gamma(alpha_iter)) - (alpha_iter - 1) * log(x) + 
      beta_iter * x + log(rho_iter/(1 - rho_iter))
    z_iter = 1/(1 + exp(-eta_iter))
    
    if (verbose) 
      cat(niter, mu_iter, sqrt(sig2_iter), alpha_iter, 
          beta_iter, rho_iter, "\n")
  }
  
  
  ll <- sum(log(rho_iter * stats::dnorm(x, mu_iter, sqrt(sig2_iter)) + 
                  (1 - rho_iter) * stats::dgamma(x, shape = alpha_iter, 
                                                 rate = beta_iter)))
  
  
  xg <- seq(0, max(x) + 1, length.out = 300)
  c1g <- rho_iter * stats::dnorm(xg, mu_iter, sqrt(sig2_iter))
  
  c2g <- (1 - rho_iter) * stats::dgamma(xg, shape = alpha_iter, 
                                        rate = beta_iter)
  fg <- rho_iter * stats::dnorm(xg, mu_iter, sqrt(sig2_iter)) + 
    (1 - rho_iter) * stats::dgamma(xg, shape = alpha_iter, 
                                   rate = beta_iter)
  
  if (plot) {
    if (hist) {
      hist(x, probability = TRUE, col = hist_col, breaks = 50, 
           main = NA, xlab = NA, ylab = "Density (zeroes removed)", 
           ylim = c(0, 0.6), xlim = c(0, 20))
    }
    if (!onlyAddCurves) {
      graphics::lines(stats::density(x, from = 0), lty = 2, 
                      lwd = 2, col = scales::alpha("darkgrey", 0.6))
    }
    graphics::lines(xg, c1g, col = scales::alpha("red", 0.6), lwd = 2)  #Normal Lines
    graphics::lines(xg, c2g, col = scales::alpha("blue", 0.6), lwd = 2)  #Gamma lines
    graphics::lines(xg, fg, col = scales::alpha("black", 0.6), lwd = 2)  #Mixture model line
    
    if (onlyAddCurves) 
      return(list(xg = xg, c1g = c1g, c2g = c2g, fg = fg))
  }
  if (calculateAreaDifference) {
    f1 <- stats::approxfun(xg, (stats::approxfun(stats::density(x, 
                                                                from = 0)))(xg) - fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = stats::integrate(f2, min(x[x != 0]), 
                                      max(x))$value
  } else {
    AreaDifference = NULL
  }
  
  if (removeZeroes) {
    z = rep(0, length(data))
    z[nonZeroInd] <- z_iter
  } else {
    z = z_iter
  }
  
  # force prob expression values above the max to stay the same
  # value
  maxdata = data[which.max(z)]
  z[which(data > maxdata)] <- max(z)
  
  
  
  if (addContextData) {
    z <- z[seq_len(nOriginal)]
  }
  if (plot) {
    if (addContextData) {
      graphics::points(data[seq_len(nOriginal)], z * 0, 
                       pch = "|", cex = 1, col = scales::alpha(grDevices::rgb(z, 
                                                                              0, 1 - z), 0.4))
    } else {
      graphics::points(data, z * 0, pch = "|", cex = 1, 
                       col = scales::alpha(grDevices::rgb(z, 0, 1 - z), 0.4))
    }
  }
  model_bic <- bic(ll, n, 5)
  model_aic <- aic(ll, 5)
  model_icl_bic <- icl_bic(ll, z, n, 5)
  return(list(probExpressed = z, propExpressed = n * rho_iter/length(data), 
              numExpressed = length(which(z > 0.5)), mu = mu_iter, 
              sd = sqrt(sig2_iter), alpha = alpha_iter, beta = beta_iter, 
              rho = rho_iter, niter = niter, loglik = ll, BIC = model_bic, 
              AIC = model_aic, ICL_BIC = model_icl_bic, AreaDifference = AreaDifference))
}



