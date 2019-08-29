#' @title Get GC content bias fitting using Expectation-Maximization algorithm for each cell
#'
#' @description Fit a Poisson generalized linear model to normalize the raw read depth data from
#'  single-cell DNA sequencing, under the case-control setting. SCOPE
#'  implements an EM algorithm with random initialization to unmask null regions. This is a
#'  bulit-in function within \code{\link{getfGC}}.
#'
#' @param gcfit.tempj vector of bin-specific GC content biases for each single cell
#' @param gctemp a vector giving values of bin-specific GC content after \code{gcfit.tempj}
#'  outlier removal (for better fitting)
#' @param Yj read depth vector for each single cell
#' @param offsetj vector of offset terms for each cell
#' @param T a vector of integers indicating number of CNV groups. Use BIC to select optimal
#'  number of CNV groups. If \code{T = 1}, assume all reads are from normal regions so that
#'  EM algorithm is not implemented. Otherwise, we assume there is always a CNV group of heterozygous
#'  deletion and a group of null region. The rest groups are representative of different duplication
#'  states.
#' @param draw.plot logical, whether to plot GC content bias fitting results using EM and the
#'  choice of optimal CNV group. Default is \code{NULL}
#' @param alphaj a vector giving values of adjusted/absolute copy number value for cell \emph{j}
#' @param minCountQC the minimum read coverage required for EM fitting. Defalut is \code{20}
#'
#' @return A list with components
#'   \item{fGCi.obj}{A list with EM estimated vector of GC content bias fitting}
#'   \item{Z.obj = Z.obj}{A list with matrix of optimized missing data using EM algorithm}
#'   \item{vec_pi.obj}{A list with vector of EM estimated incident rate for CNV events that span bin \emph{i}}
#'   \item{bin.filter}{Indices of potential outlier bins}
#'   \item{loglik}{Observed log-Likelihood}
#'   \item{BIC}{BIC for CNV group selection}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import graphics stats
#' @export
getfGCj = function(gcfit.tempj, gctemp, Yj, offsetj, T, draw.plot = NULL, alphaj, minCountQC){
  alphaj = pmax(1,round(alphaj))
  if(is.null(draw.plot)){draw.plot = FALSE}
  loglik = BIC = rep(NA,length(T))

  fGCi.obj = vector('list', length(T))
  Z.obj = vector('list', length(T))
  vec_pi.obj = vector('list', length(T))
  fGCi = fitGC(gctemp, gcfit.tempj)
  resid = abs(gcfit.tempj-fGCi)

  # bin.filter is only used for plotting to remove potential outliers
  bin.filter = which(resid > (median(resid)+5*mad(resid)) | Yj < minCountQC)
  if(length(bin.filter)==0){bin.filter = which.max(gcfit.tempj)}
  if(draw.plot){
    par(mfrow = c(5,2))
    smoothScatter(gctemp[-bin.filter], gcfit.tempj[-bin.filter], main = 'Original', xlab = 'GC content', ylab = 'Y/beta/N/exp(gxh)')
  }
  for(Ti in T){
    if(Ti==1){
      Z = matrix(nrow = length(gcfit.tempj), ncol = Ti, data = 1/Ti)
      vec_pi = 1
      loe.fit.temp = loess(gcfit.tempj[-bin.filter]~gctemp[-bin.filter])
      fGCi = predict(loe.fit.temp, newdata = gctemp, se = TRUE)$fit
      temp = min(fGCi[!is.na(fGCi) & fGCi>0])
      fGCi[fGCi <= 0 | is.na(fGCi)] = temp
    }
    if(Ti>=2){
      Z = matrix(nrow = length(gcfit.tempj), ncol = Ti, data = 0)
      mintemp = pmin(Ti,alphaj)
      Z[cbind(1:nrow(Z), mintemp)] = 1
      vec_pi = colSums(Z)/nrow(Z)

      loe.fit.temp = loess((gcfit.tempj/(Z%*%as.matrix((1:Ti)/2)))[-bin.filter]~gctemp[-bin.filter])
      fGCi = predict(loe.fit.temp, newdata = gctemp, se = TRUE)$fit
      temp = min(fGCi[!is.na(fGCi) & fGCi>0])
      fGCi[fGCi <= 0 | is.na(fGCi)] = temp

      diff.GC = Inf; diff.Z = Inf
      iter = 1
      while(iter<= 3 | diff.GC > 5e-6 | diff.Z > 5e-3){
        Mtemp = Mstep(Z, gcfit.tempj, gctemp)
        vec_pi.new = Mtemp$vec_pi
        fGCi.new = Mtemp$fGCi
        Z.new = Estep(fGCi.new, vec_pi.new, Yj, offsetj)
        diff.GC = sum((fGCi-fGCi.new)^2)/length(fGCi)
		diff.Z = sum((Z.new-Z)^2)/length(Z)
        vec_pi = vec_pi.new
        Z = Z.new
        fGCi = fGCi.new
        iter = iter+1
        if(iter >=50) break
      }
    }

    if(Ti==1){
      loe.fit.plot = loess(gcfit.tempj~gctemp)
      fGCi.plot = loe.fit.plot$fitted
      temp = min(fGCi.plot[!is.na(fGCi.plot) & fGCi.plot>0])
      fGCi.plot[fGCi.plot <= 0 | is.na(fGCi.plot)] = temp
      df = predict(loe.fit.plot, newdata = gctemp, se = TRUE)$df

      loglik[which(T==Ti)] = sum(dpois(Yj[-bin.filter], lambda = (offsetj*fGCi)[-bin.filter], log = T))
      BIC[which(T==Ti)] = 2*loglik[which(T==Ti)]-(length(gcfit.tempj)-df)*log(length(gcfit.tempj))
    } else{
      loe.fit.plot = loess((gcfit.tempj/(Z%*%as.matrix((1:Ti)/2)))~gctemp)
      fGCi.plot = loe.fit.plot$fitted
      temp = min(fGCi.plot[!is.na(fGCi.plot) & fGCi.plot>0])
      fGCi.plot[fGCi.plot <= 0 | is.na(fGCi.plot)] = temp
      df = predict(loe.fit.plot, newdata = gctemp, se = TRUE)$df

      loglik[which(T==Ti)] = sum(dpois(Yj[-bin.filter], lambda = (offsetj*fGCi*(Z%*%as.matrix((1:Ti)/2)))[-bin.filter], log = T))
      BIC[which(T==Ti)] = 2*loglik[which(T==Ti)]-(length(gcfit.tempj)-df+Ti-1)*log(length(gcfit.tempj))
    }

    if(draw.plot){
      smoothScatter(gctemp[-bin.filter],gcfit.tempj[-bin.filter], xlab = 'GC content', ylab = 'Y/beta/N/exp(gxh)',nrpoints  =  0, main = paste('T =',Ti))
      if(Ti == 1){
        points(gctemp[order(gctemp)], fGCi.plot[order(gctemp)], lty = 2, col = 2, type = 'l', lwd = 2)
        points(gctemp, gcfit.tempj, cex = 0.4, col = 2, pch = 16)
      } else{
        for(k in 1:Ti){
          points(gctemp[order(gctemp)], fGCi.plot[order(gctemp)]*k/2, lty = 2, col = k, type = 'l', lwd = 2)
          points(gctemp[which((round(Z))[,k]==1)], (gcfit.tempj)[which((round(Z))[,k]==1)], cex = 0.4, col = k, pch = 16)
        }
      }
    }

    fGCi.obj[[which(T==Ti)]] = fGCi
    Z.obj[[which(T==Ti)]] = Z
    vec_pi.obj[[which(T==Ti)]] = vec_pi
  }
  if(draw.plot){
    plot(T, loglik, type = 'b', xlab = 'T',ylab = 'loglik', main = 'Log likelihood')
    abline(v = which.max(loglik), lty = 2)
    plot(T, BIC, type = 'b', xlab = 'T',ylab = 'BIC', main = 'BIC')
    abline(v = which.max(BIC), lty = 2)
    par(mfrow = c(1,1))
  }
  return(list(fGCi.obj = fGCi.obj, Z.obj = Z.obj, vec_pi.obj = vec_pi.obj, bin.filter = bin.filter, loglik = loglik, BIC = BIC))
}
