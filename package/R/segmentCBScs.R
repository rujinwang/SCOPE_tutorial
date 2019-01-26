segmentCBScs = function(Y, Yhat, sampname, ref, lmax, mode, QCmetric, segment.CODEX2 = FALSE) {
	poolcall = NULL
	lmax = lmax - 1
	message("Cross-sample segmenting for ", ncol(Y), " samples.")
	Z = NULL
	for(sampno in 1:ncol(Y)) {
		message("Calculating scan statistic for sample ", sampno, ": ", sampname[sampno], ".")
		y = Y[, sampno]
		yhat = Yhat[, sampno]
		num = length(y)
		y = c(y, rep(0, lmax))
		yhat = c(yhat, rep(0, lmax))
		i = rep(1:num, rep(lmax, num))
		j = rep(1:lmax, num) + i
		yact = rep(0, length(i))
		lambda = rep(0, length(i))
		for (k in 1:num) {
			yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
			lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax)])[-1]
		}
		i <- i[j <= num]
		yact <- yact[j <= num]
		lambda <- lambda[j <= num]
		j <- j[j <= num]
		
		yact[lambda < 20] <- 20
		lambda[lambda < 20] <- 20
		if (mode == "integer") {
			chat <- round(2 * (yact/lambda))
		} else if (mode == "fraction") {
			chat <- 2 * (yact/lambda)
		}
		lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
		chat[chat > 2*exp(2)] <- 2*exp(2)
		Z = cbind(Z, lratio)
	}
	Z = rowSums(Z)
		
	if (sum(Z > 0) > 0) {
	  if (sum(Z > 0) >= 2) {
	    finalmat <- (cbind(i, j, Z))[Z > 0, ]
	    finalmat <- finalmat[order(-finalmat[, 3]), ]
	    s <- 1
	    while (s <= (nrow(finalmat))) {
	      rowstart <- finalmat[s, 1]
	      rowend <- finalmat[s, 2]
	      rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= rowstart)
	      rowsel[s] <- FALSE
	      finalmat <- finalmat[!rowsel, ]
	      
	      if (is.vector(finalmat)) {
	        finalmat <- t(as.matrix(finalmat))
	      }
	      s <- s + 1
	    }
	  }
	  if (sum(Z > 0) == 1) {
	    finalmat <- (cbind(i, j, Z))[Z > 0, ]
	    finalmat <- t(as.matrix(finalmat))
	  }
	  finalmat <- round(finalmat, digits = 3)
	  loglikeij <- cumsum(finalmat[, 3])
	  mBIC <- rep(NA, length(loglikeij))
	  
	  kappa1 = 3/2
	  kappa2 = 2.27
	  N = ncol(Y)
	  T = num
	  
	  for (s in 1:nrow(finalmat)){
	    tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
	    m <- length(tau) - 2
	    delta = matrix(0, nrow = m, ncol = N)
	    tau2 = tau[-c(which.min(tau), which.max(tau))]
	    for(ind in 1:ncol(Y)){
	      for (tt in 1:length(tau2)) {
	        tau.temp = tau2[tt]
	        
	        chat.idv <- Yhat[tau.temp,ind]
	        chat0.idv <- Yhat[tau.temp-1,ind]
	        delta[tt,ind] = chat.idv - chat0.idv
	      }
	    }
	    M = sum(delta!=0)
	    
	    mbic <- loglikeij[s]
	    mbic = mbic - M/2 * log(2 * loglikeij[s] / M) - log(choose(T, m))
	    if(M==N*m){
	      mbic <- mbic - M/2 - sum(log(rowSums(delta^2))) - m * (kappa1 - kappa2)
	    } else{
	      mbic <- mbic - M/2 - sum(log(rowSums(delta^2))) - m * (kappa1 - kappa2) + (M * log(M/(N*m)) + (N*m - M) * log(1 - M/(N*m)))
	    }
	    
	    mBIC[s] <- mbic
	  }
	  mBIC <- round(mBIC, digits = 3)
	  
    if (mBIC[1] > 0) {
	    finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ])
	    poolcall <- rbind(poolcall, finalmat)
	    colnames(poolcall) = c("st", "end", "Z", "mBIC")
	  }
	}

	if (is.vector(poolcall)) {
		poolcall <- t(as.matrix(poolcall))
	}
	
	finalcall = NULL
	image.orig=log(pmax(0.001,Y)/Yhat)
	image.orig[image.orig<=-2]=-2
	image.orig[image.orig>2]=2
	
	image.seg = matrix(data=0, nrow = nrow(image.orig), ncol = ncol(image.orig))
	colnames(image.seg) = colnames(Y)
	
	for (i in 1:nrow(poolcall)) {
		st_bin = poolcall[i,"st"]
		ed_bin = poolcall[i,"end"]
		st_bp <- start(ref)[as.numeric(poolcall[i,"st"])]
		ed_bp <- end(ref)[as.numeric(poolcall[i,"end"])]
		
		yact.ind = colSums(Y[st_bin:ed_bin,])
		lambda.ind = colSums(Yhat[st_bin:ed_bin,])
		if (mode == "integer") {
			chat.ind <- round(2 * (yact.ind/lambda.ind))
		} else if (mode == "fraction") {
			chat.ind <- 2 * (yact/lambda)
		}
		chat.ind[chat.ind > 2*exp(2)] = 2*exp(2)
		image.seg[poolcall[i,1]:poolcall[i,2], ] = matrix(nrow = (poolcall[i,2] - poolcall[i,1] + 1), 
														  ncol = ncol(image.seg), data = log(chat.ind / 2), byrow = TRUE)
		temp = cbind(colnames(Y), rep(st_bin, ncol(Y)), rep(ed_bin, ncol(Y)), rep(st_bp, ncol(Y)), 
					rep(ed_bp, ncol(Y)), chat.ind)
		temp = temp[chat.ind != 2, ]
		finalcall = rbind(finalcall, temp)
		rownames(finalcall) = NULL
		colnames(finalcall) = c("sample_name", "st_bin", "ed_bin", "st_bp", "ed_bp", "cnv_no")
	}
	
	finalcall = as.data.frame(finalcall)
	
	finalcall$st_bin = as.numeric(paste(finalcall$st_bin))
	finalcall$ed_bin = as.numeric(paste(finalcall$ed_bin))
	finalcall$st_bp = as.numeric(paste(finalcall$st_bp))
	finalcall$ed_bp = as.numeric(paste(finalcall$ed_bp))
	finalcall$cnv_no = as.numeric(paste(finalcall$cnv_no))

	
	finalcall_CODEX2 = NULL
	if(segment.CODEX2){
	  # CODEX2 segmenting by individuals
	  for(sampno in 1:ncol(Y)) {
	    message("Segmenting sample ", sampno, ": ", sampname[sampno], ".")
	    y = Y[, sampno]
	    yhat = Yhat[, sampno]
	    num = length(y)
	    y = c(y, rep(0, lmax))
	    yhat = c(yhat, rep(0, lmax))
	    i = rep(1:num, rep(lmax, num))
	    j = rep(1:lmax, num) + i
	    yact = rep(0, length(i))
	    lambda = rep(0, length(i))
	    for (k in 1:num) {
	      yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
	      lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax)])[-1]
	    }
	    i <- i[j <= num]
	    yact <- yact[j <= num]
	    lambda <- lambda[j <= num]
	    j <- j[j <= num]
	    yact[lambda < 20] <- 20
	    lambda[lambda < 20] <- 20
	    if (mode == "integer") {
	      chat <- round(2 * (yact/lambda))
	    } else if (mode == "fraction") {
	      chat <- 2 * (yact/lambda)
	    }
	    lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
	    chat[chat > 2*exp(2)] <- 2*exp(2)
	    if (sum(lratio > 0) > 0) {
	      if (sum(lratio > 0) >= 2) {
	        finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 0, ]
	        finalmat <- finalmat[order(-finalmat[, 6]), ]
	        s <- 1
	        while (s <= (nrow(finalmat))) {
	          rowstart <- finalmat[s, 1]
	          rowend <- finalmat[s, 2]
	          rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= rowstart)
	          rowsel[s] <- FALSE
	          finalmat <- finalmat[!rowsel, ]
	          if (is.vector(finalmat)) {
	            finalmat <- t(as.matrix(finalmat))
	          }
	          s <- s + 1
	        }
	      }
	      if (sum(lratio > 0) == 1) {
	        finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 0, ]
	        finalmat <- t(as.matrix(finalmat))
	      }
	      finalmat <- round(finalmat, digits = 3)
	      loglikeij <- cumsum(finalmat[, 6])
	      mBIC <- rep(NA, length(loglikeij))
	      for (s in 1:nrow(finalmat)) {
	        tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
	        P <- length(tau) - 2
	        mbic <- loglikeij[s]
	        mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - tau[1:(length(tau) - 1)]))
	        mbic <- mbic + (0.5 - P) * log(num)
	        mBIC[s] <- mbic
	      }
	      mBIC <- round(mBIC, digits = 3)
	      if (mBIC[1] > 0) {
	        finalmat <- cbind(rep(sampname[sampno], nrow(finalmat)), finalmat)
	        finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ])
	        finalcall_CODEX2 <- rbind(finalcall_CODEX2, finalmat)
	      }
	    }
	  }
	  
	  if (is.vector(finalcall_CODEX2)) {
	    finalcall_CODEX2 <- t(as.matrix(finalcall_CODEX2))
	  }
	  colnames(finalcall_CODEX2) <- c("sample_name", "st_bin", "ed_bin", "raw_cov", 
	                                  "norm_cov", "copy_no", "lratio", "mBIC")
	  rownames(finalcall_CODEX2) <- rep("", nrow(finalcall_CODEX2))
	  lratio = as.numeric(finalcall_CODEX2[, "lratio"])
	  rownames(finalcall_CODEX2) = paste("cnv", 1:nrow(finalcall_CODEX2), sep = "")
	  finalcall_CODEX2 = as.data.frame(finalcall_CODEX2)
	  finalcall_CODEX2$copy_no = as.numeric(paste(finalcall_CODEX2$copy_no))
	  finalcall_CODEX2$raw_cov = as.numeric(paste(finalcall_CODEX2$raw_cov))
	  finalcall_CODEX2$norm_cov = as.numeric(paste(finalcall_CODEX2$norm_cov))
	  finalcall_CODEX2$st_bin = as.numeric(paste(finalcall_CODEX2$st_bin))
	  finalcall_CODEX2$ed_bin = as.numeric(paste(finalcall_CODEX2$ed_bin))
	  finalcall_CODEX2$lratio = as.numeric(paste(finalcall_CODEX2$lratio))
	  finalcall_CODEX2$mBIC = as.numeric(paste(finalcall_CODEX2$mBIC))	  
	}

	return(list(poolcall = poolcall, finalcall = finalcall, finalcall_CODEX2 = finalcall_CODEX2, image.orig = image.orig, image.seg = image.seg))
}