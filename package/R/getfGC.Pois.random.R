getfGC.Pois.random = function(gcfitj, gctemp, Yj, Nj, betatemp, numGroup, verbose.plot = FALSE, gctemp.keep, min.prop){
	BIC = rep(NA,length(numGroup))
	logL = rep(NA,length(numGroup))
	fGCi.hat = vector('list',length(numGroup))
	fGCi.keep.hat = vector('list',length(numGroup))
	Z.hat = vector('list',length(numGroup))
	vec_pi.hat = vector('list',length(numGroup))
	
	if(verbose.plot){
	  par(mfrow=c(5,2))
	  smoothScatter(gctemp, gcfitj, main = 'Original')
	}
	
	for (K in numGroup){
		cat(paste0("K = ", K, ' '))
		if(K==1){
			Z = matrix(nrow = length(gcfitj), ncol = K, data = 1/K)
			spl = smooth.spline(gctemp[Yj>0], gcfitj[Yj>0], spar = 0.9)
			wt = 1
			fGCi = pmax(1e-20, predict(spl, gctemp)$y * wt)
			vec_pi = rep(1/K, K)

			if(verbose.plot){
			  smoothScatter(gctemp, gcfitj, xlab = 'GC content', ylab = 'Y/beta/N', nrpoints = 0)
			  for(k in 1:K){
			    points(gctemp[order(gctemp)], fGCi[order(gctemp)], lty = 2, col = k, type = 'l', lwd = 2)
			    points(gctemp[which((round(Z))[,k]==1)], gcfitj[which((round(Z))[,k]==1)], cex = 0.4, col = k, pch = 16)
			  }
			}
		} else {
			spl = smooth.spline(gctemp[Yj>0], gcfitj[Yj>0], spar = 0.9)
			wt = runif(1, 3/5, 1.2)
			fGCi = pmax(1e-20, predict(spl, gctemp)$y * wt)
			vec_pi = rep(1/K, K)
			Z = matrix(nrow = length(gcfitj), ncol = K, data = NA)
		
			diff.GC = Inf
			diff.Z = Inf
			diff_LL = Inf
			iter = 0
			max_iter = 200
			break_iter = 1e3
			while(iter <= max_iter){
			  if( iter >= break_iter ){
			    break
			  }
				Etemp = Estep.Pois(Yj, Nj, betatemp, fGCi, vec_pi, min.prop)
				obs_LL = Etemp$obs_LL
				Z.new = Etemp$Z
				
				Mtemp = Mstep.Pois(Z.new, gcfitj, gctemp)
				vec_pi.new = Mtemp$vec_pi
				fGCi.new = Mtemp$fGCi
				
				if( !Etemp$keep_going ){
					break
				}
				
				if(iter > 0){
					diff.GC = sum((fGCi-fGCi.new)^2)
					diff.Z = sum((Z-Z.new)^2)
				}
				
				# Check convergence
				if(iter > 0){
					diff_LL = obs_LL - last_LL
					if(diff_LL > 0 && diff_LL < 1e-5){
						message('Converged! ', 'diff.LL = ', diff_LL, '\n')
						break
					} 
					else if(diff_LL <= 0){
						break
					}
				}
				
				vec_pi = vec_pi.new
				Z = Z.new
				fGCi = fGCi.new
				
				if(iter > 0 && verbose.plot){
					cat('Iter', iter, ':', 'diff.GC = ', diff.GC, 'diff.Z = ', diff.Z, 'diff.LL = ', diff_LL, '\n')
				}
				
				last_LL = obs_LL
				iter = iter + 1
			}
			
			membership = apply(Z, 1, which.max)

			if(verbose.plot){
			  if(Etemp$keep_going){
			    smoothScatter(gctemp, gcfitj, xlab = 'GC content', ylab = 'Y/beta/N', nrpoints = 0)
				spl = smooth.spline(gctemp[Yj>0], (gcfitj / (Z %*% as.matrix((1:K)/2)))[Yj>0], spar = 0.9)
				fGCi = pmax(1e-20, predict(spl, gctemp)$y)
			    for(k in 1:K){
			      points(gctemp[order(gctemp)], fGCi[order(gctemp)]*k/2, lty = 2, col = k, type = 'l', lwd = 2)
			      points(gctemp[membership==k], gcfitj[membership==k], cex = 0.4, col = k, pch = 16)
			    }
			  }
			}
		}
		cat('\n')
		
		if(K==1){
			BIC[which(numGroup==K)]=2*sum(dpois(Yj, lambda=Nj*betatemp*fGCi, log=T))-(length(gcfitj)-spl$df)*log(length(gcfitj))
			logL[which(numGroup==K)]=sum(dpois(Yj, lambda=Nj*betatemp*fGCi, log=T))
		} else {
			if(Etemp$keep_going){
				obs_LL = Estep.Pois(Yj, Nj, betatemp, fGCi, vec_pi, min.prop)$obs_LL
				BIC[which(numGroup==K)] = 2 * obs_LL - (length(gcfitj)-spl$df+K-1)*log(length(gcfitj))
				logL[which(numGroup==K)] = obs_LL
			}
		}
		
		# fit f(GC) given all gctemp, including outliers
		if(K==1){
			spl = smooth.spline(gctemp[Yj>0], gcfitj[Yj>0], spar = 0.9)
			fGCi.keep = pmax(1e-20, predict(spl, gctemp.keep)$y)
		} else{
			spl = smooth.spline(gctemp[Yj>0], (gcfitj / (Z %*% as.matrix((1:K)/2)))[Yj>0], spar = 0.9)
			fGCi.keep = pmax(1e-20, predict(spl, gctemp.keep)$y)
		}
		
		
		fGCi.hat[[which(numGroup==K)]] = fGCi
		fGCi.keep.hat[[which(numGroup==K)]] = fGCi.keep
		Z.hat[[which(numGroup==K)]] = Z
		vec_pi.hat[[which(numGroup==K)]] = vec_pi

	}

	if(verbose.plot){
	  plot(numGroup, BIC, type = 'b', xlab = 'K', ylab = 'BIC')
	  abline(v = numGroup[which.max(BIC)], lty = 2)
	  plot(numGroup, logL, type = 'b', xlab = 'K', ylab = 'logLikelihood')
	  abline(v = numGroup[which.max(logL)], lty = 2)
	}

	K = numGroup[which.max(BIC)]

	return(list(BIC = BIC, logL = logL, fGCi = fGCi, Z = Z, vec_pi = vec_pi, K = K, fGCi.keep = fGCi.keep))
}
