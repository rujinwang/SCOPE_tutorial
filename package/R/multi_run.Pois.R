multi_run.Pois = function(gcfitj, gctemp, Yj, Nj, betatemp, numGroup, rerun, verbose.plot = FALSE, qc.thres = 5e-5, min.prop = 2e-3){
	gcfitj.keep = gcfitj
	gctemp.keep = gctemp
	Yj.keep = Yj
	betatemp.keep = betatemp
	
	resid.index = which((median(gcfitj)-5*mad(gcfitj)) < gcfitj & gcfitj< (median(gcfitj)+5*mad(gcfitj)) & gcfitj > qc.thres)

	gcfitj = gcfitj[resid.index]
	gctemp = gctemp[resid.index]
	Yj = Yj[resid.index]
	betatemp = betatemp[resid.index]
	
	out=list()
	if(verbose.plot){
		par(mfrow=c(5,2))
		smoothScatter(gctemp, gcfitj, main = 'Original')
	}
	for(Groupi in numGroup){
		if(Groupi==1){
			out[[which(numGroup==Groupi)]] = getfGC.Pois.random(gcfitj, gctemp, Yj, Nj, betatemp, numGroup = Groupi, verbose.plot = FALSE, gctemp.keep, min.prop = min.prop)
			if(verbose.plot){
				fGCi = out[[which(numGroup==Groupi)]]$fGCi
				smoothScatter(gctemp, gcfitj, xlab = 'GC content', ylab = 'Y/beta/N', nrpoints = 0)
				points(gctemp[order(gctemp)], fGCi[order(gctemp)], lty = 2, col = 1, type = 'l', lwd = 2)
				points(gctemp, gcfitj, cex = 0.4, col = 1, pch = 16)
			}
		} else{
			aa = list()
			for (r in seq(rerun)){
				aa[[r]] = getfGC.Pois.random(gcfitj, gctemp, Yj, Nj, betatemp, numGroup = Groupi, gctemp.keep = gctemp.keep, min.prop = min.prop)
			}
			if(sum(is.na(sapply(aa, function(z) z$BIC)))==rerun){
				out[[which(numGroup==Groupi)]] = list(BIC = NA, logL = NA, fGCi = NA, Z = NA, vec_pi = NA, K = Groupi, fGCi.keep = NA)
			} else{
				out[[which(numGroup==Groupi)]] = aa[[which.max(sapply(aa, function(z) z$BIC))]]
			}
			
			if(verbose.plot){
				if(!is.na(out[[which(numGroup==Groupi)]]$BIC)){
					smoothScatter(gctemp, gcfitj, xlab = 'GC content', ylab = 'Y/beta/N', nrpoints = 0)
					fGCi = out[[which(numGroup==Groupi)]]$fGCi
					Z = out[[which(numGroup==Groupi)]]$Z
					membership = apply(Z, 1, which.max)
					for(k in 1:Groupi){
						points(gctemp[order(gctemp)], fGCi[order(gctemp)]*k/2, lty = 2, col = k, type = 'l', lwd = 2)
						points(gctemp[membership==k], gcfitj[membership==k], cex = 0.4, col = k, pch = 16)
					}
				}
			}
		}
	}
	if(verbose.plot){
		plot(numGroup, sapply(out, function(z) z$BIC), type = 'b', xlab = 'K', ylab = 'BIC')
		abline(v = numGroup[which.max(sapply(out, function(z) z$BIC))], lty = 2)
		plot(numGroup, sapply(out, function(z) z$logL), type = 'b', xlab = 'K', ylab = 'logLikelihood')
		abline(v = numGroup[which.max(sapply(out, function(z) z$logL))], lty = 2)
	}
	return(out = out[[which.max(sapply(out, function(z) {z$BIC}))]])
}