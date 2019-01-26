Mstep.Pois = function(Z, gcfitj, gctemp){
	vec_pi = apply(Z, 2, mean)
	gcfit.temp = gcfitj / (Z %*% as.matrix((1:ncol(Z))/2))
	spl <- smooth.spline(gctemp[gcfitj>0], gcfit.temp[gcfitj>0], spar = 0.9)
	fGCi = pmax(1e-20, predict(spl, gctemp)$y)
	return(list(vec_pi = vec_pi, fGCi = fGCi))
}