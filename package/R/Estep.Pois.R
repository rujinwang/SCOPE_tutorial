Estep.Pois = function(Yj, Nj, betatemp, fGCi, vec_pi, min.prop){
	keep_going = TRUE
	if (any(vec_pi<=min.prop)){
		message("A mixture proportion has gone to zero!!")
		keep_going = FALSE
	}
	if (vec_pi[2]<5e-2){
		message("Normal proportion has gone to zero!!")
		keep_going = FALSE
	}
	Z = matrix(nrow = length(fGCi), ncol = length(vec_pi))
	for (k in 1:length(vec_pi)){
		Z[,k] = log(vec_pi[k]) + dpois(Yj, lambda = Nj * betatemp * fGCi * (k/2), log = T)
	}
	obs_LL = 0
	for (ii in 1:length(fGCi)){
		tmp_num = logsumexp(Z[ii,])
		obs_LL = obs_LL + tmp_num
		Z[ii,] = exp(Z[ii,] - tmp_num)
	}
	return(list(Z = Z, obs_LL = obs_LL, keep_going = keep_going))
}