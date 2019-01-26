logsumexp = function(xx){
	if(length(xx)==1){
		xx
	}else{
		max_val = max(xx)
		xx2 = xx - max_val
		log(sum(exp(xx2))) + max_val
	}
}