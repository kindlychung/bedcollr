
Plotcoll$methods(
		snppairs = function() {
			minpidx = apply(pvals, 1, which.min)
			minpmat = cbind(1:nrow(pvals), minpidx)
			retpairs = cbind(snp[, 1], snp2[minpmat])
			retpairs
		}
)