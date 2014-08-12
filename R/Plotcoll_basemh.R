
Plotcoll$methods(
		basemh = function(chrfilter=NULL, bplower=NULL, bpupper=NULL, pvallower=NULL, pvalupper=NULL) {

			filter = rep(TRUE, nsnp)
			if(! is.null(chrfilter)) {
				message(paste("Filtering by chromosome, as you requested..."))
				filter = chr[, 1] %in% chrfilter
			}
			if(! is.null(bplower)) {
				message(paste("Filtering by bp lower limit, as you requested..."))
				filter = filter & bp[, 1] >= bplower
			}
			if(! is.null(bpupper)) {
				message(paste("Filtering by bp upper limit, as you requested..."))
				filter = filter & bp[, 1] <= bpupper
			}
			if(! is.null(pvallower)) {
				message(paste("Filtering by pval lower limit, as you requested..."))
				filter = filter & pvals[, 1] > pvallower
			}
			if(! is.null(pvalupper)) {
				message(paste("Filtering by pval upper limit, as you requested..."))
				filter = filter & pvals[, 1] < pvalupper
			}
			filter = which(filter)

			baseplotObj = Mhplot(chr[filter, 1], bp[filter, 1], pvals[filter, 1])
		}
)