Plotcoll$methods(
		contrastplot = function(chrfilter=NULL, bplower=NULL,
				bpupper=NULL, pvallower=NULL, pvalupper=NULL, minpcorrect=TRUE
		) {
			
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
			
			message(paste("Combining single-SNP and QCDH results..."))
			bothChr = c(chr[filter, 1], chr[filter, 1])
			bothBp = c(bp[filter, 1], bp[filter, 1])
			if(minpcorrect == TRUE) {
				bothP = c(pvals[filter, 1], minPvalsBonfer[filter])
			} else {
				bothP = c(pvals[filter, 1], minPvals[filter])
			}
			colorvec = rep(c("Single SNP", "QCDH"), each=length(filter))
			
			# get the order right
			message(paste("Reordering data by chr and bp..."))
			posorder = order(bothChr, bothBp)
			bothChr = bothChr[posorder]
			bothBp  = bothBp[posorder]
			bothP   = bothP[posorder]
			colorvec = colorvec[posorder]
			
			contrastPlotObj = Mhplot(bothChr, bothBp, bothP, colorvec=colorvec)
			contrastPlotObj
		}
)