# TODO: Add comment
# 
# Author: kaiyin
###############################################################################



shiftDataByCol = function(dat) {
	nr = nrow(dat)
	nc = ncol(dat)
	message("Shifting data col by col...")
	message(paste("    dat has", nr, "rows"))
	message(paste("    dat has", nc, "cols"))
	
	message(paste("Head of dat before shifting: "))
	print(head(dat))
	message(paste("Tail of dat before shifting: "))
	print(tail(dat))

	message("initialize a matrix newdat with same size as dat...")
	newdat = matrix(NA, nr, nc)
	for(i in 1:nc) {
		colShiftN = as.integer(colnames(dat)[i])
		message(paste("    This col should be shifted by ", colShiftN, "according to colname"))
		newdat[1:(nr-colShiftN), i] = dat[(colShiftN+1):nr, 1]
	}
	
	message(paste("Head of dat after shifting: "))
	print(head(newdat))
	message(paste("Tail of dat after shifting: "))
	print(tail(newdat))
	
	message(paste("Setting colnames for new dat..."))
	colnames(newdat) = colnames(dat)
	print(head(newdat))
	print(tail(newdat))

	newdat
}
