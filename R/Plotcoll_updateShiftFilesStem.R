# TODO: Add comment
# 
# Author: kaiyin
###############################################################################



Plotcoll$methods(
		updateShiftFilesStem = function() {
			shiftFilesStem <<- Sys.glob(paste(bedstem, "_shift_*.bed", sep=""))
			shiftFilesStem <<- getstem(shiftFilesStem)
			names(shiftFilesStem) <<- NULL
			nshiftStrs <<- getNshiftStr(shiftFilesStem)
			nshift <<- as.integer(nshiftStrs)
		}
)
