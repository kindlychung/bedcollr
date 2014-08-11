# TODO: Add comment
# 
# Author: kaiyin
###############################################################################



getNshiftStr = function(shiftpath) {
	nshiftStr = gsub(".*_shift_(\\d{4}).*", "\\1", shiftpath)
	nshiftStr
}

getNshift = function(shiftpath) {
	nshiftRes = as.integer(getNshiftStr(shiftpath))
	nshiftRes
}

getShiftStem = function(shiftpath) {
	shiftstem = paste(shiftstemCommon, getNshiftStr(shiftpath), sep="")
	shiftstem
}

getstem = function(pathname) {
	gsub("(.*?)\\.bed$", "\\1", pathname)
}
