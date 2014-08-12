# TODO: Add comment
#
# Author: kaiyin
###############################################################################




require(devtools)
load_all("/Users/kaiyin/personal_config_bin_files/workspace/bedcollr", reset=TRUE)
#load_all("/Users/kaiyin/personal_config_bin_files/workspace/manqq", reset=TRUE)
setwd("/Volumes/wdDataTransfer/data/sskn_regions_from_fan/AgeSexSskn")
o = Plotcoll("sskn_reg")
o$shiftstemCommon
o$shiftFilesStem
o$nshiftStrs
o$nshift
o$readout("assoc.linear")
o.contrast = o$contrastplot()$mhplot()
print(o.contrast)

#f = function(...) {
#	x = match.call(expand.dots = FALSE)$`...`
#	x
#}
#x = f(a=1, b=2, c=3)
#names(x) = paste("--", names(x), sep="")
#names(x)
#x



require(bedcollr)
o = Plotcoll("RS123_1kg")
o$shiftstemCommon
o$shiftFilesStem
o$nshiftStrs
o$nshift
o$readout("assoc.linear")
o.contrast = o$contrastplot()$mhplot()
ggsave(filename = "/tmp/height_pval5e-3.png", plot = o.contrast, height=5, width=8)