
Plotcoll$methods(
		snp2info = function() {
			chr2 <<- shiftDataByCol(chr)
			snp2 <<- shiftDataByCol(snp)
			bp2 <<-shiftDataByCol(bp)
		}
)