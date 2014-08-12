
Plotcoll$methods(
		initialize = function(bedstem) {
			bedstem <<- bedstem
			bedpath <<- paste(bedstem, "bed", sep=".")
			fampath <<- paste(bedstem, "fam", sep=".")
			bimpath <<- paste(bedstem, "bim", sep=".")
			message(paste("Original bed file:", bedpath))
			message(paste("Original fam file:", fampath))
			message(paste("Original bim file:", bimpath))
			
			shiftstemCommon <<- paste(bedstem, "_shift_", sep = "")
			message(paste("Commen prefix of shifted bed files: ", shiftstemCommon))
			updateShiftFilesStem()
			nsnp <<- 0
			nindiv <<- 0
			nTotalShifts <<- length(shiftFilesStem)
			pvals <<- matrix(1, 1, 1)
			chr <<- matrix(, 0, 0)
			snp <<- matrix(, 0, 0)
			bp <<- matrix(, 0, 0)
			chr2 <<- matrix(, 0, 0)
			snp2 <<- matrix(, 0, 0)
			bp2 <<- matrix(, 0, 0)
			chrunique <<- 0
			minPvalsBonfer <<- 1
			minPvals <<- 1
		}
)