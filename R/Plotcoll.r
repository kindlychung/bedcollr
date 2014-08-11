pminNoNA = function(...) {
	pmin(..., na.rm = TRUE)
}

Plotcoll = setRefClass("Plotcoll",
		fields = list(
				bedstem="character",
				shiftstemCommon = "character",
				bedpath="character",
				fampath="character",
				bimpath="character",
				shiftFilesStem="character",
				nsnp="numeric",
				nindiv="numeric",
				nshift = "numeric",
				nTotalShifts = "numeric",
				nshiftStrs = "character",
				chr = "matrix",
				snp = "matrix",
				bp = "matrix",
				chr2 = "matrix",
				snp2 = "matrix",
				bp2 = "matrix",
				pvals = "matrix",
				minPvals = "numeric",
				minPvalsBonfer = "numeric",
				chrunique = "numeric"
		))


Plotcoll$methods(
		updateShiftFilesStem = function() {
			shiftFilesStem <<- Sys.glob(paste(bedstem, "_shift_*.bed", sep=""))
			shiftFilesStem <<- getstem(shiftFilesStem)
			names(shiftFilesStem) <<- NULL
			nshiftStrs <<- getNshiftStr(shiftFilesStem)
			nshift <<- as.integer(nshiftStrs)
		}
)



Plotcoll$methods(
		snp2info = function() {
			chr2 <<- shiftDataByCol(chr)
			snp2 <<- shiftDataByCol(snp)
			bp2 <<-shiftDataByCol(bp)
		}
)

Plotcoll$methods(
		readout = function(tagname) {
			
			if (nsnp == 0) {
				nsnp <<- R.utils::countLines(paste(bedstem, "_shift_0000.", tagname, sep="")) - 1
				message(paste("Number of SNPs: ", nsnp))
			}

			outRdata = paste(bedstem, "RData", sep = ".")
			if(file.exists(outRdata)) {
				message("readout:   ")
				message("Reading from previously generated data...")
				load(outRdata)
				currentEnv = environment()
				
				# make sure dimensions are right
				dim1 = dim(chr)
				dim2 = dim(snp)
				dim3 = dim(bp)
				dim4 = dim(pvals)
				if(
						length(unique(c(dim1[1], dim2[1], dim3[1], dim4[1]))) > 1 |
						length(unique(c(dim1[2], dim2[2], dim3[2], dim4[2]))) > 1 
						) {
					stop(paste(outRdata, "is corrupt, you had better remove it!", sep = " "))
				}
				
				# figure out what files to read
				# add new files
				
				
				message("Reading colnames from chr matrix...")
				hasColnames = colnames(chr)
				addIdx = which(!(nshiftStrs %in% hasColnames))
				addColnames = nshiftStrs[addIdx]
				addFileStems = shiftFilesStem[addIdx]
				
				if(length(addFileStems) > 0) {
					fileAdded = TRUE
					message("It looks like you have added some new files, let me read them: ")
					print(addFileStems)
					
					pcounter = ncol(pvals) + 1
					for(i in 1:length(addFileStems)) {
						addf = addFileStems[i]
						addf = paste(addf, tagname, sep = ".") 
						if(!file.exists(addf)) {
							stop("You should analyze before you can plot!")
						}
						addcoln = addColnames[i]
						message(paste("   Name of the new col: ", addcoln))

						message("Reading ", addf, "...")
						datshift = readplinkoutr(addf)
						currentEnv$chr = cbind(currentEnv$chr, datshift$CHR)
						currentEnv$snp = cbind(currentEnv$snp, datshift$SNP)
						currentEnv$bp = cbind(currentEnv$bp, datshift$BP)
						currentEnv$pvals = cbind(currentEnv$pvals, datshift$P)
						colnames(currentEnv$chr)[pcounter] = addcoln
						colnames(currentEnv$snp)[pcounter] = addcoln
						colnames(currentEnv$bp)[pcounter] = addcoln
						colnames(currentEnv$pvals)[pcounter] = addcoln
						pcounter = pcounter + 1
					}
				} else {
					fileAdded = FALSE
				}
				
				
				# remove colums if corresponding files are removed
				# e.g. if x_shift_0002.bed is removed, then remove column "2"
				# first you need to update hasColnames!
				hasColnames = colnames(chr)
				selectIdx = which(hasColnames %in% nshiftStrs)
				rmIdx = which(!(hasColnames %in% nshiftStrs))
				if(length(selectIdx) != length(hasColnames)) {
					fileRemoved = TRUE
					rmColname = hasColnames[rmIdx]
					rmFiles = sapply(rmColname, function(col) {
								sprintf("%s%04d", shiftstemCommon, as.integer(col))
							})
					names(rmFiles) = NULL
					message("It looks like you have removed some file(s), let me delete corresponding data:")
					print(rmFiles)
					currentEnv$chr = currentEnv$chr[, selectIdx]
					currentEnv$snp = currentEnv$snp[, selectIdx]
					currentEnv$bp = currentEnv$bp[, selectIdx]
					currentEnv$pvals = currentEnv$pvals[, selectIdx]
				} else {
					fileRemoved = FALSE
				}
				
				# assign into the class!
				chr <<- chr
				snp <<- snp
				bp <<- bp
				pvals <<- pvals
				
				if(fileAdded | fileRemoved) {
					message("Saving changes back to RData file...")
					save(chr, snp, bp, pvals, file=outRdata)
				} else {
					message("I did not see any change, keep the RData file as it is.")
				}
			} else {
				# set up matrices
				message("Setting up matrices for pvals, chr, snp and bp...")
				pvalsNcols = length(nshiftStrs)
				message(paste("Need", nsnp, "rows"))
				message(paste("Need", pvalsNcols, "columns"))
				pvals <<- matrix(NA, nsnp, pvalsNcols)
				chr <<- pvals
				snp <<- pvals
				bp <<- pvals
				
				# read shifted results
				message("Get list of plink output files...")
				outfiles = sapply(shiftFilesStem, function(eachstem) {
							paste(eachstem, tagname, sep = ".")
						})
				outfiles = setNames(outfiles, NULL)
				print(outfiles)
				
				for(i in 1:length(outfiles)) {
					message(paste("Processing", outfiles[i], "..."))
					outfile = outfiles[i]
					coln = nshiftStrs[i]
					message(paste("Assign col number for it: ", coln))
					message("Reading ", outfiles[i], "...")
					datshift = readplinkoutr(outfiles[i])
					chr[, i] <<- datshift$CHR
					snp[, i] <<- datshift$SNP
					bp[, i] <<- datshift$BP
					pvals[, i] <<- datshift$P
				}
				
				message("Naming columns by nshift...")
				print(nshiftStrs)
				colnames(chr) <<- colnames(snp) <<- colnames(bp) <<- colnames(pvals) <<- nshiftStrs
				print(head(chr))
				print(head(snp))
				print(head(bp))
				print(head(pvals))
				message("Tail of pvals: ")
				print(tail(pvals))

				message("Saving to ", outRdata, "...")
				save(chr, snp, bp, pvals, file=outRdata)
			}
			
			
			if(chrunique[1] == 0) {
				message("Nothing in the set of unique CHRs, calculating it...")
				chrunique <<- unique(chr[, 1])
			} 
			message("Unique set of CHRs:")
			print(chrunique)
			
			message("Calculating minimal p values...")
			minPvals <<- do.call(pminNoNA, as.data.frame(pvals))
			print(head(minPvals))
			
			message("Calculating number of tests...")
			ntests = apply(pvals, 1, function(rowIter) {sum(!is.na(rowIter))})
			print(head(ntests, 50))
			
			message("Bonferroni correction for minimal p valus...")
			minPvalsBonfer <<- minPvals * ntests
			print(head(minPvalsBonfer, 25))
			message(paste("Prune p values greater than 1 to 1..."))
			minPvalsBonfer <<- ifelse(minPvalsBonfer > 1, 1, minPvalsBonfer)
			print(head(minPvalsBonfer, 25))
			
			message("Update info for SNP2...")
			snp2info()
		}
)

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

Plotcoll$methods(
		minpmh = function(
				chrfilter=NULL, bplower=NULL, bpupper=NULL,
				pvallower=NULL, pvalupper=NULL, minpcorrect=TRUE
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

			if(minpcorrect == TRUE) {
				minpplotObj = Mhplot(chr[filter, 1], bp[filter, 1], minPvalsBonfer[filter])
			} else {
				minpplotObj = Mhplot(chr[filter, 1], bp[filter, 1], minPvals[filter])
			}

			minpplotObj
		}
)

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

Plotcoll$methods(
		basepMinp = function(logpvals=TRUE, minpcorrect=TRUE) {
			df = data.frame(basep = pvals[, 1])
			if(minpcorrect == TRUE) {
				df$minp = minPvalsBonfer
			} else {
				df$minp = minPvals
			}
			
			if(logpvals == TRUE) {
				df$minp = -log10(df$minp)
				df$basep = -log10(df$basep)
			}
			
			bmplot = ggplot(df, aes(basep, minp)) + geom_point(alpha=0.4)
		}
)

Plotcoll$methods(
		snppairs = function() {
			minpidx = apply(pvals, 1, which.min)
			minpmat = cbind(1:nrow(pvals), minpidx)
			retpairs = cbind(snp[, 1], snp2[minpmat])
			retpairs
		}
)

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
			
#			readout("assoc.linear")
		}
)

