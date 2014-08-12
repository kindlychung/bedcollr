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