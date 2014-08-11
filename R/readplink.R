readplinkoutr = function(filename, colnameSelect=c("CHR", "SNP", "BP", "P"), covarHiden=FALSE) {
	# allow globbing
	filename = Sys.glob(filename)

	# extract colnames
    pilotColnames = read.table(filename, header=TRUE, nrows=1, stringsAsFactors=FALSE)
    cnames = colnames(pilotColnames)

	if(!covarHiden) {
		message("Covariates has not been hidden, need to sort it out...")
		# get index of selected cols
		if(! all(colnameSelect %in% cnames)) {
			stop("Some requested columns unavailable!")
		}
		colSelectIdx = which(cnames %in% colnameSelect)

		# calculate number of variables from the TEST col
		testIdx = which(cnames == "TEST")
		alltests = readcols(filename, testIdx, 1, 1)[,1];
		message("Take a look at the TEST column of the plink output file: ")
		print(head(alltests))

		# assuming nobody would cram 500 variables into a GWAS model
		if (length(alltests) > 500) {
			ncovar = length(unique(alltests[1:500]))
		} else {
			ncovar = length(unique(alltests))
		}
	}
	else {
		message("Covariates has been hidden, cool.")
		ncovar = 1
	}

	message(paste("You have", ncovar, "variables in this output file."))

	plinkout = as.data.frame(readcols(filename, colSelectIdx, 1, ncovar), stringsAsFactors=FALSE)
	plinkoutColnames = cnames[colSelectIdx]
    plinkout = setNames(plinkout, plinkoutColnames)

	# correct data types
	for (i in 1:ncol(plinkout)) {
		if (plinkoutColnames[i] == "CHR" | plinkoutColnames[i] == "NMISS" | plinkoutColnames[i] == "BP") {
			message(paste("Setting", plinkoutColnames[i], "to integer..."))
			plinkout[, i] = as.integer(plinkout[, i])
		}
		else if(plinkoutColnames[i] == "BETA" | plinkoutColnames[i] == "STAT" | plinkoutColnames[i] == "P") {
			message(paste("Setting", plinkoutColnames[i], "to numeric..."))
			plinkout[, i] = as.numeric(plinkout[, i])
		}
	}

    plinkout
}



