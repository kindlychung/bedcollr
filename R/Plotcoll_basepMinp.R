
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
