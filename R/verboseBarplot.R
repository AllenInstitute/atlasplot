
# WGCNA's verboseBarplot
# copied here to avoid installing entire package to use local
# script
.verboseBarplot <- function (x, g, main = "", xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, 
          cex.lab = 1.5, cex.main = 1.5, color = "grey", numberStandardErrors = 1, 
          KruskalTest = TRUE, AnovaTest = FALSE, two.sided = TRUE, 
          addCellCounts = FALSE, horiz = FALSE, ...) 
{
    stderr1 = function(x) {
        sqrt(var(x, na.rm = TRUE)/sum(!is.na(x)))
    }
    SE = tapply(x, factor(g), stderr1)
    err.bp = function(dd, error, two.sided = FALSE, numberStandardErrors, 
                      horiz = FALSE) {
        if (!is.numeric(dd)) {
            stop("All arguments must be numeric")
        }
        if (is.vector(dd)) {
            xval = (cumsum(c(0.7, rep(1.2, length(dd) - 1))))
        }
        else {
            if (is.matrix(dd)) {
                xval = cumsum(array(c(1, rep(0, dim(dd)[1] - 
                                                 1)), dim = c(1, length(dd)))) + 0:(length(dd) - 
                                                                                        1) + 0.5
            }
            else {
                stop("First argument must either be a vector or a matrix")
            }
        }
        MW = 0.25 * (max(xval)/length(xval))
        NoStandardErrors = 1
        ERR1 = dd + numberStandardErrors * error
        ERR2 = dd - numberStandardErrors * error
        if (horiz) {
            for (i in 1:length(dd)) {
                segments(dd[i], xval[i], ERR1[i], xval[i])
                segments(ERR1[i], xval[i] - MW, ERR1[i], xval[i] + 
                             MW)
                if (two.sided) {
                    segments(dd[i], xval[i], ERR2[i], xval[i])
                    segments(ERR2[i], xval[i] - MW, ERR2[i], xval[i] + 
                                 MW)
                }
            }
        }
        else {
            for (i in 1:length(dd)) {
                segments(xval[i], dd[i], xval[i], ERR1[i])
                segments(xval[i] - MW, ERR1[i], xval[i] + MW, 
                         ERR1[i])
                if (two.sided) {
                    segments(xval[i], dd[i], xval[i], ERR2[i])
                    segments(xval[i] - MW, ERR2[i], xval[i] + MW, 
                             ERR2[i])
                }
            }
        }
    }
    if (is.na(ylab)) 
        ylab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(xlab)) 
        xlab = as.character(match.call(expand.dots = FALSE)$g)
    Means1 = tapply(x, factor(g), mean, na.rm = TRUE)
    if (length(unique(x)) > 2) {
        p1 = signif(kruskal.test(x ~ factor(g))$p.value, 2)
        if (AnovaTest) 
            p1 = signif(anova(lm(x ~ factor(g)))$Pr[[1]], 2)
    }
    else {
        p1 = tryCatch(signif(fisher.test(x, g, alternative = "two.sided")$p.value, 
                             2), error = function(e) {
                                 NA
                             })
    }
    if (AnovaTest | KruskalTest) 
        main = paste(main, "p =", p1)
    ret = barplot(Means1, main = main, col = color, xlab = xlab, 
                  ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, 
                  cex.main = cex.main, horiz = horiz, ...)
    if (addCellCounts) {
        cellCountsF = function(x) {
            sum(!is.na(x))
        }
        cellCounts = tapply(x, factor(g), cellCountsF)
        mtext(text = cellCounts, side = if (horiz) 
            2
            else 1, outer = FALSE, at = ret, col = "darkgrey", las = 2, 
            cex = 0.8, ...)
    }
    abline(h = 0)
    if (numberStandardErrors > 0) {
        err.bp(as.vector(Means1), as.vector(SE), two.sided = two.sided, 
               numberStandardErrors = numberStandardErrors, horiz = horiz)
    }
    attr(ret, "height") = as.vector(Means1)
    attr(ret, "stdErr") = as.vector(SE)
    invisible(ret)
}