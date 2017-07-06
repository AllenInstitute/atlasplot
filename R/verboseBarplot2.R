
#------------------------------HELPER FUNCTIONS--------------------------------#

# Adapted from WGCNA verboseBarplot by Jeremy Miller<jemerym@alleninstitute.org>
# creates a verbose barplot similar to the one at https://www.rdocumentation.org/packages/WGCNA/versions/1.51/topics/verboseBarplot
.verboseBarplot2 <- function (x, g, main = "", xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5,
                              cex.lab = 1.5, cex.main = 1.5, color = "grey", numberStandardErrors = 1,
                              KruskalTest = TRUE, AnovaTest = FALSE, two.sided = TRUE,
                              ...){
    
    if(is.factor(g)) { l = levels(g)
    } else {l = levels(factor(g))}
    stderr1 = function(x) {
        sqrt(var(x, na.rm = F)/sum(!is.na(x)))
    }
    SE = tapply(x, factor(g,levels=l), stderr1)
    err.bp = function(dd, error, two.sided = FALSE, numberStandardErrors) {
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
        for (i in 1:length(dd)) {
            segments(xval[i], dd[i], xval[i], ERR1[i])
            segments(xval[i] - MW, ERR1[i], xval[i] + MW, ERR1[i])
            if (two.sided) {
                segments(xval[i], dd[i], xval[i], ERR2[i])
                segments(xval[i] - MW, ERR2[i], xval[i] + MW,
                         ERR2[i])
            }
            if (is.na(dd[i])) segments(xval[i], -10^10, xval[i], 10^10, lty="dashed", col="grey")
        }
    }
    if (is.na(ylab))
        ylab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(xlab))
        xlab = as.character(match.call(expand.dots = FALSE)$g)
    p1 = signif(kruskal.test(x, factor(g))$p.value, 2)
    Means1 = tapply(x, factor(g,levels=l), mean, na.rm = TRUE)
    p1 = signif(kruskal.test(x ~ factor(g))$p.value, 2)
    if (AnovaTest)
        p1 = signif(anova(lm(x ~ factor(g)))$Pr[[1]], 2)
    if (AnovaTest | KruskalTest)
        main = paste(main, " p=", as.character(p1))
    ret = barplot(Means1, main = main, col = color, xlab = xlab,
                  ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab,
                  cex.main = cex.main, ...)
    abline(h = 0)
    if (numberStandardErrors > 0) {
        err.bp(as.vector(Means1), as.vector(SE), two.sided = two.sided,
               numberStandardErrors = numberStandardErrors)
    }
    attr(ret, "height") = as.vector(Means1)
    invisible(ret)
}
