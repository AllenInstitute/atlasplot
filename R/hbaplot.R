#    Copyright (C) 2017 Allen Institute
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#' hba_subregions_plot
#'
#' Plot gene expression values for each structure in each donor. This function
#' requires the hbadata package (https://github.com/AllenInstitute/hbadata), and
#' is a thin wrapper around a modified version of WGCNA's verboseBarplot
#'
#' @export
#' @param gene Gene in the `hbadata::datHBA` data set; character string
#' @param log_transform Boolean for plotting of linear or log data
#' @param im_height PDF image height
#' @param im_width PDF image width
#' @param save_pdf Output pdf or plot object
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_HBA_subregionsPlot.pdf`
#' @examples
#' # call on a given gene
#' hba_subregions_plot("FOXP2")
hba_subregions_plot <- function(gene, log_transform = FALSE, im_height = 20,
                                im_width = 32, save_pdf=TRUE){

    # ensure the user has hbadata installed; inform them to download it if not
    if (!requireNamespace("hbadata", quietly = TRUE)){
        stop("hbadata needed for this function to work. Please install it.",
             call. = FALSE)
    }

    # check if hbadata is loaded; if not load it and note that
    if ("package:hbadata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("hbadata")  # technically breaking the rules
    }

    # make sure the gene they're asking for is in the HBA
    if (!any(is.element(gene, genesHBA))){
        stop(paste(gene, "does not appear to be in the Human Brain Atlas"))
    }

    # save old par
    opar <- par(no.readonly = TRUE)

    # set a logical max y-lim for the gene we're looking at
    vlim <- sapply(donor_framesHBA, function(x){
        brain <- get(x)
        bool_select <- (gene == brain$gene)
        as.numeric(quantile(brain$value[bool_select], 0.995))
        })

    # set ylim and create pdf for plotting
    if (!log_transform) {
        ylim <- c(0, 2 ^ max(vlim))
    } else {
        ylim <- c(0, max(vlim))
    }

    # tryCatch block to ensure the pdf is closed for unexpected errors
    tryCatch({

        mar <- c(1, 1, 1, 1)
        if (save_pdf) {
            pdf(paste(gene, "HBA_subregionsPlot.pdf", sep = "_"),
                height = im_height, width = im_width)
            mar <- c(5, 10, 4, 2)
        }

        par(mfrow = c(6, 1), mar = mar)

        # iterate through the brains and plot the result for each of them
        i <- 1
        for (d in donor_framesHBA) {
            brain <- get(d)

            if (!log_transform) {
                brain$value <- 2 ^ brain$value
            }

            bool_select <- (gene == brain$gene)
            .verboseBarplot2(brain$value[bool_select],
                            factor(brain$brain_structure[bool_select],
                            levels = subregionsHBA),
                            main = paste(gene, "- Brain:", i), las = 2,
                            xlab = "", ylab = "", ylim = ylim, color = colsHBA)
            i <- i + 1
        }}, finally = {
            if (save_pdf) {
                dev.off()
            }
        })

    # unload hbadata to keep it from affecting global environment
    if (!loaded){
        detach("package:hbadata", unload = TRUE)
    }

    # restore par settings
    suppressWarnings(par(opar))
}
