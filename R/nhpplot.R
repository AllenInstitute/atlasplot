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

#' nhp_cortex_series_plot
#'
#' Plot times series for a single gene in the non-human primate brain atlas from the Allen
#' institute. This function relies on the nhpdata package, located at (GITHUB GOES HERE).
#' 
#' @export
#' @param gene Gene in the `nhpdata::exprl2` data set; character string
#' @param save_pdf Save to disk or return to console
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' nhp_cortex_series_plot("MKI67")
#'
nhp_cortex_series_plot <- function(gene, save_pdf=TRUE){

    # ensure the user has nhpdata installed; inform them to download it if not
    if (!requireNamespace("nhpdata", quietly = TRUE)){
        stop("nhpdata needed for this function to work. Please install it.",
             call. = FALSE)
    }
    
    # check if ggplot2 is installed
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        msg <- "ggplot2 required for this function to work\n Install with `install.packages(\"ggplot2\")`"
        stop(msg, call.=FALSE)
    }
    
    # check if nhpdata is loaded; if not load it and note loading
    if ("package:nhpdata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("nhpdata")  # technically breaking the rules
    }
    
    # get plotting parameters
    num.timepoints <- length(table(nhpdata::exprl2$age))
    pd <- ggplot2::position_dodge(width = 0.4)
    
    i <- match(gene, nhpdata::probes$macaque_genesymbol)
    if (is.na(i)) {
        # confirms a match is found
        msg <- "No matches for the given gene; Check spelling and capitalization" 
        stop(msg, call.=FALSE)
    }
    
    # Add human ortholog symbol
    m.gene <- nhpdata::probes$macaque_genesymbol[i]
    m.gene.title <- ifelse(m.gene == "LOC713917", "NGEF (LOC713917)", m.gene)
    m.probe <- nhpdata::probes$probeid[i]
    m.id <- nhpdata::probes$macaque_entrezid[i]
    m.ortho <- nhpdata::probes$human_genesymbol[i]
    if (is.na(m.ortho)) {
        m.ortho <- "_"
    }

    # Select gene data
    exprl2.subset <- subset(nhpdata::exprl2, gene == m.gene)
    exprl2.subset <- droplevels(exprl2.subset)

    f.name <- paste(m.gene, m.ortho, m.id, m.probe, sep="_")

    suppressWarnings(
        p1 <- ggplot2::ggplot(data = exprl2.subset, 
               ggplot2::aes(x = as.numeric(age), y = expr, 
                   color=area, fill=area, shape=area)) + 
    ggplot2::facet_grid(layer ~ .) + 
    ggplot2::stat_summary(fun.y = mean, geom = "line", position=pd) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_sdl, geom = "pointrange", mult=1, position=pd) +
    ggplot2::xlab("Age") + 
    ggplot2::ylab("log2(Expression)") +
    ggplot2::ggtitle(m.gene.title) +
    ggplot2::scale_x_continuous(breaks = 1:num.timepoints, 
                     labels = levels(exprl2$age)) +
    ggplot2::scale_color_manual(values=c("#f8766d", "#00ba38", "#619cff")) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45),
                       axis.text = ggplot2::element_text(size=9),
                       panel.grid.minor = ggplot2::element_blank())
    )

    if (save_pdf) {
        ggplot2::ggsave(plot=p1, file = paste0(f.name, ".pdf"), width=4, height=12)
    }

    if (!loaded){
        detach("package:nhpdata", unload = TRUE)  
    }
    p1
}



#' species_expression_time_series
#' 
#' Plot gene expression values for a time series of ages for human, non-human primate, mouse,
#' and rat. Will only include data that is present in plots. This function requires the
#' nhpdata package (GITHUB GOES HERE)/
#' 
#' @export
#' @param gene Gene in the `nhpdata::dev.expr2` data set; character string
#' @param colormap Color pallete used to choose plot colors; default rainbow
#' @param save_pdf Save to disk or return to console
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' species_expression_time_series("EMX2", colormap = heat.colors)
#'
species_expression_time_series<- function(Gene, col_map=rainbow, save_pdf=TRUE) {
     # ensure the user has nhpdata installed; inform them to download it if not
    if (!requireNamespace("nhpdata", quietly = TRUE)){
        stop("nhpdata needed for this function to work. Please install it.",
             call. = FALSE)
    }

    # check if ggplot2 is installed
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        msg <- "ggplot2 required for this function to work\n Install with `install.packages(\"ggplot2\")`"
        stop(msg, call.=FALSE)
    }

    # check if nhpdata is loaded; if not load it and note loading
    if ("package:nhpdata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("nhpdata")  # technically breaking the rules
    }

    # subset to the appropriate gene
    dev.expr2.subset <- subset(dev.expr2, gene == Gene)

    # Remove redundant species from plots
    dev.expr2.subset <- subset(dev.expr2.subset, ! (gene == "LGALS1" & species == "mouse") &
                                 ! (gene == "EMX2" & species == "human_bc") & 
                                 ! (gene == "CNTN1" & species == "human_bc") &
                                 ! (gene == "BMP3" & species == "human_bc") & 
                                 ! (gene == "CNTN2" & species == "human_bc"))
    dev.expr2.subset <- droplevels(dev.expr2.subset)

    # Reorder levels to change panel arrangement
    dev.expr2.subset$gene <- factor(dev.expr2.subset$gene, 
                                    levels=c("EMX2", "BMP3", "LGALS1",
                                             "CNTN1", "CNTN2", "LIN7A"))

    pal1 <- col_map(5)[c(1, 5, 2:4)]

    g3 <- ggplot2::ggplot(dev.expr2.subset, ggplot2::aes(x=escore, y=exprz, shape=species,
                                       col=species, fill=species)) +
        ggplot2::geom_point(alpha=0.5) +
        ggplot2::geom_smooth(method="loess", se=FALSE, span=1, size=2) +
        #   #   geom_point() + geom_line() +  # More clear legend
        ggplot2::facet_wrap( ~ gene, ncol=3, scales="free_y") +
        ggplot2:: theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_color_manual(values = pal1) +
        ggplot2::scale_fill_manual(values = pal1)

    if (save_pdf) {
        ggplot2::ggsave(g3, file= paste0(Gene,"_comparison_time_series.pdf"), width=10, 
                    height=10)
    }

    if (!loaded){
        detach("package:nhpdata", unload = TRUE)  
    }
    g3
}


#' nhp_cortex_expression2D_plot
#'
#' Creates a boxed heatmap of the expression values for a single gene in the non-human
#' primate atlas. 
#' 
#' @export
#' @param gene Gene in the `nhpdata::datNHP` data set; character string
#' @param colVec List of colors to create linear color map; c("white", "red") default
#' @param save_pdf Save to disk or return to console
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' nhp_cortex_expression2D_plot("PAX6", c("green", "skyblue"))
#'
nhp_cortex_expression2D_plot <- function(gene, colVec=c("white", "red"), save_pdf=TRUE) {
    
    if (!requireNamespace("nhpdata", quietly = TRUE)){
        stop("nhpdata needed for this function to work. Please install it.",
             call. = FALSE)
    }
    
    # check if nhpdata is loaded; if not load it and note loading
    if ("package:nhpdata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("nhpdata")  # technically breaking the rules
    }
    
    # ensure gene is in datNHP data set; keeps from having weird messages
    if (!any(gene == datNHP$gene)) {
        msg <- paste(gene, "not found in datNHP; check spelling and capitalization")
        stop(msg, call.=FALSE)
    }
    
    # select the rows for the given gene
    geneDat <- datNHP[datNHP$gene == gene,]
    
    # log the expression and create a named vector
    expr_value <- log(geneDat$value,b=2)
    names(expr_value) <- geneDat$id_string
    
    # convert layer, subregions, and age to character
    age <- as.character(geneDat$age)
    layer <- as.character(geneDat$layer)
    subregions <- as.character(geneDat$subregion)
    
    # create label name
    label_id <- strsplit(as.character(geneDat$id_string[[1]]), "_")[[1]][[3]]
    
    
    # tryCatch to ensure proper resource handling
    tryCatch({
        fn <- paste(gene, "_NHP_2Dexpression.pdf", sep="")
        if (save_pdf) {
            pdf(fn, width=18, height=9)
        }
        .plotMacaqueCortex(expr_value, NULL, layer, subregions, age, layerPositions, 
                           regionPositions, ageOffsets, paste(gene,label_id, sep=" - "), 
                           quantileScale=c(0.1,0.95), colVec=colVec)},
        finally = { if (save_pdf) { dev.off() }
    })
    
    if (!loaded){
        detach("package:nhpdata", unload = TRUE)  
    }
}


#' nhp_cortex_expression2D_small_plot
#'
#' Creates a boxed heatmap of the expression values for a single gene in the non-human
#' primate atlas. This plot is limited to only the V1 region.
#' 
#' @export
#' @param gene Gene in the `nhpdata::datNHP` data set; character string
#' @param colVec List of colors to create linear color map; c("white", "red") default
#' @param save_pdf Save to disk or return to console
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' nhp_cortex_expression2D_plot("PAX6", c("green", "skyblue"))
#'
nhp_cortex_expression2D_small_plot<- function(gene, colVec=c("white", "red"), save_pdf=TRUE) {
    
    if (!requireNamespace("nhpdata", quietly = TRUE)){
        stop("nhpdata needed for this function to work. Please install it.",
             call. = FALSE)
    }
    
    # check if nhpdata is loaded; if not load it and note loading
    if ("package:nhpdata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("nhpdata")  # technically breaking the rules
    }
    
    # ensure gene is in the datNHP data set; keeps from having weird error messages
    if (!any(gene == datNHP$gene)) {
        msg <- paste(gene, "not found in datNHP; check spelling and capitalization")
        stop(msg, call.=FALSE)
    }

    # select the rows for the given gene that are in the V1 subregion
    geneDat <- datNHP[datNHP$gene == gene & datNHP$subregion == "V1",]
        
    # isolate the expression value and create a named vector
    expr_value <- log(geneDat$value, b=2)
    names(expr_value) <- geneDat$id_string
    
    # convert layer, subregions, and age to character
    age <- as.character(geneDat$age)
    layer <- as.character(geneDat$layer)
    subregions <- as.character(geneDat$subregion)
    
    # change layers to be correct for plotting
    layer[layer == "L2-3"] <- "CPo"
    layer[is.element(layer,c("L4Cb","L4Ca","L4B","L4A"))] <- "L4"
    
    
    # create id
    label_id <- strsplit(as.character(geneDat$id_string[[1]]), "_")[[1]][[3]]
    
    # tryCatch to ensure proper resource handling
    tryCatch({
        fn <- paste(gene,"_NHP_small_expression2D_cortex.pdf", sep="")
        legendPos  = c(8,-11.5,-8.5)
        if (save_pdf) {
            pdf(fn,width=7,height=7)
        }
        .plotMacaqueCortexSmall(expr_value, layer, age, layerPositionsS, agePositionsS,
                            paste(gene, label_id, sep=" - "), isLog2=TRUE,
                            combineFn=".meanNA", quantileScale=c(0.1,0.95),
                            linearOrLog="linear", bgPar="grey95", displayLayers=FALSE,
                            legendPos=legendPos, colVec=colVec)
        abline(v=c(0,6,10),lwd=2)
        abline(h=c(0,-8,-12),lwd=2)
    },
            finally={ if (save_pdf) {dev.off()}
    })

    if (!loaded){
        detach("package:nhpdata", unload = TRUE)  
    }
}
