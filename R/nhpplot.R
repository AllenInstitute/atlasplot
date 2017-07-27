
#' nhp_subregions_plot
#' 
#' Plot gene expression values for each structure in each donor for the nhpal human brain atlas. This function
#' requires the nhpdata package (GITHUB GOES HERE), and is a thin wrapper
#' around a modified version of WGCNA's verboseBarplot
#' 
#' @export
#' @param gene Gene in the `nhpdata::datnhp` data set; character string
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_nhpalHuman_structureBarPlotAll.pdf`
#' @examples
#' # call on a given gene
#' nhp_subregions_plot("SHH")
nhp_subregions_plot <- function(gene){

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

    ggplot2::ggsave(plot=p1, file = paste0(f.name, ".pdf"), width=4, height=12) 
}


#'@export
species_expression_time_series<- function(Gene, col_map=rainbow) {
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

    ggplot2::ggsave(g3, file= paste0(Gene,"_comparison_time_series.pdf"), width=10, 
                    height=10)
}