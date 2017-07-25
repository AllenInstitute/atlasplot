
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