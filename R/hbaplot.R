
#' hba_subregions_plot
#' 
#' Plot gene expression values for each structure in each donor. This function
#' requires the hbadata package (GITHUB GOES HERE), and is a thin wrapper
#' around a modified version of WGCNA's verboseBarplot
#' 
#' @export
#' @param gene Gene in the `hbadata::datHBA` data set; character string
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_HBA_subregionsPlot.pdf`
#' @examples
#' # call on a given gene
#' hba_subregions_plot("FOXP2")
hba_subregions_plot <- function(gene){

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
        stop(paste(gene,"does not appear to be in the Human Brain Atlas"))
    }
    
    # set a logical max y-lim for the gene we're looking at
    vlim <- sapply(donor_framesHBA, function(x){
        brain <- get(x)
        bool_select <- (gene == brain$gene)
        as.numeric(quantile(brain$value[bool_select], 0.995))
        })
    
    # set ylim and create pdf for plotting    
    ylim <- c(0,2^max(vlim))
    pdf(paste(gene, "HBA_subregionsPlot.pdf", sep="_"), height=20, width=32)
    par(mfrow=c(6,1), mar=c(5,10,4,2))

    # iterate through the brains and plot the result for each of them
    i <- 1
    for (d in donor_framesHBA) {
        brain <- get(d)
        bool_select <- (gene == brain$gene)
        .verboseBarplot2(2^brain$value[bool_select],
                         factor(brain$brain_structure[bool_select],levels=subregionsHBA),
                         main=paste(gene,"- Brain:",i),las=2,
                         xlab="",ylab="",ylim=ylim, color=colsHBA)
        i <- i + 1
    }
    dev.off()

    # unload hbadata to keep it from affecting global environment
    if (!loaded){
        detach("package:hbadata", unload = TRUE)  
    }
}
