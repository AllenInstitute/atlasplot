
#' fet_subregions_plot
#' 
#' Plot gene expression values for each structure in each donor for the fetal human brain atlas. This function
#' requires the fetdata package (GITHUB GOES HERE), and is a thin wrapper
#' around a modified version of WGCNA's verboseBarplot
#' 
#' @export
#' @param gene Gene in the `fetdata::datFET` data set; character string
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_fetalHuman_structureBarPlotAll.pdf`
#' @examples
#' # call on a given gene
#' fet_subregions_plot("SHH")
fet_subregions_plot <- function(gene){

    # ensure the user has fetdata installed; inform them to download it if not
    if (!requireNamespace("fetdata", quietly = TRUE)){
        stop("fetdata needed for this function to work. Please install it.",
             call. = FALSE)
    }
    
    # check if fetdata is loaded; if not load it and note loading
    if ("package:fetdata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("fetdata")  # technically breaking the rules
    }
    
    # make sure the gene they're asking for is in the FET
    if (!any(is.element(gene, genesFET))){
        stop(paste(gene,"does not appear to be in the Fetal Human Brain Atlas"))
    }
    
    # set a logical max y-lim for the gene we're looking at
    vlim <- sapply(donor_framesFET, function(x){
        brain <- get(x)
        bool_select <- (gene == brain$gene)
        as.numeric(quantile(brain$value[bool_select], 0.995))
    })
    
    # set ylim and create pdf for plotting  
    ylim <- c(0,2^max(vlim))
    
    # ensures the pdf is closed even after encountering an error
    tryCatch({
        pdf(paste(gene, "fetalHuman_structureBarPlotAll.pdf", sep="_"), height=10,width=50)
        par(mfrow=c(4,1),mar=c(5,10,4,2))
        
        # iterate through the brains and plot the result for each of them
        i <- 1
        for (d in donor_framesFET) {
            brain <- get(d)
            donorID <- as.character(brain$donorID[1])
            bool_select <- (gene == brain$gene)
            .verboseBarplot2(2^brain$value[bool_select],
                             factor(brain$brain_structure[bool_select],levels=subregionsFET),
                             main=paste(gene,"- Brain:",i, '-', donor_to_age[donorID]),las=2,
                             xlab="",ylab="",ylim=ylim, color = colsFET)
            i <- i + 1
        }}, finally = {dev.off()})
    
    # unload fetdata to keep it from affecting global environment
    if (!loaded){
        detach("package:fetdata", unload = TRUE)  
    }
    
}
