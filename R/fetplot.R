
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


#' fet_expression_2D_plot
#'
#' Plot a 2D heatmap for cortex layers and structures in each donor in the fetal human brain atlas. This function
#' requires the fetdata package, available at (GITHUB GOES HERE).
#' 
#' @export
#' @param gene Gene in the `fetdata::datFET` data set; character string
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_fetalHuman_2DPlotInNeocortex.pdf.pdf`
#' @examples
#' # call on a given gene
#' fet_expression2D_plot("SHH")
#'
fet_expression2D_plot <- function(gene, colbox="red") {

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
    
    # major lobes of the brain; used in subsetting
    lobes <- c("f", "o", "p", "t")
    ontology <- .fetch_ontology("16")  # graph_id 16 is devhuman
    onto_c <- c("acronym", "graph_order")
    # for each brain create the plot; tryCatch for safe resource allocation
    tryCatch({
        pdf(paste(gene,"fetalHuman_2DPlotInNeocortex.pdf",sep="_"),height=10,width=20)
        par(mfrow=c(2,2))
        
        for (d in donor_framesFET) {
            brain <- get(d)
            brain <- merge(brain, ontology[onto_c],
                          by.x="brain_structure", by.y="acronym")
            brain <- brain[order(brain$graph_order),]
            donorID <- as.character(brain$donorID[1])
            
            # title using age and brainID
            title <- paste(gene," - BrainID", donorID, " - ", donor_to_age[donorID])
            
            # select the appropriate fields to plot and create plot regions
            bool_select <- (gene==brain$gene & substr(brain$brain_structure,1,1) %in% lobes)
            brain <- brain[bool_select,]
            lls <- .format_group(brain$brain_structure) # lls for lobe_layer_str

            # remove some layers for comparibility
            lls_keep <- !(lls[,2] %in% c("CP", "SZ"))
            brain <- brain[lls_keep,]
            lls <- lls[lls_keep,]
            
            # remove region for comparibility
            lls_keep <- !(lls[,3] == 'Z')
            brain <- brain[lls_keep,]
            lls <- lls[lls_keep,]

            .plotExpressionMap2D(2^brain$value, factor(lls[,1], levels=c("f","p","t","o")),
                                 lls[,2], main=title, minIs0=TRUE, pch=22,
                                 sampleLabel=lls[,3], bgPar="lightgrey", 
                                 sizeRange=c(3,5), sizeText=1.5, sizeLabel=0.8,
                                 colBox=colbox)
        }
    }, finally = {dev.off()})

    if (!loaded){
        detach("package:fetdata", unload = TRUE)  
    }
}



#------------------------------HELPER FUNCTIONS-----------------------------------------#
.get_all_structs <- function() {
    # debugging helper function
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
    
    structs <- c()
    for (d in donor_framesFET) {
        brain <- get(d)
        structs <- c(structs, unique(brain$brain_structure))
    }
    unique(structs)
}

.format_group <- function(structures) {

    # returns the lobes
    layers_str_lobes <- sapply(structures, function(x) {
        lobe <- substr(x,1,1)
        layer <- substr(x,2,3)
        c_last <- substr(x,nchar(x), nchar(x))

        if (c_last %in% c("i", "o")) {
            layer <- paste(layer, c_last, sep="")
            substruc <- substr(x, 4,nchar(x)-1)
        } else {
            substruc <- substr(x, 4, nchar(x))
        }
        
        if (substruc %in% c("dm", "mi")) {
            substruc <- paste(substruc,lobe,sep="-")
        }
        
        if (layer %in% c("SGi", "SPi", "MZi")) {
            layer <- substr(layer, 1,2)
        }
        c(lobe, layer, substruc)
    })
    t(layers_str_lobes)
}