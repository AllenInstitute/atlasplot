#' brsp_subregions_plot
#' 
#' 
#' Plot subregion expression levels for the human brainspan atlas. Provides
#' coloring by both age and ontology
#' 
#' @export
#' @param gene Gene in the `brspdata` package
#' @param technique Either `rna`, for RNAseq, or `array`, for MicroArray
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' brsp_subregions_plot("SCARF1", "rna")
brsp_subregions_plot <- function(gene, technique){

    # ensure the user has brspdata installed; if not inform them of the need
    if (!requireNamespace("brspdata", quietly = TRUE)){
        stop("brspdata needed for this function to work. Please install it",
             call. = FALSE)
    }

    # check if brspdata is loaded; if not load it and note load status
    if ("package:brspdata" %in% search()){
        loaded <- TRUE
    } else {
        loaded <- FALSE
        library("brspdata") # this is cheating
    }

    # get the correct data file and ensure that the user chooses a real option
    group <- "gene"
    techniques <- c("rna", "array")
    if (technique %in% techniques) {
        data_id <- paste("datBRSP", group, technique, sep=".")
    } else {
        message <- paste("Technique is not a valid selection\n","techniques:")
        for (t in techniques) {
            message <- paste(message, t)
        }
        stop(message, call. = FALSE)
    }

    data <- get(data_id)

    # ensure that the gene is a valid selection
    if(!any(is.element(gene, data$gene))) {
        stop(paste(gene, "does not appear to be in", data_id))
    }

    # reduce data down to the necessary elements, easier for plotting later
    bool_select <- data$gene == gene
    data <- data[bool_select,]
    ymax <- 2^quantile(data$value, 0.995)

    title <- .create_title(gene, group, technique)

    tryCatch({
        pdf( title, height=25, width=75)
        par(mfrow=c(4,1), mar=c(5,10,4,2))

        brain_age <- paste(data$brain_structure, data$age)

        # ontology for developing human brain
        # http://api.brain-map.org/api/v2/data/Structure/query.json?criteria=[graph_id$eq16],[acronym$eqA1C]
        
        # plot organized by structure
        .verboseBarplot(2^data$value, brain_age, main=gene,
                        las=2,xlab="",ylab="Expression Level", ylim=c(0,ymax))
        plot.new()

        # assign sortable age for graphing; see .assign_value
        order_ages <- order(sapply(brain_age, .assign_value))
        unq_ages <-unique(brain_age[order_ages])
        color_map <- .create_age_colormap(unq_ages, heat.colors)
  
        # plot organized by age
        .verboseBarplot(2^data$value, factor(brain_age, unq_ages), main=gene,
                    las=2,xlab="",ylab="Expression Level", 
                    ylim=c(0,ymax), color = color_map)
        }, finally = { dev.off() })

    # unload brspdata to keep it from affecting global environment
    if (!loaded){
        detach("package:brspdata", unload = TRUE)
    }
}


#------------------------------HELPER FUNCTIONS--------------------------------#

.create_title <- function(gene, group, technique){
    # title helper function; formats the title for the technique and gene used
    if (technique == "rna") {
        type <- paste(group, "RNAseq", sep = "")
    } else {
        type <- paste(group, "Microarray", sep = "")
    }
    
    return(paste(gene, "BrainSpan", type, "brainRegionaAndAgeBarplot.pdf", sep = "_"))
}


.create_age_colormap <- function(ages, cmap){
    # creates colors for the age plot
    cols <- cmap(length(ages))
    col_list <- list()
    for (i in 1:length(ages)) {
        col_list[ages[i]] <- cols[i]
    }
    
    colors <- c()
    for (a in ages) {
        colors <- c(colors, col_list[[a]])
    }

    colors
}


.assign_value <- function(astring) {
    # helper function to assign sortable integer for age classes
    # transfers age strign into sortable integer
    state_to_id <- list("yrs" = 1000, "mos" = 100, "pcw" = 0)

    # ensure is a string and remove the second part
    astring <- as.character(astring)
    astring <- strsplit(astring, ' ')[[1]][[2]]
    
    # split on '_' to get the age code and the number of that unit
    splt_str <- strsplit(astring, '_')[[1]]
    age_value <- as.integer(splt_str[[1]])
    unit_code <- state_to_id[[splt_str[[2]]]]
    age_value + unit_code
}
