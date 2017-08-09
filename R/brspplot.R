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

#' brsp_subregions_plot
#' 
#' 
#' Plot subregion expression levels for the human brainspan atlas. Provides
#' coloring by both age and ontology in a single plot.
#' 
#' @export
#' @param gene Gene in the `brspdata` package
#' @param technique Either `rna`, for RNAseq, or `array`, for MicroArray
#' @param cmp R color-map function for ages; default heat.colors
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' brsp_subregions_plot("SCARF1", "rna")
brsp_subregions_plot <- function(gene, technique = "array", cmp = heat.colors, 
                                 save_pdf=TRUE) {

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
    
    # download ontology for color and plotting order
    ont_c <- c("acronym", "color_hex_triplet", 'graph_order')
    onto <-  .fetch_ontology("16")[ont_c]
    data <- merge(data, onto, by.x = "brain_structure", by.y = "acronym")
    
    # set some plot universal features
    ymax <- 2^quantile(data$value, 0.995)
    title <- .create_title(gene, group, technique)

    # plot; tryCatch is to ensure proper closing of resources
    tryCatch({
        mar <- c(1,1,1,1)
        if (save_pdf) {
            pdf( title, height=25, width=75)
            mar <- c(5,10,4,2)
        }
        par(mfrow=c(4,1), mar=mar)
        
        # create plot order and change how things are added
        p_ord <- order(data$graph_order)
        data <- data[p_ord,]    
        brain_age <- paste(data$brain_structure, data$age)
        
        # create a color map according to the ontology
        str_col <- paste("#", data$color_hex_triplet, sep="")
        color_map <- .create_str_colormap(brain_age, str_col)
        
        # plot organized by structure
        .verboseBarplot(2^data$value, factor(brain_age, unique(brain_age)), main=gene,
                        las=2,xlab="",ylab="Expression Level", ylim=c(0,ymax),
                        color = color_map)
        plot.new()

        # assign sortable age for graphing; see .assign_value
        order_ages <- order(sapply(brain_age, .assign_value))
        unq_ages <-unique(brain_age[order_ages])
        color_map <- .create_age_colormap(unq_ages, cmp)
  
        # plot organized by age
        .verboseBarplot(2^data$value, factor(brain_age, unq_ages), main=gene,
                    las=2,xlab="",ylab="Expression Level", 
                    ylim=c(0,ymax), color = color_map)
        }, finally = { 
            if (save_pdf) {
                dev.off()
            }
    })

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


.create_str_colormap <- function(ages, colors) {
    # helper function to create color vector for ontology ordering
    brain_color <- unique(cbind(ages, colors))
    brain_color[,2]
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
