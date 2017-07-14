#'
#'
#'
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

        # plot organized by structure
        .verboseBarplot(2^data$value, brain_age, main=gene,
                        las=2,xlab="",ylab="Expression Level", ylim=c(0,ymax))
        plot.new()

        # assign sortable age for graphing; see .assign_value
        order_ages <- order(sapply(brain_age, .assign_value))

  
        # plot organized by age
        .verboseBarplot(2^data$value, factor(brain_age, unique(brain_age[order_ages])),
                    main=gene,las=2,xlab="",ylab="Expression Level", 
                    ylim=c(0,ymax))
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


.create_colormap <- function(ages, cmap){
    split_ages <- strsplit(gsub("([0-9]+) ", "\\1_", ages), " ")
    age_id <- c()
    for (age in split_ages){
        age_id <- c(age_id, age[[2]])
    }
    age_id <- unique(age_id)
    colormap <- list()
    colors <- cmap(length(age_id))
    i <- 1
    for (age in age_id) {
        colormap[age] = colors[i]
        i <- i + 1
    }

    colormap
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
