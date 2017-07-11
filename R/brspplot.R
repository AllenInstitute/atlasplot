#'
#'
#'
brsp_subregions_plot <- function(gene, group, technique){

    # ensure the use has brspdata installed; if not inform them of the need
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
    groups <- c("gene", "exon")
    techniques <- c("rna", "array")
    if (group %in% groups & technique %in% techniques) {
        data_id <- paste("datBRSP", group, technique, sep=".")
    } else {
        message <- paste("Either group or technique is not a valid selection\n",
                         "\tgroups:", groups, "\n", "\ttechniques:", techniques)
        stop(message, call. = FALSE)
    }
    
    # check that gene is in brspdata
    # TODO

    data <- get(data_id)
    pdf(paste(gene, "BrainSpanRNAseq_brainRegionAndAgeBarplot.pdf",sep="_"),height=25,width=75)
    par(mfrow=c(4,1), mar=c(5,10,4,2))
    # verboseBarplot(dat[orderBS2],factor(samplesBS2,levels=unique(samplesBS2)),color=colorsBS2[isFirst(samplesBS2)],
    #                main=gn,las=2,xlab="",ylab="Expression Level",ylim=ylim)
    # plot.new()
    # verboseBarplot(dat,factor(samplesBS,levels=unique(samplesBS)),color=numbers2colors(ageBS[isFirst(samplesBS)]),
    #                main=gn,las=2,xlab="",ylab="Expression Level",ylim=ylim)
    dev.off()
}