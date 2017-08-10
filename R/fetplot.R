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

#' fet_subregions_plot
#' 
#' Plot gene expression values for each structure in each donor for the fetal human brain atlas. This function
#' requires the fetdata package (GITHUB GOES HERE), and is a thin wrapper
#' around a modified version of WGCNA's verboseBarplot
#' 
#' @export
#' @param gene Gene in the `fetdata::datFET` data set; character string
#' @param save_pdf Save to disk or output to console
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_fetalHuman_structureBarPlotAll.pdf`
#' @examples
#' # call on a given gene
#' fet_subregions_plot("SHH")
fet_subregions_plot <- function(gene, save_pdf=TRUE){

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

    # save old par setting to reset after 
    opar <- par(no.readonly = TRUE) 

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

        mar <- c(1,1,1,1)
        if (save_pdf) {
            pdf(paste(gene, "fetalHuman_structureBarPlotAll.pdf", sep="_"), height=10,width=50)
            mar <- c(5,10,4,2)
        }
        par(mfrow=c(4,1), mar=mar)

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
        }}, finally = {
            if (save_pdf) {
                dev.off()
            }
        })
    
    # unload fetdata to keep it from affecting global environment
    if (!loaded){
        detach("package:fetdata", unload = TRUE)  
    }
    
    # restore par settings
    suppressWarnings(par(opar))
}


#' fet_expression_2D_plot
#'
#' Plot a 2D heatmap for cortex layers and structures in each donor in the fetal human brain atlas. This function
#' requires the fetdata package, available at (GITHUB GOES HERE).
#' 
#' @export
#' @param gene Gene in the `fetdata::datFET` data set; character string
#' @param save_pdf Save to disk or output to console
#' @param colbox Heatmap box color; default red
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_fetalHuman_2DPlotInNeocortex.pdf.pdf`
#' @examples
#' # call on a given gene
#' fet_expression2D_plot("SHH")
#'
fet_expression2D_plot <- function(gene, colbox="red", save_pdf=TRUE) {

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
    
    # record old par
    opar <- par(no.readonly = TRUE)
    
    # major lobes of the brain; used in subsetting
    lobes <- c("f", "o", "p", "t")
    ontology <- .fetch_ontology("16")  # graph_id 16 is devhuman
    onto_c <- c("acronym", "graph_order")
    # for each brain create the plot; tryCatch for safe resource allocation
    tryCatch({
        
        if (save_pdf) {
            pdf(paste(gene,"fetalHuman_2DPlotInNeocortex.pdf",sep="_"),height=10,width=20)
        }
        par(mfrow=c(2,2))
        
        for (d in donor_framesFET) {
            brain <- get(d)
            brain_c <- brain$gene == gene  # subsetting before merge gives speed boost
            brain <- merge(brain[brain_c,], ontology[onto_c],
                          by.x="brain_structure", by.y="acronym")
            brain <- brain[order(brain$graph_order),]
            donorID <- as.character(brain$donorID[1])
            
            # title using age and brainID
            title <- paste(gene," - BrainID", donorID, " - ", donor_to_age[donorID])
            
            # select the appropriate fields to plot and create plot regions
            bool_select <- (substr(brain$brain_structure,1,1) %in% lobes)
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
    }, finally = {
        if (save_pdf) {
            dev.off()
        } else {
            print("Recommended saving dimensions: 1700 x 1000")
        }
    })

    if (!loaded){
        detach("package:fetdata", unload = TRUE)  
    }
    
    # restore par settings
    suppressWarnings(par(opar))
}


#' fet_brain_location_expression_plot
#'
#' Creates a matrix of brain images. Each column corresponds to a brain in the fetal human
#' brain atlas, while each row is a specific brain layer.
#'
#' @export
#' @param gene Gene in the `fetdata::datFET` data set; character string
#' @return Creates a new plot of the gene in the current working directory; label as `GENE_fetalHuman_BrainLocationExpressionPlotsForEachLayer.pdf`
#' @examples
#' # call on a given gene
#' fet_expression2D_plot("SHH")
#'
fet_brain_location_expression_plot <- function(gene) {

    # ensure the user has fetdata installed; inform them to download it if not
    if (!requireNamespace("fetdata", quietly = TRUE)) {
        stop("fetdata needed for this function to work. Please install it.",
             call. = FALSE)
    }

    # check that the pixmap requirement is satistified
    if (!requireNamespace("pixmap", quietly = TRUE)) {
        msg <- "pixmap is required for this function to work; Please install it with `install.packages(\"pixmap\")`"
        stop(msg, call. = FALSE)
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
    
    # record old par
    opar <- par(no.readonly = TRUE)

    lobes <- c("f", "o", "p", "t")
    save_pdf <- TRUE  # incase this can be changed later
    tryCatch({
        
        mar <- c(2,2,2,2)
        if (save_pdf) {
            pdf(paste(gene,"fetalHuman_BrainLocationExpressionPlotsForEachLayer.pdf",sep="_")
                ,width=20,height=30,version= "1.4")
            mar <- rep(2,4)
        }
        par(cex=1.5,mar = mar,las = 1)
        op <- par(mar = rep(0, 4)) # we want it to fill the file
        suppressWarnings( x <- .construct_image() )
        pixmap::plot(x)
        lays = c("MZ","CPo","CPi","SP","IZ","SZo","SZi","VZ")
        par(mfrow=c(8,4), new=TRUE)

        i <- 0  # i gives column for the brain
        for (d in donor_framesFET) {
            i <- i + 1  # increment i to get 1 indexed column position
            brain <- get(d)
            donorID <- as.character(brain$donorID[1])

            # select the appropriate fields to plot and create plot regions
            bool_select <- (
                gene==brain$gene & substr(brain$brain_structure,1,1) %in% lobes
            )
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

            j <- 0  # j will give the correct row for the brain
            for (lay in  lays) {
                j <- j + 1 # increment j to get 1 indexed row position

                # keep the layers that match lay
                lls_keep <- lls[,2] == lay
                brain_lay <- brain[lls_keep,]  # keep the full data frame from overwrite
                lls_lay <- lls[lls_keep,]  # see above

                # get position of layer
                par(cex=0.5,mar = rep(2, 4),las = 1,mfg=c(j,i),new=TRUE)
                m <- paste(gene,"-", lay," - BrainID", donorID, " - ", donor_to_age[donorID])

                # select colors and positions for each of the subregions
                col <- yzcFET[lls_lay[,3],4]
                y <- yzcFET[lls_lay[,3], 2]
                z <- yzcFET[lls_lay[,3],3]

                .plotExpressionCoordinates2Db(2^brain_lay$value, z, y, lls_lay[,3], textCol=col,
                                            main=m, bgPar="white", xlim=c(-38,30),
                                            ylim=c(-25,20), minIs0=FALSE)
            }
        }
    }, finally = {
        if (save_pdf) {
            dev.off() 
        }
    })

    if (!loaded){
        detach("package:fetdata", unload = TRUE)  
    }

    # restore par settings
    suppressWarnings(par(opar))
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


# the following two functions reconstruct the image for the brain plots
.construct_image <- function() {    
    red.matrix <- .make_cmatrix(im_comp$red)
    green.matrix <- .make_cmatrix(im_comp$green)
    blue.matrix <- .make_cmatrix(im_comp$blue)
   pixmap::pixmapRGB(c(red.matrix, green.matrix, blue.matrix), nrow=900, ncol=593)
}


.make_cmatrix <- function(clist) {
    # undoes the svd decomp
    clist$u %*% diag(clist$d) %*% t(clist$v)
}