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

.download_adult_mouse <- function(gene, atlas_id, graph_id, experiment=NULL) {
    # function to download adult mouse
    if (is.null(experiment) & atlas_id == "1") {
        experiment <- .pick_mouse_experiment(gene, atlas_id)
    }

    .fetch_mouse_unionize(experiment, graph_id)
}

#--------------------------------API FUNCTIONS------------------------------------------#
.pick_mouse_experiment <- function(gene, atlas_id) {
    # function for user to pick experiment if one isn't supplied
    # Arguments
    # ---------
    # gene -- gene acronym used in allen institutes for a mouse
    experiments <- .fetch_mouse_experiments(gene, atlas_id)
    if (nrow(experiments) == 1) {
        return(experiments$id)
    }

    print("There are multiple experiments for this gene.")
    print(experiments)

    # loop until user has valid selection
    while (TRUE) {
        experiment <- readline("Please an experiment id -> ")

        if (experiment == "quit") {
            stop("User quit", call. = FALSE)
        }

        # change sting to integer and check that it is valid
        exp_input <- as.integer(experiment)
        if (exp_input %in% experiments$id) {
            return(exp_input)
        } else if (!is.na(exp_input) & exp_input <= nrow(experiments)) {
            print("Please use experiment id *not* row number")
        } else {
            print("NOT A VALID SELECTION")
        }
    }
}


.fetch_mouse_experiments <- function(gene, atlas_id) {
    # retrieves all of the experiments for a given gene in the mouse atlas
    # Arguments
    # ---------
    # gene -- a single gene acronym; spelling and capitalization important

    plane <- list("1" = "coronal", "2" = "sagital")

    # define query and url
    set <- "SectionDataSet"
    # eq right now, il for fuzzy searching or find table fo human mouse homology
    query <- paste("query.json?criteria=products[id$eq", atlas_id,
                   "],genes[acronym$eq'", gene,
                   "'],section_images[failed$eqfalse]", sep = "")

    # pull URL and save locally; result is a bool if failed
    URL <- .construct_api_url(set, query)
    result <- .safe_api_call(URL)

    # see if the api call was succesful; if not stop
    if (is.logical(result)) {
        msg <- paste(gene, "is not a valid gene input")
        stop(msg, call. = FALSE)
    }

    # check if the gene returned a valid message; if not may not exist for mouse
    if (length(result$msg) == 0) {
        msg <- paste("API Call successful but no content delivered for ", gene, ".\nCheck spelling and capitalization; there may be no experiments with ", gene, sep = "")
        stop(msg, call. = FALSE)
    }

    # extract the experiment numbers and thier plane
    experiments <- result$msg[c("id", "plane_of_section_id")]
    experiments$plane <- sapply(experiments$plane_of_section_id,
                                function(x) {
                                    plane[[as.character(x)]]
                                })

    # remove id and return; this format is easy to scan
    experiments[c("id", "plane")]
}


.fetch_mouse_unionize <- function(exp_id, graph_id) {
    # fetch the structure unionization for an experiment
    # Arguments
    # ---------
    # exp_id -- id for an experiment with a gene; from fetch_mouse_experiments
    set <- "StructureUnionize"
    query <- paste("query.json?criteria=[section_data_set_id$eq", exp_id,
                  "],structure[graph_id$eq", graph_id, "]", sep = "")
    URL <- .construct_api_url(set, query)

    result <- .safe_api_call(URL)

    # check that experiment ID is valid
    if (is.logical(result)) {
        msg <- paste(exp_id, "is not a valid experiment ID")
        stop(msg, call. = FALSE)
    }

    # ensure there is content in the message
    if (length(result$msg) == 0) {
        msg <- paste("API Call successful but no content delivered for", exp_id,
                    ".\nCheck to ensure this is a valid experiment")
        stop(msg, call. = FALSE)
    }
    result$msg[c("structure_id", "expression_energy")]
}


#------------------------------PLOT FUNCTIONS-------------------------------------------#
.plot_mouse_substructures <- function(onto_str, atlas, atlas_id, gene,
                                      struct_depth, im_height, im_width,
                                      save_pdf) {
    # function to plot the mouse atals substructures
    ylim <- c(0, quantile(onto_str$expression_energy, 0.995))
    p_ord <- order(onto_str$graph_order)

    if (is.null(im_height)) {
        im_height <- 9
    }

    if (is.null(im_width)) {
        im_width <- 20
    }

    # create pdf; tryCatch to ensure proper resource closing
    tryCatch({

        mar <- c(1, 1, 1, 1)
        if (save_pdf) {
                pdf(paste(gene, atlas, struct_depth,
                    "brainRegionBarplot_ExpressionEnergies.pdf", sep = "_"),
            height = im_height, width = im_width)
            mar <- c(5, 10, 4, 2)
        }
        # save old par settings
        par(mfrow = c(2, 1), mar = mar)

        # plot the structures
        .verboseBarplot2(onto_str$expression_energy[p_ord],
                        factor(onto_str$name[p_ord],
                        unique(onto_str$name[p_ord])), main = gene,
                        color = paste("#", onto_str$color_hex_triplet,
                        sep = ""), las = 2, xlab = "",
                        ylab = "Mean Expression Energy", ylim = ylim,
                        KruskalTest = FALSE)
    }, finally = {
        if (save_pdf) {
            dev.off()
        }
    })
}
