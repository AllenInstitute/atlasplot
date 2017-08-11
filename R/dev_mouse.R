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

.download_dev_mouse <- function(gene, atlas_id, graph_id, experiment) {
    if (!is.null(experiment)) {
        msg <- "Experiment IDs are not used to search the Developing Mouse Atlas and as such will be ignored"
        warning(msg)
    }

    .fetch_dev_mouse_unionize(gene, graph_id)
}

#--------------------------------API FUNCTIONS------------------------------------------#
.fetch_dev_mouse_unionize <- function(gene, graph_id) {
    set <- "StructureUnionize"
#    query <- paste("query.json?criteria=rma::criteria,section_data_set[delegate$eqfalse](genes[acronym$eq'",
#                  gene, "'],specimen(donor(age[name$in'E11.5','E13.5','E15.5','E18.5','P4','P14','P28']))),structure[graph_id$eq17]&include=section_data_set(specimen(donor(age)))",
#                  sep = "")
    query <- paste("query.json?criteria=section_data_set[delegate$eqtrue](genes[acronym$eq'",
                  gene, "']),structure[graph_id$in'4','13','17']&include=section_data_set(specimen(donor(age)))",
                  sep = "")
    
    URL <- .construct_api_url(set, query)

    result <- .safe_api_call(URL)

    # check that experiment ID is valid
    if (is.logical(result)) {
        msg <- paste("This may not be a valid search")
        stop(msg, call. = FALSE)
    }

    # ensure there is content in the message
    if (length(result$msg) == 0) {
        msg <- paste("API Call successful but no content delivered for", gene,
                     ".\nCheck to ensure this is a valid gene in the DevMouse atlas")
        stop(msg, call. = FALSE)
    }
    ages <- result$msg$section_data_set$specimen$donor$age$name
#    saveRDS(ages, "age.rda")
#    stop()
    cbind( result$msg[c("structure_id", "expression_energy")], ages )
}


#-------------------------------PLOT FUNCTIONS------------------------------------------#
.plot_devmouse_substructures <- function(onto_str, atlas, atlas_id, gene,
                                        struct_depth, im_height, im_width,
                                        save_pdf) {

    ages <- c("E11.5", "E13.5", "E15.5", "E18.5", "P4", "P14", "P28")

    # function to plot the mouse atals substructures
    ylim <- c(0, quantile(onto_str$expression_energy, 0.995))

     if (is.null(im_height)) {
        im_height <- 20
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
        par(mfrow = c(8, 1), mar = mar)
        # plot the structures
        for (age in ages) {
            i_age <- onto_str$ages == age
            onto_str_i <- onto_str[i_age, ]
            if (length(unique(onto_str_i$ages)) <= 1) {
                next
            }
            p_ord <- order(onto_str_i$graph_order)
            .verboseBarplot2(onto_str$expression_energy[p_ord],
                            factor(onto_str$name[p_ord],
                            unique(onto_str$name[p_ord])), main = paste(gene,
                            age, sep = " - "), color = paste("#",
                            onto_str$color_hex_triplet, sep = ""), las = 2,
                            xlab = "", ylab = "Mean Expression Energy",
                            ylim = ylim, KruskalTest = FALSE)
        }
    }, finally = {
        if (save_pdf) {
            dev.off()
        }
    })
}
