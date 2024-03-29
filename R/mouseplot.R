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

#' mouse_subreginos_plot
#'
#'
#' Plot gene expression values for strucutures at a given ontology depth for
#' the mouse or developing mouse atlas. This function pull the mean expression
#' energies directly from the Allen institute API. Utilizes WGCNA's
#' verboseBarplot function for visualization.
#'
#' @export
#' @param gene Gene in the Allen Mouse or Developing Mouse atlas; character string
#' @param atlas Atlas to plot the gene for; either "DevMouse" or "Mouse"
#' @param log_transform Boolean for log vs linear data
#' @param experiment Experiment ID string; Experiment ID for the adult mouse atlas. The devmouse uses all experiments, so the ID is ignored. If not ID is provided, the program will ask you to choose one from a list.
#' @param struct_depth ontology structure depth; default of 3.  It is often useful creating multiple plots with different structure depths for comparison.
#' @param im_width Output image width; optional.  Note that if the plot is too narrow, not all brain structure names will show up.
#' @param im_height Output image height; optional
#' @param save_pdf Save plot to disk as pdf or return to console
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' mouse_structureplot("Scarf1", atlas="Mouse")
#'
#' # call a given gene with a specific experiment
#' mouse_structureplot("Shh", atlas="DevMouse", struct_depth = 4, im_width = 25)
mouse_subregions_plot <- function(gene, atlas, log_transform = FALSE, experiment=NULL,
                                  struct_depth = 3, im_width = NULL,
                                  im_height = NULL, save_pdf=TRUE) {
    # current choices
    atlases <- list(
        "Mouse" = "1",
        "DevMouse" = "3"
    )

    graphs <- list(
        "Mouse" = "1",
        "DevMouse" = "17"
    )

    # ensure an atlas is picked before executing code
    if (missing(atlas)) {
        msg <- "Please supply an atlas to query:\n"
        for (n in names(atlases)) {
            msg <- paste(msg, "\t* ", n, "\n", sep = "")
        }
        stop(msg, call. = FALSE)
    }

    # save old par
    opar <- par(no.readonly = TRUE)

    # check atlas is in the list; if not stop the execution
    if (atlas %in% names(atlases)) {
        atlas_id <- atlases[[atlas]]
        graph_id <- graphs[[atlas]]
    } else {
        stop("Not A Valid Atlas", call. = FALSE)
    }

    # format data frame for the correct mouse call
    if (atlas_id == "1") {
        structure_unionize <- .download_adult_mouse(gene, atlas_id, graph_id,
                                                    experiment)
    } else if (atlas_id == "3") {
        warning("Different ages have different ontology depths")
        structure_unionize <- .download_dev_mouse(gene, atlas_id, graph_id,
                                                  experiment)
    } else {
        stop("An Error Has Occured with Altas ID's", call. = FALSE)
    }

    # get experiment number for a gene; can be multiple experiments per gene
    if (atlas_id == "1") {
        ontology <- .fetch_ontology(graph_id)
    } else if (atlas_id == "3") {
        ontology <- .fetch_ontology("17")
    }

    # ensure struct_depth makes sense
    cond <- struct_depth > min(ontology$depth) & struct_depth <= max(ontology$depth)
    if (!cond) {
        msg <- "Structure depth is too deep or too shallow; default of 3"
        stop(msg, call. = FALSE)
    }

    # select important variables and merge to create a new frame
    ont_c <- c("id", "name", "acronym",
               "color_hex_triplet", "depth", "structure_id_path", "graph_order")
    onto_str <- merge(
        ontology[ont_c],
        structure_unionize,
        by.x = "id", by.y = "structure_id"
    )

    # exlude objects that are not at our structure depth
    include  <- onto_str$depth == struct_depth
    onto_str <- onto_str[include, ]

    if (log_transform) {
        onto_str$expression_energy <- log(onto_str$expression_energy, base = 2)
    }

    if (atlas_id == "1") {
        .plot_mouse_substructures(onto_str, atlas, altas_id, gene, struct_depth,
                                  im_height, im_width, save_pdf)
    } else if (atlas_id == "3") {
        .plot_devmouse_substructures(onto_str, atlas, atlas_id, gene,
                                    struct_depth, im_height, im_width, save_pdf)
    } else {
        stop("An Error Has Occured with Atlas ID's", call. = FALSE)
    }

    # restore par settings
    par(opar)
}
