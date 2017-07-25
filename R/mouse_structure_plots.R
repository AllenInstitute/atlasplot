#' mouse_structureplot
#' 
#' 
#' Plot gene expression values for strucutures at a given ontology depth for
#' the mouse or developing mouse atlas. This function pull the mean expression
#' energies directly from the Allen institute API. Utilizes WGCNAs 
#' verboseBarplot function for visualization.
#' 
#' @export
#' @param gene Gene in the Allen Mouse or Developing Mouse atlas; character string
#' @param atlas Atlas to plot the gene for; either "DevMouse" or "Mouse"
#' @param experiment Experiment ID string; Experiment ID for the adult mouse atlas. The devmouse uses all experiments, so the ID is ignored. If not ID is provided, the program will ask you to choose one from a list.
#' @param struct_depth ontology structure depth; default of 3
#' @param im_width Output image width; optional
#' @param im_height Output image height; optional
#' @return Creates a new plot of the gene in the current working directory
#' @examples
#' # call on a given gene
#' mouse_structureplot("Scarf1", atlas="Mouse")
#' 
#' # call a given gene with a specific experiment
#' mouse_structureplot("Shh", atlas="DevMouse", struct_depth = 4, im_width = 25)
mouse_structureplot <- function(gene, atlas, experiment=NULL, struct_depth = 3, 
                                im_width = NULL, im_height = NULL) {
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
            msg <- paste(msg,"\t* ", n, "\n", sep="")   
        }
        stop(msg, call. = FALSE)
    }
    
    # check atlas is in the list; if not stop the execution
    if (atlas %in% names(atlases)) {
        atlas_id <- atlases[[atlas]]
        graph_id <- graphs[[atlas]]
    } else {
        stop("Not A Valid Atlas", call. = FALSE)
    }
    
    # format data frame for the correct mouse call
    if (atlas_id == "1") {
        structure_unionize <- .download_adult_mouse(gene, atlas_id, graph_id, experiment)
    } else if (atlas_id == "3") {
        structure_unionize <- .download_dev_mouse(gene, atlas_id, graph_id, experiment)
    } else {
        stop("An Error Has Occured with Altas ID's", call. = FALSE)
    }
    
    ########################################################################
    # get experiment number for a gene; can be multiple experiments per gene
    if (atlas_id == "1") {
        ontology <- .fetch_ontology(graph_id)
    } else if (atlas_id == "3") {
        ontology <- .fetch_ontology("1,17")
    }

    # ensure struct_depth makes sense
    cond <- struct_depth > min(ontology$depth) & struct_depth <= max(ontology$depth)
    if (!cond) {
        msg <- "Structure depth is too deep or too shallow; default of 3" 
        stop(msg, call. = FALSE)
    }

    
    # select important variables and merge to create a new frame
    ont_c <- c("id", "name", "acronym",
               "color_hex_triplet", "depth", "structure_id_path", 'graph_order')
    onto_str <- merge(
        ontology[ont_c], 
        structure_unionize, 
        by.x="id", by.y="structure_id"
    )
    
    # exlude objects that are not at our structure depth
    include  <- onto_str$depth == struct_depth
    onto_str <- onto_str[include,]
        
    if (atlas_id == "1") {
        .plot_mouse_substructures(onto_str, atlas, altas_id, gene, struct_depth,
                                  im_height, im_width)
    } else if (atlas_id == "3") {
        .plot_devmouse_substructures(onto_str, atlas, atlas_id, gene, struct_depth,
                                    im_height, im_width)
    } else {
        stop("An Error Has Occured with Atlas ID's", call. = FALSE)
    }
}
