#'
#'
#'
mouse_subregions_plot <- function(gene, experiment=NULL, struct_depth = 3) {
    
    # get experiment number for a gene; can be multiple experiments per gene
    if (is.null(experiment)) {
        experiment <- .pick_experiment(gene)
    }

    # pull ontology and structure energy from API
    structure_unionize <- .fetch_structure_unionize(experiment)
    ontology <- .fetch_ontology()

    # ensure struct_depth makes sense
    stopifnot(struct_depth > min(ontology$depth) & struct_depth < max(ontology$depth))

    # # get a factor variable with the names
    # struct_bool <- ontology$depth == struct_depth
    # struct_names <- ontology[c("name", "structure_id_path")][struct_bool,]

    # select important variables and merge to create a new frame
    ont_c <- c("id", "name", "acronym",
               "color_hex_triplet", "depth", "structure_id_path") 
    str_c <- c("structure_id", "expression_energy")
    onto_str <- merge(
        ontology[ont_c], 
        structure_unionize[str_c], 
        by.x="id", by.y="structure_id"
        )

    # exlude objects above our granularity
    include  <- onto_str$depth == struct_depth
    onto_str <- onto_str[include,]
# 
#     # truncate the paths to merge with struct_names
#     onto_str$structure_id_path <- sapply(onto_str$structure_id_path,
#         function(x){
#             splt_path <- strsplit(x, '/')[[1]]
#             rtrn <- ""
#             for (i in 2:(struct_depth+2)) {
#                 rtrn <- paste(rtrn, splt_path[i], sep="/" )
#             }
#             paste(rtrn,"/", sep="")
#         })
# 
#     # merge our structure names
#     onto_str <- merge(onto_str, struct_names, by="structure_id_path")

    ylim <- c(0, quantile(onto_str$expression_energy, 0.995))

    # create pdf; tryCatch to ensure proper resource closing
    tryCatch({
        pdf(paste(gene, "adultMouse_brainRegionBarplot_ExpressionEnergies_depth", 
                struct_depth, "experiment", experiment, ".pdf", sep="_"),
                height=9,width=20)
        par(mfrow=c(2,1),mar=c(5,10,4,2))
        
        # plot the structures
        .verboseBarplot2(onto_str$expression_energy, onto_str$name, main=gene,
                         color=paste("#",onto_str$color_hex_triplet, sep=""),
                         las=2,xlab="",ylab="Mean Expression Energy",ylim=ylim)
        }, finally = {
            dev.off()
        })
}

#------------------------------HELPER FUNCTIONS--------------------------------#
.fetch_mouse_experiments <- function(gene) {
    # retrieves all of the experiments for a given gene in the mouse atlas
    # Arguments
    # ---------
    # gene -- a single gene acronym; spelling and capitalization important

    plane <- list("1" = "coronal", "2" = "sagital")
    
    # define query and url
    set <- "SectionDataSet"
    # eq right now, il for fuzzy searching or find table fo human mouse homology
    query <- paste("query.json?criteria=products[id$eq1],genes[acronym$eq'",
                gene,
                "'],section_images[failed$eqfalse]&include=genes,section_images",
                sep="")
    
    # pull URL and save locally; result is a bool if failed
    URL <- .construct_api_url(set, query)
    result <- .safe_api_call(URL)

    # see if the api call was succesful; if not stop
    if (is.logical(result)) {
        msg <- paste(gene, "is not a valid gene input")
        stop(msg, call.=FALSE)
    }
    
    # check if the gene returned a valid message; if not may not exist for mouse
    if (length(result$msg) == 0) {
        msg <- paste("API Call successful but no content delivered for", gene, ".\nCheck spelling and capitalization; there may be no experiments with", gene)
        stop(msg, call.=FALSE)    
    }
    
    # extract the experiment numbers and thier plane
    experiments <- result$msg[c("id", "plane_of_section_id")]
    experiments$plane <- sapply(experiments$plane_of_section_id,
                                          function(x){plane[[as.character(x)]]})

    # remove id and return; this format is easy to scan
    experiments[c("id", "plane")]
}


.fetch_structure_unionize <- function(exp_id) {
    # fetch the structure unionization for an experiment
    # Arguments
    # ---------
    # exp_id -- id for an experiment with a gene; from fetch_mouse_experiments
    set <- "StructureUnionize"
    query <- paste("query.json?criteria=[section_data_set_id$eq",exp_id,"],structure[graph_id$eq1]", sep="")
    URL <- .construct_api_url(set, query)
    
    result <- .safe_api_call(URL)
    
    # check that experiment ID is valid
    if (is.logical(result)) {
        msg <- paste(exp_id, "is not a valid experiment ID")
        stop(msg, call.=FALSE)
    }
    
    # ensure there is content in the message
    if (length(result$msg) == 0) {
        msg <- paste("API Call successful but no content delivered for", exp_id, ".\nCheck to ensure this is a valid experiment")
        stop(msg, call.=FALSE)    
    }
    result$msg
}


.fetch_ontology <- function() {
    set <- "Structure"
    query <- "query.json?criteria=[graph_id$eq1]&num_rows=2000"
    URL <- .construct_api_url(set, query)
    
    result <- .safe_api_call(URL)
    result$msg
}


.pick_experiment <- function(gene) {
    # function for user to pick experiment if one isn't supplied
    # Arguments
    # ---------
    # gene -- gene acronym used in allen institutes for a mouse
    experiments <- .fetch_mouse_experiments(gene)
    if (nrow(experiments) == 1) {  # if only one return it
        return(experiments$id)
    }
    
    print("There are multiple experiments for this gene.")
    print(experiments)
    
    # loop until user has valid selection
    while (TRUE) {
        experiment <- readline("Please pick one -> ")
    
        if (as.integer(experiment) %in% experiments$id) {
            return(experiment)
        } else {
            print("NOT A VALID SELECTION")
        }
    }
}