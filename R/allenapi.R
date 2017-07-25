.safe_api_call0 <- function(URL) {
     # Call API with less boiler plate and ensure error handling
     tryCatch( {
         URL <- paste(URL, "&num_rows=all", sep="")
         result <- jsonlite::fromJSON(URL)
     }, error = function(e) {
         msg <- paste("API URL not available or changed\nJSON error:", e)
         stop(msg, call.=FALSE )
     })
     
     if (!result$success) {
         return(FALSE)
     } 
     
     return(result)
 }


.json_cache <- function(f) {
    # decorator function that returns a URL cacheing version of an API tracker
    # uses saveRDS and readRDS to save and load cached items
    function(...) {
        # create cache directory
        if (!("./.json_cache" %in% list.dirs(recursive=FALSE))) {
            print("creating cache directory")
            dir.create("./.json_cache")
        }

        # create unique URL name for cacheing
#        split_URL <- strsplit(URL,"/")[[1]]
#        s_URL <- split_URL[[length(split_URL)]]
#        f_URL <- ""
#        for (c in serialize(s_URL, connection=NULL)) {
#            f_URL <- paste(f_URL, c, sep="")
#        }        
        f_hash <- digest::sha1(...)
        cache_f <- paste(f_hash,".cache.rda", sep="")
        cache_path <- paste("./.json_cache/", cache_f, sep="")
        
        if (cache_f %in% list.files("./.json_cache", recursive=FALSE)) {
            print("Fetching cache")
            result <- readRDS(file = cache_path)
            
        } else {
            print("No cache, requesting")
            result <- f(...)
            saveRDS(result, file = cache_path)
        }

        result
    }
}


# verson of safe_api using caching
.safe_api_call <- .json_cache(.safe_api_call0)


.construct_api_url <- function(set, query) {
    # pastes the url for query and ensures that the api base is the same for all
    # searches (only needs to be changed once if the site changes)
    api_url <- "http://api.brain-map.org/api/v2/data"
    URL <- paste(api_url,set,query, sep="/")
    URL
}

#------------------------------HELPER FUNCTIONS--------------------------------#
.fetch_ontology <- function(graph_id) {
    # smaller helper function to fetch ontology via graph ID
    set <- "Structure"
    query <- paste("query.json?criteria=[graph_id$eq", graph_id,"]", sep="")
    URL <- .construct_api_url(set, query)
    result <- .safe_api_call(URL)
    result$msg
}


#-----------------------------------LEGACY-------------------------------------#
## Helper functions for the allen institute API
#.safe_api_call0 <- function(URL, n_row = 2000) {
#    # Call API with less boiler plate and ensure error handling
#    all_results <- list()
#    rows_left <- TRUE
#    start_row <- 0
#    while (rows_left) {
#        URL_next <- paste(URL, "&start_row=", start_row,"&num_rows=", n_row, sep="")    
#        tryCatch( {
#            result <- jsonlite::fromJSON(URL_next)
#        }, error = function(e) {
#            msg <- paste("API URL not available or changed\nJSON error:", e)
#            stop(msg, call.=FALSE )
#        })
#        
#        if (!result$success) {
#            return(FALSE)
#        } 
#        
#        if (length(result$msg) == 0) {
#            rows_left <- FALSE
#        } 
#        
#        all_results$msg <- rbind(all_results$msg, result$msg)
#        start_row <- start_row + n_row
#    }
#    
#    return(all_results)
#}
#        illegal <- c("\\$", "\\?", "\\%", "\\&", "\\{", "\\}","\\=", "\\[", "\\]", 
#                    "\\!", "\\'", "\\:", "\\@", "\\.", "\\(", "\\)")
