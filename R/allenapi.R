# Helper functions for the allen institute API
.safe_api_call <- function(URL, n_row = 2000) {
    # Call API with less boiler plate and ensure error handling
    all_results <- list()
    rows_left <- TRUE
    start_row <- 0
    while (rows_left) {
        URL_next <- paste(URL, "&start_row=", start_row,"&num_rows=", n_row, sep="")    
        tryCatch( {
            result <- jsonlite::fromJSON(URL_next)
        }, error = function(e) {
            msg <- paste("API URL not available or changed\nJSON error:", e)
            stop(msg, call.=FALSE )
        })
        
        if (!result$success) {
            return(FALSE)
        } 
        
        if (length(result$msg) == 0) {
            rows_left <- FALSE
        } 
        
        all_results$msg <- rbind(all_results$msg, result$msg)
        start_row <- start_row + n_row
    }
    
    return(all_results)
}


.construct_api_url <- function(set, query) {
    # pastes the url for query and ensures that the api base is the same for all
    # searches (only needs to be changed once if the site changes)
    api_url <- "http://api.brain-map.org/api/v2/data"
    URL <- paste(api_url,set,query, sep="/")
    URL
}



#-----------------------------------LEGACY-------------------------------------#

# .safe_api_call <- function(URL, n_row = "4000") {
#     # Call API with less boiler plate and ensure error handling
#     tryCatch( {
#         URL <- paste(URL, "&num_rows=", n_row, sep="")
#         result <- jsonlite::fromJSON(URL)
#     }, error = function(e) {
#         msg <- paste("API URL not available or changed\nJSON error:", e)
#         stop(msg, call.=FALSE )
#     })
#     
#     if (!result$success) {
#         return(FALSE)
#     } 
#     
#     return(result)
# }
