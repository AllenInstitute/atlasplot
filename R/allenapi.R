# Helper functions for the allen institute API
.safe_api_call <- function(URL, n_row = "1000000") {
    # Call API with less boiler plate and ensure error handling
    tryCatch( {
        URL <- paste(URL, "&num_rows=", n_row, sep="")
        result <- jsonlite::fromJSON(URL)
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
    
    
    return(result)
}


.construct_api_url <- function(set, query) {
    # pastes the url for query and ensures that the api base is the same for all
    # searches (only needs to be changed once if the site changes)
    api_url <- "http://api.brain-map.org/api/v2/data"
    URL <- paste(api_url,set,query, sep="/")
    URL
}
