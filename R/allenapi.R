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

#-------------------------------API FUNCTIONS-------------------------------------------#
.safe_api_call0 <- function(URL) {
     # Call API with less boiler plate and ensure error handling
     tryCatch( {
         URL <- paste(URL, "&num_rows=all", sep = "")
         result <- jsonlite::fromJSON(URL)
     }, error = function(e) {
         msg <- paste("API URL not available or changed\nJSON error:", e)
         stop(msg, call. = FALSE )
     })

     if (!result$success) {
         return(FALSE)
     }

     return(result)
 }


.construct_api_url <- function(set, query) {
    # pastes the url for query and ensures that the api base is the same for all
    # searches (only needs to be changed once if the site changes)
    api_url <- "http://api.brain-map.org/api/v2/data"
    URL <- paste(api_url, set, query, sep = "/")
    URL
}

#------------------------------HELPER FUNCTIONS--------------------------------#
.fetch_ontology <- function(graph_id) {
    # smaller helper function to fetch ontology via graph ID
    set <- "Structure"
    query <- paste("query.json?criteria=[graph_id$in", graph_id, "]", sep = "")
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
