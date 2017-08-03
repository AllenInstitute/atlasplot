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

#--------------------------------USER OPTIONS-------------------------------------------#
#' set_cache_size
#' 
#' set_cache_size is a atlasplot helper function that allows the user to decide the amount
#' of storage they'd like to allocate to the on disk file cache. All of the functions
#' using the Allen Institute API make use of the cache. If cache size is exceeded, the
#' oldest item are removed.
#' 
#' @export
#' @param size number of bytes allocated to the cache; default 2e9. size = FALSE disables cacheing
#' @examples
#' # set cache to 10 kilobyte
#' set_cache_size(10000)
#'
#' # turn off cache
#' set_cache_size(FALSE)
#'
set_cache_size <- function(size) {

    if (missing(size)) {
        stop("Please specify a size in bytes", call. = FALSE)
    }

    if (is.numeric(size)) {
        size <- as.integer(size)
    }

    CACHE_HOME <- path.expand("~/.json_cache")
    manager_file <- paste(CACHE_HOME, "manager.rda", sep="/")

    # create home directoy and manager if it doesn't exist; otherwise edit manager
    if (!( CACHE_HOME %in% list.dirs("~", recursive=FALSE))) {
        print("creating cache directory")
        dir.create(CACHE_HOME)
        saveRDS(list("cache_size" = size), manager_file)
    } else {
        manager <- readRDS(manager_file)
        manager$cache_size <- size
        saveRDS(manager, manager_file)
    }
}


#' clear_cache
#'
#' Deletes all current files in the cache. This is implemented in the case that an item
#' on the API updates. Clearing the cache will make it so the new item is downloaded at
#' next run. Not necessary for short term use.
#'
#' @export
#' @examples
#' # delete cache items since a new ontology has been released!
#' clear_cache()
#'
clear_cache <- function() {
    # get system cache location and manager location
    CACHE_HOME <- path.expand("~/.json_cache")
    manager_file <- paste(CACHE_HOME, "manager.rda", sep="/")

    # ensure that a cache exists to clear
    if (!( CACHE_HOME %in% list.dirs("~", recursive=FALSE))) {
        msg <- paste("No cache in ", CACHE_HOME, sep="")
        stop(msg, call. = FALSE)
    }

    # read manager file and delete all files (manager included)
    manager <- readRDS(manager_file)
    size <- manager[["cache_size"]]
    files <- names(manager)[names(manager) != "cache_size"]
    for (file in files) {
        rm_file <- paste(CACHE_HOME, "/", file, ".cache.rda", sep="")
        .remove_file(rm_file)
    }

    saveRDS(list("cache_size" = as.integer(size)), manager_file)   
}


#------------------------------CACHE MANAGEMENT-----------------------------------------#
.json_cache <- function(f) {
    # decorator function that returns a URL cacheing version of an API tracker
    # uses saveRDS and readRDS to save and load cached items
    function(...) {
        # create cache directory and management file
        CACHE_HOME <- path.expand("~/.json_cache")
        manager_file <- paste(CACHE_HOME, "manager.rda", sep="/")
        if (!( CACHE_HOME %in% list.dirs("~", recursive=FALSE))) {
            print("creating cache directory")
            dir.create(CACHE_HOME)
            saveRDS(list("cache_size" = as.integer(2e9)), manager_file)
        }

        # hash input parameters to create unique file names
        f_hash <- digest::sha1(...)

        # manage cache updating
        cache_status <- .cache_manager(f_hash, CACHE_HOME, manager_file)


        if (cache_status) {
            # create write to path
            cache_f <- paste(f_hash,".cache.rda", sep="")
            cache_path <- paste("~/.json_cache/", cache_f, sep="")

            # check if the file is in the cache
            if (cache_f %in% list.files(CACHE_HOME, recursive=FALSE)) {
                print("Fetching cache")
                result <- readRDS(file = cache_path)

            } else {
                print("No cache, requesting")
                result <- f(...)
                saveRDS(result, file = cache_path)
            }
        } else {
            result <- f(...)
        }
        result
    }
}


.cache_manager <- function(f_hash, CACHE_HOME, manager_file) {
    # .cache_manager function; responsible for keeping track of accesses
    manager <- readRDS(manager_file) 

    if (manager[["cache_size"]] == FALSE) {
        return(FALSE)
    }

    max_size <- manager[["cache_size"]]
    current_size <- .cache_size(manager, CACHE_HOME)
    if (current_size > max_size) {
        manager <- .flush_cache(manager, CACHE_HOME, max_size)
    }
    
    t <- as.integer(as.POSIXct(Sys.time()))
    if (f_hash %in% names(manager)) {
        manager[f_hash] <- t
    } else {
        manager[f_hash] <- t
    }

    
    saveRDS(manager, manager_file)
    return(TRUE)
}


.cache_size <- function(manager, CACHE_HOME) {
    # helper function computes size of current cache
    files <- names(manager)[names(manager) != "cache_size"]
    if (length(files) == 0) {
        return(0)
    }

    files <- paste(CACHE_HOME,"/", files, ".cache.rda", sep="")
    f_size <- file.size(files)
    sum(f_size, na.rm=TRUE)
}


.flush_cache<- function(manager, CACHE_HOME, max_size) {
    # maintains and cleans the cache
    files <- unlist(manager)
    files <- files[names(files) != "cache_size"]
    size <- .cache_size(files, CACHE_HOME)
    files <- sort(files)

    # remove least used files until correct size or no files
    while (size > max_size & length(files) != 0) {
        rm_file <- paste(CACHE_HOME, "/", names(files)[[1]], ".cache.rda", sep="")

        # remove the least used file
        .remove_file(rm_file)

        # recompute file sizes
        files <- files[-1]
        size <- .cache_size(files, CACHE_HOME)
    }

    # recreate manager without the removed files
    manager <- list("cache_size" = max_size)
    for (i in 1:length(files)) {
        key <- names(files)[i]
        manager[key] <- files[i]
    }
    manager[ !sapply(manager, is.na) ]
}


.remove_file <- function(rm_file) {
    # small function to add some error prevention boiler plate
    if (!file.exists(rm_file)) {
        msg <- paste("Cache file", rm_file, "does not exist")
        stop(msg, call. = FALSE)
    }

    rm_bool <- file.remove(rm_file)
    if (!rm_bool) {
        msg <- "Cannot remove files"
        stop(msg, call. = FALSE)
    }
}


# cacheing version of .safe_api_call
.safe_api_call <- .json_cache(.safe_api_call0)
