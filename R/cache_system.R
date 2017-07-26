# cacheing system for atlasplot

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
        .cache_manager(f_hash, manager_file)
        
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

        result
    }
}


.cache_manager <- function(f_hash, manager_file) {
    # .cache_manager function; responsible for keeping track of accesses
    manager <- readRDS(manager_file)

    if (f_hash %in% names(manager)) {
        manager[f_hash] <- manager[[f_hash]] + 1
    } else {
        manager[f_hash] <- 1
    }
    
    saveRDS(manager, manager_file)
}

