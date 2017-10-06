# atlasplot
Simple Reproducible Plots for the Allen Institute's Brain Atlases

This package aims to simplify the process of producing graphics for the expression data available in the Allen Institute's Brain Atlases. Additionally, it will be coupled with individual data packages, allowing anyone to easily access this vast data resource.

## Install
To install the `atlasplot` package it is recommended to use the 
[`devtools`](https://github.com/hadley/devtools) package. This is a well trusted and reputable package from Hadley Wickham intended to streamline package development and installation.

In the `R` command line perform the following commands
```
# Install devtools package
install.package("devtools")

# Install atlasplot
# build_vignettes = TRUE creates the html vignette
devtools::install_github("AllenInstitute/atlasplot", build_vignettes = TRUE)
```

## Usage
Here is a list of the current accompanying data packages:
* [`hbadata`](https://github.com/AllenInstitute/hbadata)
* [`fetdata`](https://github.com/AllenInstitute/fetdata)
* [`brspdata`](https://github.com/AllenInstitute/brspdata)
* [`nhpdata`](https://github.com/AllenInstitute/nhpdata) 

`atlasplot` additionally has an Allen Institute API wrapper for:
* Adult Mouse Brain Atlas
* Developing Mouse Brain Atlas

#### Caching
`atlasplot` makes use of a simple caching system to speed up access to commonly used resources.
By caching common API calls locally, `atlasplot` can significantly reduce network communication time.
Caching is initialized the first time you call a function that uses to web API. The cache
file will always be located in `~/.json_cache`. Files are saved as `.rda` objects, and
are compressed to decrease their footprint.

There are two user functions associated with caching:
* `set_cache_size`
* `clear_cache`

By default, `atlasplot` limits its cache to 2 Gb. `set_cache_size` can be used to change
this default

```
> set_cache_size(10000) # change cache size to 10 kb
```

Additionally, caching can be turned off by selecting a size of `FALSE`.

`clear_cache` deletes all current files in the users cache. It is called with no arguments
and is only inteneded to be used when updated data has been added to the API.
