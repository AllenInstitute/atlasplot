# Maintaining `atlasplot`
The `atlasplot` package and all its associated data packages were created using the 
[`devtools`](https://github.com/hadley/devtools). As such, it is highly recommended that
all package edits and updates are done through `devtools`.

## Basic Structure
The package contains a few important components

* `R` directory -- contains all R functions and logic
* `data-raw` directory -- appropriate directory for helper scripts and dev information
* `man` directory -- help documentation directory; maintained via `devtools::document()`
* `atlasplot.Rproj`, `NAMESPACE` -- package text files maintained via `devtools::document()`
* `LICENSE` -- license file; not to be changed without legal consultation
* `README.md` -- project README file; first thing people see
* `function_list.md` -- useful `.md` file that lists functions by location

## Workflow
The basic workflow of `atlasplot` is as follows

1) Create new R script and place in `R`

2) Open `R` and `setwd` to the location of `atlasplot`
```
> setwd("/path_to_package/atlasplot")
```

3) Load all functions via `devtools::load_all()` and test new functions. Repeat this step often.
```
> # example of testing new cool_plots function
> devtools::load_all()
> cool_plots()  # A REALLY COOL PLOT APPEARS
> # after adding a new parameter `cool_factor`
> devtools::load_all()
> cool_plots(cool_factor = ">9000")

```

4) Once satisfied with your new addition, add [ROxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html) docstrings in front of user functions. Then document with devtools
```
> devtools::document()
```

5) Finally, run `create_function_list.sh` to update `function_list.md`. This file will be
located in `data-raw`. To run, open a `bash` shell and go to the directory above `atlasplot`. Then run
```
$ bash create_function_list.sh atlasplot
```

This is the basic process of updating `atlasplot`. Additionally, I highly suggest commiting
early and often. For a good resource on using `git` see [Pro Git](https://git-scm.com/book/en/v2).

## Imporant Notes

#### Avoid using `library` in any functions 
Reference them directly using `::` notation

Eg.
```
> jsonlite::fromJSON("my_json.json")
> atlasplot::mouse_subregions_plot("Shh", "DevMouse", save_pdf = FALSE)
```
#### Export user functions
Functions are only available to users if you specify `@export` in the `ROxygen2` docstring.
In general, user functions should be snake_case and package functions should begin with a dot.

Eg. 
```
# file fun_function.R

#' fun_function
#' 
#' @param name name of the person who will have fun!
#' @export
fun_function <- function(name) {
    msg <- .my_hidden_function(name, "do some science!")
    print(msg)
}

.my_hidden_function <- function(name, activity) {
    paste("Hi", name, "lets", activity)
}

```
`fun_function` will be exported and user accessible in the final package while
`.my_hidden_function` will only be internal accessible.

#### `sysdata.rda` for small internal objects
`sysdata.rda` is located in the `R` directory and hold important local pieces of data.

*Be careful not to overwrite it when adding new internal data!*

# Additional Resources
I highly recommend Hadley Wickham's book [R Packages](http://r-pkgs.had.co.nz/). It is comprehensive and very readible.