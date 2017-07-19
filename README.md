# atlasplot
Simple Reproducible Plots for the Allen Institutes Brain Atlases

This package aim to simplify the process of producing graphics for the expression data available in the Allen Institutes Brain Atlases. Additionally, it will be coupled with individual data packages, allowing anyone to easily access this vast data resource.

## Install
To install the `atlasplot` package it is recommended to use the 
[`devtools`](https://github.com/hadley/devtools) package. This is a well trusted and reputable package from Hadley Wickham intended to streamline package development and installation.

In the `R` command line perform the following commands
```
# Install devtools package
install.package("devtools")

# Install hbadata
devtools::install_github("edelsonc/atlasplot")
```

## Usage
Here is a list of the current accompanying data packages:
* [`hbadata`](https://github.com/edelsonc/hbadata)
* [`fetdata`](https://github.com/edelsonc/fetdata)
* [`brspdata`](https://github.com/edelsonc/brspdata)

`atlasplot` additionally has an allen institute API wrapper for:
* Adult Mouse Brain Atlas
* Developing Mouse Brain Atlas
