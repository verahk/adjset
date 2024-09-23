
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adjset

<!-- badges: start -->
<!-- badges: end -->

This package contains code for identifying valid adjustment sets in a
DAG, including the o-set, the minimal o-set and the minimal parent set.

## Installation

You can install the development version of adjset from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("verahk/adjset")
```

## Example

``` r
library(adjset)

# specify DAG row-wise:
dag <- rbind(Z1  = c(0, 0, 0, 1, 0, 0, 0),
             Z2  = c(0, 0, 0, 1, 0, 0, 0),
             L   = c(0, 1, 0, 0, 1, 0, 0),
             X   = c(0, 0, 0, 0, 1, 0, 0),
             M   = c(0, 0, 0, 0, 0, 1, 0),
             Y   = c(0, 0, 0, 0, 0, 0, 0),
             U   = c(0, 0, 0, 0, 0, 1, 0))
colnames(dag) <- rownames(dag)


# compute adjustment set w.r.t. X and Y
x <- 4
y <- 6
adjsets <- c("anc", "pa","pa-min", "o", "o-min")
sets <- lapply(adjsets,
               function(a) adjset(dag, x-1, y-1, a)+1)
names(sets) <- adjsets
sets
#> $anc
#> [1] 1 2 3 7
#> 
#> $pa
#> [1] 1 2
#> 
#> $`pa-min`
#> [1] 2
#> 
#> $o
#> [1] 3 7
#> 
#> $`o-min`
#> [1] 3

# compare adjustment sets 
g <- graph::graphAM(dag, edgemode = "directed")

color <- c("grey", "lightblue",  "blue", "pink", "red")
names(color) <- adjsets
par(mfrow = c(1, 5),
    mar = c(5, 1, .1, .1))

for (a in names(sets)) {
  z <- sets[[a]]
  tmp <- rep(color[a], length(z))
  names(tmp) <- colnames(dag)[z]
  Rgraphviz::plot(g,
                  nodeAttrs = list(fillcolor = tmp),
                  main = a)
}
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r



# compare runtimes with alternative routines
microbenchmark::microbenchmark(pcalg::optAdjSet(t(dag), x, y),
                               adjset(dag, x-1, y-1, "o")+1,
                               check = "equivalent")
#> Unit: microseconds
#>                                expr      min        lq       mean    median
#>      pcalg::optAdjSet(t(dag), x, y) 2691.269 2737.7595 2911.78811 2787.0810
#>  adjset(dag, x - 1, y - 1, "o") + 1   13.207   15.0745   23.11478   23.1685
#>         uq      max neval cld
#>  2857.1695 5357.698   100   b
#>    27.5525   67.641   100  a

microbenchmark::microbenchmark(which(dag[, x] == 1),
                               adjset(dag, x-1, y-1, "pa")+1,
                               check = "equivalent")
#> Unit: microseconds
#>                                 expr   min     lq    mean median    uq    max
#>                 which(dag[, x] == 1) 1.730 1.7985 1.94152  1.846 1.886  8.221
#>  adjset(dag, x - 1, y - 1, "pa") + 1 3.036 3.1340 3.57394  3.192 3.282 36.834
#>  neval cld
#>    100  a 
#>    100   b
```
