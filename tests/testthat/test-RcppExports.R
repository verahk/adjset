

test_that("find_adjset returns the same sets as pcalg", {
  
  ## Check that adjustment sets returned by bida:::adjset_from_dag
  ## coincide with those from pcalg::optAdjSet (o-set) and pcalg::adjustment (minimal sets)
  
  adjsets <- c("o", "o-min", "pa-min", "pa", "anc")
  n <- 10
  ngraphs <- 10
  verbose <- T
  
  for (g in 1:ngraphs){
    set.seed(g)
    dag <- as(pcalg::randDAG(n, 3, weighted = FALSE), "matrix")
    colnames(dag) <- rownames(dag) <- paste0("X", seq_len(n))
    dmat  <- sign(round(solve(diag(n)-dag), 1))
    tdag <- t(dag)
    
    #pcalg::plot(as(dag, "graphNEL"))
    for (x in 1:n){
      for (y in seq_len(n)[-x]){
        if (verbose) cat("\n dag:", g, "x:", x, "y:", y)
        adjsets <- c("anc", "o", "o-min", "pa-if", "pa-min")
        sets <- lapply(adjsets, function(a) adjset(dag, x-1, y-1, a)+1)
        names(sets) <- adjsets
        
        if (dmat[x, y] == 0) {
          expect_equal(unique(sets), list(y))
        } else {
          expect_equal(sets$anc, 
                       unname(which(replace(dmat[, x] | dmat[, y] & !dmat[x, ],  x, FALSE))))
          expect_equal(sets$o, sort(pcalg::optAdjSet(tdag, x, y)))
          expect_equal(sets$`pa-if`, unname(which(dag[, x] == 1)))
          
          # compute all minimum sets
          min <- pcalg::adjustment(tdag, "dag", x, y, "minimal") 
          # loop through all minimal sets and compare with minimal o/pa
          o_equal <- pa_equal <- FALSE
          for (z in min) {
            z <- sort(z)
            if (!o_equal)  o_equal <- length(sets$`o-min`) == length(z) && all(sets$`o-min`== z)
            if (!pa_equal) pa_equal <- length(sets$`pa-min`) == length(z) && all(sets$`pa-min` == z)
            if (o_equal && pa_equal) break
          }
          
          expect_true(o_equal)
          expect_true(pa_equal)
          
          if (FALSE) {  # debugging
            # plot each adjustment set
            color <- c("grey", "lightblue",  "blue", "pink", "red")
            names(color) <- adjsets
            par(mfrow = c(1, 5),
                mar = c(5, 1, .1, .1))
            
            for (a in names(sets)) {
              z <- sets[[a]]
              tmp <- rep(color[a], length(z))
              names(tmp) <- colnames(dag)[z]
              Rgraphviz::plot(graph::graphAM(dag, edgemode = "directed"),
                              nodeAttrs = list(fillcolor = tmp),
                              main = a)
            }
            
            bdag <- dag
            bdag[x, ] <- 0
            Z0 <- sets$anc
            
            cat("\nZ0")
            Z0-1
            cat("\nAncestors of Z0")
            which(areAncestors(bdag, Z0-1))-1
            cat("\nNodes in Z0 reachable from y")
            Z0[areReachable(bdag, y-1, Z0-1, dmat[, x] | dmat[, y])[Z0]]
            find_nearest_adjset(bdag, y-1, sets$anc-1, dmat[, x] | dmat[, y])+1
            adjset(dag, x-1, y-1, "o")
          }
        }
      }
    }
  }
})
