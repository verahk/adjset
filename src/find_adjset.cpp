#include <Rcpp.h>
using namespace Rcpp;


//' Identify adjustment set
//'
//' Identify different sets of nodes in a DAG valid for adjustment with respect to 
//' a cause node X and a effect node Y. Based on Zander and Liśkiewicz (2020) algorithm 
//' for finding (minimal) separators and a modified Algorithm 3.1 in Koller (2009, p. 75) 
//' for finding reachable nodes.
//' 
//' @name adjset
//' @param G a n-by-n adjacency matrix of a graph. `G[x, y] = 1` if there is an edge from node `x` to `y`, zero otherwise.
//' @param x,y (integer)
//'   column position of the cause and effect variable, respectively (zero-based index, from 0 to n-1).
//' @param Z (integer vector)
//'    column positions of conditioning variables (zero-based index, from 0 to n-1).
//' @param A (logical vector)
//'   indicates a subset of nodes paths must be contained in.
//'   The search along a path ends if visiting a node `i` such that `A[i] = FALSE`.
//' @param name (string) name of adjustment set, see details.
//' @param dmat (optional) a n-by-n adjacency matrix such that `dmat[x, y] = 1` if `x` is an ancestor of `y` in `G`.
//'   Added as an optional argument for speed. If not specified, `areDescendants` and `areAncestors` are called to 
//'   identify the relevant relationships.
//' @param Gx the backdoor graph of `G` w.r.t. `x`.
//' @details
//' The function `areReachable` implements Algorithm 3.1 in Koller (2009, p. 75) and
//' identifies nodes in `G` reachable from `x` given `Z` via an active trail in `G`, with the restriction
//' that all nodes on the trail are in `A`.
//' 
//' The function `find_nearest_adjset` is used to identify the `o-set` and minimal separators. 
//' It assumes `Z` is a valid adjustment set (a set of nodes that includes no descendants of `x` and separates `x` and `y` in the backdoor graph),
//' and returns the subset of `Z` reachable from `x` given `Z`. That is the subset of `Z` closest to `x`. 
//' A minimal separator is found by applying `find_nearest_adjset` again w.r.t `y` and the reduced set.
//' 
//' The function `find_adjset` applies these algorithms to identify different adjustment set: 
//' 
//' - "o": The optimal adjustment set (o-set). Given `A = Anc(x, y)`, the o-set is the set of nodes in `Z = A\Desc(x)` that can be reached from `y` via active paths in `A`.
//' - "o_min": The minimal O-set. The subset of the o-set reachable from `x`, i.e. a minimal separator of `x` and `y` in the backdoor graph.
//' - "pa_min": The minimal parent set. The subset of the parents of `x` reachable from `y`, i.e. a minimal separator of `x` and `y` in the backdoor graph.
//'   
//' These adjustment sets are defined only when `x` is ancestor of `y`.
//' If that is not the case, the function returns `y` indicating that $P(y|do(x)) = P(y)$, such that `x` has no causal effect on `y`.
//' Additionally, 
//' 
//' - "pa": the parents of `x`, read of `G` directly.
//' - "pa_if": if `y` is a descendant of `x`, then the set of parents is returned. Otherwise `y` is returned.


// [[Rcpp::export]]
bool isParent(NumericMatrix G, int x, int y) {
   return G(x, y) == 1;
 }
// [[Rcpp::export]]
bool isChild(NumericMatrix G, int x, int y) {
   return G(y, x) == 1;
}

// Recursive helper function
LogicalVector areFam_recursive(NumericMatrix G, int x,  bool (*isFun)(NumericMatrix, int, int), int n, LogicalVector out) {
   out[x] = 1;
   for (int i = 0; i < n; i++) {
      if (!out[i] && isFun(G, i, x)) {
         out[i] = 1;
         out = areFam_recursive(G, i, isFun, n, out);
      }
   }
   return out;
}

LogicalVector areFam(NumericMatrix G, IntegerVector nodes,  bool (*isFun)(NumericMatrix, int, int)) {
   int n = G.ncol();
   LogicalVector out(n);
   for (int j = 0; j < nodes.length(); j++) {
      int x = nodes[j];
      if (!out[x]) {
         out = areFam_recursive(G, x, isFun, n, out);
      }
   }
   return out;
}

// [[Rcpp::export]]
LogicalVector areDescendants(NumericMatrix G, IntegerVector nodes) {
   return areFam(G, nodes, isChild);
}

// [[Rcpp::export]]
LogicalVector areAncestors(NumericMatrix G, IntegerVector nodes) {
   return areFam(G, nodes, isParent);
}


// [[Rcpp::export]]
bool isDescendant(NumericMatrix G, int x, int y) {
   IntegerVector yvec(1, y);
   LogicalVector descendants = areDescendants(G, yvec);
   return descendants[x] == 1;
}



//' @rdname adjset
//' @export
//' @returns
//' - `areReachable`: a logical vector indicating which nodes in `G` can be reached from `x`
//'    through active paths that traverses only nodes in `A` conditional on `Z`.
// [[Rcpp::export]]
LogicalVector areReachable(NumericMatrix G,
                          int x,
                          IntegerVector Z,
                          LogicalVector A) {
   
   int n = G.ncol();
   
   // Find all ancestors of Z 
   LogicalVector AncZ = areAncestors(G, Z); 
   
   // Starting from x, traverse trails  ----
   IntegerVector Ly(1, x);       // nodes (start from x)..
   IntegerVector Ld(1, 0);       // .. and directions to be visited (start from below)
   LogicalMatrix V(2, n); // visisted nodes and dir
   LogicalVector R(n);    // reachable nodes
   
   while (Ly.length() > 0){
      Rcpp::checkUserInterrupt();
      
      int y = Ly[0];
      int d = Ld[0];
      
      Ly.erase(0);
      Ld.erase(0);
   
      if (A[y]) {       // only care about paths traversing nodes in A
         if (!V(d, y)){ // if not already visited
            // Rcout << "The value of y : " << y << "\n";
            // Rcout << "The value of d : " << d << "\n";
            
            R[y] = true;    // mark y as reachable
            V(d, y) = true; // track visited nodes
            
              
            // check if y is in conditioning set Z
            bool inZ = false;
            for (int i = 0; i<Z.length(); i++) {
             if (Z[i] == y) {
                inZ = true;
                break;
             }
            }
            
            // add nodes along active paths to be visited next
            if (d == 0 && !inZ) {
               // if arrived y from the bottom (from a child) and Y is not in Z, 
               // active paths continue both via y's children and parents
               for (int i = 0; i<n; i++) {
                  if (G(y, i) == 1 && !V(1, i)) {
                      // continue via children of y from the top
                      Ly.push_back(i);
                      Ld.push_back(1);
                  } else if (G(i, y) == 1 && !V(0, i)) {
                      // continue via parents of y from the bottom
                      Ly.push_back(i);
                      Ld.push_back(0);
                } 
             }
            } else if (d == 1) {
               // if arrived y from the top (from a parent), 
               // paths via y's children are active if y is not observed.
               // Paths via y's parents are blocked, unless a descendant of y is an ancestor of Z.
               for (int i = 0; i<n; i++) {
                  if (!inZ && G(y, i) == 1 && !V(1, i)) {
                      // continue via children of y from the top
                      Ly.push_back(i);
                      Ld.push_back(1);
                  } else if (AncZ[y] && G(i, y) == 1 && !V(0, i)) {
                     // continue via parents of y from the bottom
                     Ly.push_back(i);
                     Ld.push_back(0);
                  } 
               }
            }
         }
      }
   }
   return R;
}

NumericMatrix backdoorgraph(NumericMatrix G, int x) {
   int n = G.ncol();
   NumericMatrix bG = clone(G);
   for (int i = 0; i < n; i++) {
      if (G(x, i) == 1) {
         bG(x, i) = 0;
      }
   }
   return bG;
}


//' @rdname adjset
//' @returns
//' - `find_nearest_adjset`: The subset of Z nearest `x` valid for adjustment (assuming Z is a valid adjusment set).
//' @export
// [[Rcpp::export]]
IntegerVector find_nearest_adjset(NumericMatrix Gx, int x, IntegerVector Z, LogicalVector A) {

   if (Z.length() == 0 || na_omit(Z).length() == 0) {
      return Z;
   } else {
      LogicalVector R = areReachable(Gx, x, Z, A);
      LogicalVector indx = R[Z];
      IntegerVector out = Z[indx];
      return out;
   }
}

//' @rdname adjset
//' @export
//' @return
//' - `adjset`: an integer vector with the (zero-indexed) column positions of an adjustment set as given by `name`.
// [[Rcpp::export]]
IntegerVector find_adjset(NumericMatrix G, 
                           int x, 
                           int y, 
                           String name, 
                           NumericMatrix dmat = NumericMatrix(0),
                           NumericMatrix Gx = NumericMatrix(0)) {

   int n = G.ncol();
   IntegerVector Z = IntegerVector::create(); // init vector for adjustment set 
   LogicalVector De(n);                       // init indicator for descendants
   
   if (name != "pa") {
      if (dmat.ncol() == n) {
         De =  dmat(x, _) > 0;
      } else {
         De = areDescendants(G, IntegerVector(1, x));
      }
      if (De(y) == 0) {
        // Rcout << "Y is not descendant of X \n";
        // return `y` if `y` is not a descendant of `x`
        return IntegerVector(1, y);
      }
   }
   

   if (name == "pa" || name == "pa_if" || name == "pa_min") {
      // list all parents of `x`
      for (int i = 0; i < n; i++) {
         if (G(i, x) == 1) {
            Z.push_back(i);
         }
      }
      if (name != "pa_min" || Z.length() == 0) {
         return Z; 
      }
   }
   
   // compute backdoor graph
   if (Gx.ncol() != n) {
      Gx = backdoorgraph(G, x);
   }
   
   // find common ancestors of x and y, potential confounders
   LogicalVector A(n);     
   if (dmat.ncol() == n) {
      A  = dmat(_, x) + dmat(_, y) > 0;
   } else {
      IntegerVector nodes = {x, y};
      A = areAncestors(G, nodes);
   }
 
   if (name == "pa_min") {
      return find_nearest_adjset(Gx, y, Z, A);
   } else {
      
      // list set of common ancestors minus descendants - potential confounders
      for (int i = 0; i < n; i++) {
         if (A(i) && !De(i)) {
            Z.push_back(i);
         }
      }
      
      if (Z.length() == 0) {
         return Z;
      } else if (name == "o") {
         return find_nearest_adjset(Gx, y, Z, A);
      } else if (name == "o_min") {
         Z = find_nearest_adjset(Gx, y, Z, A);     // o-set
         return find_nearest_adjset(Gx, x, Z, A);  // minimal o-set
      } else {
         stop("Invalid value of argument `name`");
      }
   }
}


/*** R

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
x <- 4-1
y <- 6-1

bdag <- dag
bdag[x+1, ] <- 0
areDescendants(dag, x)
A <- areAncestors(dag, x) | areAncestors(dag, y)
Z0 <- which(A & !areDescendants(dag, x))-1
areAncestors(dag, Z0)
areReachable(bdag, y, Z0, A)
colnames(dag)[Z0[areReachable(bdag, y, Z0, A)[Z0+1]]+1]
find_nearest_adjset(dag, y, Z0, A)

find_nearest_adjset(dag, x, find_nearest_adjset(dag, x, Z0, A), A)

find_adjset(dag, y, x, "o_min")


*/
