# found at: http://augix.com/wiki/Make%20trees%20in%20R,%20test%20its%20stability%20by%20bootstrapping.html
# modified slightly to return tree and bp in list, and default to NJ method.

boot.tree <- function(data, B = 100, tree = "nj") {
    library(phangorn)
    if (tree == "upgma") {
        func <- function(x) upgma(dist(x, method = "euclidean"), method="average")
    }
    if (tree == "nj") {
        func <- function(x) nj(dist(x, method = "euclidean"))
    }
    if (tree == "hclust") {
        func <- function(x) {
            tr = hclust(dist(x, method = "euclidean"))
            tr = as.phylo(tr)
            return(tr)
        }
    }
    tr_real = func(data)
    plot(tr_real)
    library(ape)
    bp <- boot.phylo(tr_real, data, FUN=func, B=B)
    nodelabels(bp)
    
    ret = list();
    ret$bp = bp;
    ret$tree = tr_real;

    return(ret)
}
 
#data = t(USArrests) # columns are the branches
#boot.tree(data, B=1000, tree='hclust')


# get the newick format for the tree
# write.tree(ret$tree, file="nj.tree")
 

 
# Description:
#   This function builds a upgma or nj tree and tests its stability by bootstrapping. It returns a tree, and bootstrap result.
#
# Usage:
#   boot.tree(data, B = 100, tree = "upgma")
#
# Arguments:
#   data: a numeric matrix, data fram
#   
#   B: the number of bootstrap replicates. (100 by default).
#
#   tree: tree type. It could be "upgma", or "nj". ("upgma" by default.)
#
# Values:
#   It return a numeric vector which _i_th element is the times that we observe the _i_th node of the real tree.
#
# Examples:
#   # compare it with the hclust
#   par(ask = TRUE)
#   plot(hclust(dist(t(USJudgeRatings))))
#   boot.tree(t(USJudgeRatings), 100) 
