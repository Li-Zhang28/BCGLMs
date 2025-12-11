#' Similarity matrix calculation
#' @description We provide phylogenetic similarity choice to capture the phylogenetic relatedness among taxa.
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param dist method to calculate dissimilarity matrix, default is Bray-curtis
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @author Li Zhang
#' @references \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9946}
similarity <- function(dis.taxa) {
  dis.taxa = as.matrix(dis.taxa); dim(dis.taxa)

  # --- CHECK: must be square and symmetric ---
  if (!is.matrix(dis.taxa) || nrow(dis.taxa) != ncol(dis.taxa)) {
    stop("dis.taxa must be a square matrix.")
  }
  if (!all.equal(dis.taxa, t(dis.taxa), check.attributes = FALSE)) {
    stop("dis.taxa must be symmetric.")
  }

  ###similarity matrix
  simi.mat <- function(dis.mat)
  {
    n <- ncol(dis.mat)
    D2 <- dis.mat^2
    I <- diag(1, n)
    II <- array(1, c(n,n))
    K <- -0.5 * (I - II/n) %*% D2 %*% (I - II/n)
    # to ensure a positive semi-definite matrix
    # method 1
    ev <- eigen(K)
    v <- ev$vectors
    e <- abs(ev$values)
    K <- v %*% diag(e) %*% t(v)
    K
  }

  K.taxa = simi.mat(dis.taxa); dim(K.taxa)

  taxa = as.numeric(K.taxa)
  node1 = node2 = w = NULL
  for (i in 1:(ncol(K.taxa)-1))
    for (j in (i+1):ncol(K.taxa))
    {
      if(abs(K.taxa[i,j]) > 0.18) ###0.18
      {
        node1 = c(node1,i);
        node2 = c(node2,j);
        w = c(w, K.taxa[i,j])
      }
    }
  w = sqrt(abs(w)); length(w)

  return(list(node1=node1,
            node2=node2,
            w=w))

}
