#' Bayesian Compositional GLMs
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param y outcome (binary or continuous)
#' @param family gaussian or bernoulli
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @examples dat=sim_c(200,100)
#'           fit=bcglm(dat$x,dat$y)
#'           summary(fit)
#'           fixef(fit)
#'           mcmc_plot(fit)
#'
#' @author Li Zhang
#' @references \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9946}
bcglm <- function(x,y,family=gaussian,similarity=NULL) {

  log_x=log(x)

  fm = bf(y ~ ., center=T)

  dat=data.frame(y,log_x)
  ### In real data analysis, the number of identified taxa in univariate generalized linear models can serve as a reasonable
  ### prior guess for 𝑚0
 if(is.null(similarity)){
   bp4 = set_prior("horseshoe(df=3, df_global=3)", class="b")
   bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)

   f4= brm(fm, data=dat,family=family,prior=bp4, control = list(adapt_delta = 0.99,max_treedepth= 18),
           chains=2, iter=2000)

 }

  if(!is.null(similarity)){
    node1=similarity$node1
    node2=similarity$node2
    w=similarity$w

    bp4 = set_prior("horseshoe(df=3, df_global=3)", class="b")
    bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)
    bp4 = bp4 + set_prior("target += -0.5*dot_self(w .* (log(hs_local[node1])-log(hs_local[node2])))", check=F)

    ln = length(node1)
    stanvars =stanvar(x=ln, name="ln", scode="int ln;", block="data") +
      stanvar(x=node1, name="node1", scode="int node1[ln];", block="data") +
      stanvar(x=node2, name="node2", scode="int node2[ln];", block="data") +
      stanvar(x=w, name="w", scode="vector[ln] w;", block="data")

    f4= brm(fm, data=dat,family=family,prior=bp4, stanvars=stanvars,control = list(adapt_delta = 0.99,max_treedepth= 18),
            chains=2, iter=2000)

  }
return(f4)
}



