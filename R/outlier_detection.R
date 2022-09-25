outliers_tau_cutoff <- function(cutoff, mscale_bp) {
  function(ktau) {
    thr <- qchisq(cutoff, df = ktau$p)
    robust_scale <- mscale(u = ktau$di,
     c = normal_consistency_constants(ktau$p),
                            b = mscale_bp
                            )
    
    ktau$outliers <- which(ktau$di ^ 2 > thr * robust_scale ^ 2)
    
    return(ktau)
  }
  
}