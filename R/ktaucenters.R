#' ktaucenters
#'
#' A robust and efficient version of kmeans algorithm
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param centers either the number of clusters, say *k*, or a set of initial (distinct) cluster centers.
#' @param max_iter a maximum number of iterations used for the algorithm stopping rule
#' @param tolerance a tolerance parameter used for the algorithm stopping rule
#' @param n_runs the number of trials that the base algorithm must run. If it is greater than 1 and center is not set as NULL, a random set of (distinct) rows in \code{x} will be chosen as the initial centers.
#' @param flag_outliers optional argument for outliers detection - quantiles
#'     of chi-square to be used as a threshold for outliers detection,
#'     defaults to 0.999

#' @return A list including the estimated K centers and labels for the observations
##' \itemize{
##'  \item{\code{centers}}{:   matrix of size K x p, with the estimated K centers.}
##'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
##'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
##'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights
##' associated to each observation.}
##'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some
##' iteration there were clusters totally empty}
##'  \item{\code{niter}}{: number of iterations until convergence is achieved
##' or maximun number of iteration is reached}
##'  \item{\code{di}}{: distance of each observation to its assigned cluster-center}
##'  \item{\code{outliers}}{: indices observation that can be considered as outliers}
##' }

#' @examples
#' # Generate Sinthetic data (three cluster well separated)
#' Z <- rnorm(600);
#' mus <- rep(c(-3, 0, 3), 200)
#' X <-  matrix(Z + mus, ncol=2)
#'
#' # Generate 60 sinthetic outliers (contamination level 20%)
#' X[sample(1:300,60), ] <- matrix(runif(40, 3 * min(X), 3 * max(X)),
#'                                 ncol = 2, nrow = 60)
#'
#' ### Applying the algorithm ####
#' sal <- ktaucenters(
#'      X, K=3, centers=X[sample(1:300,3), ],
#'      tolerance=1e-3, max_iter=100)
#'
#' ### plotting the clusters ###
#'
#' oldpar = par(mfrow = c(1, 2))
#'
#' plot(X, type = 'n', main = 'ktaucenters (Robust) \n outliers: solid black dots')
#' points(X[sal$cluster == 1, ], col = 2);
#' points(X[sal$cluster == 2, ], col = 3);
#' points(X[sal$cluster == 3, ], col = 4)
#' points(X[sal$outliers, 1], X[sal$outliers, 2], pch = 19)
#'
#' ### Applying a classical (non Robust) algortihm ###
#' sal <- kmeans(X, centers = 3, n_runs = 100)
#'
#' ### plotting the clusters ###
#' plot(X, type = 'n', main = 'kmeans (Classical)')
#' points(X[sal$cluster == 1, ], col = 2);
#' points(X[sal$cluster == 2, ], col = 3);
#' points(X[sal$cluster == 3, ], col = 4)
#'
#' par(oldpar)
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198.
#'
#' @importFrom stats kmeans dist qchisq
#' @export
ktaucenters <- function(x,
                        centers,
                        iter_max = 100L,
                        tolerance = 1e-6,
                        n_runs = 1L,
                        init_centers = list(quote(init_kmeans), quote(init_robin)),
                        flag_outliers = outliers_tau_cutoff(0.999,
                                                            0.5)) {
    # Parameters check                                                    
    x <- as.matrix(x)
    n_clusters <- ifelse(is.list(centers), length(centers), centers)
    
    # Set up center initialization
    add_init_custom = NULL
    if (is.list(centers))
        add_init_custom = quote(init_custom)
    add_init_random = replicate(n_runs, init_random, simplify = FALSE)
    init_centers = append(init_centers, c(add_init_random, add_init_custom), 1)
    
    # Runs
    best_tau = Inf
    for (iter in seq_along(init_centers)) {
        current_run = ktaucenters_run(x, 
                                      eval(init_centers[[iter]])(x, n_clusters),
                                      tolerance, iter_max)
        if (current_run$tau < best_tau) {
            best_run = current_run
            best_tau = best_run$tau
        }
    }
    
    # Outlier detection
    best_run = flag_outliers(best_run)
    
    return(best_run)
}