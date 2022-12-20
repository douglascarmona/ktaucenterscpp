KMEANS_K_DEFAULT <- 20
ROBIN_DEFAULT_METHOD <- "density"

init_kmeans <- function(data, K) {
  ret <- kmeans(data, K, nstart = KMEANS_K_DEFAULT)$centers
  return(ret)
}

init_random <- function(data, K) {
  ret = data[sample(1:nrow(data), K), ]
  return(ret)
}

init_robin <- function(data, K) {
  ret_rob <- robinden(distance(data), K, 10)
  r_centers <- ret_rob$centers + 1
  ret <- data[r_centers, ]
  return(ret)
}

init_custom <- function(data, K) {
  ret <- K
  return(ret)
}