nhppCDF <- function(tottime, N, param, inten){
    1 - exp(-integrate(Vectorize(function(X){
        sum(inten(X, N, param))}), lower = 0, upper = tottime)$value)
}
