exponentialINT <- function(t, X, param){
    with(as.list(c(param)),{
        bint <- X * births(b, env(t))
        dint <- X * deaths(d, env(t))
        c(bint, dint)})
}
