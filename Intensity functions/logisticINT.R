logisticINT <- function(t, X, param){
    with(as.list(c(param)),{
        bint <- X * (births1(b1, env(t)))
        dint <- X * (deaths1(d1, env(t))) + (X ^ 2) * (deaths2(d2, env(t)))
        c(bint, dint)})
}
