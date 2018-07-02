## Function gillespie: SSA function (grid version)
## Performs one run of the homogenous SSA, outputting population sizes at different times.
## We refer to this as a `grid` version because it does not record all events, rather it
## records the state of the system across a vector of monotonically increasing time points.
## The vector must begin with 0, but can otherwise be any size and need not be regular.
## The grid of times does not affect the accuracy of the simulation and can be easier to
## work with than a function that records all events.
##
## Function `gillespie` requires the following arguments:
## `init`=Array containing initial conditions (e.g. initial population size)
## `times`=Vector of times over which to record the system state
## `param`=Array containing all parameters and functions required for calculating intensities
## `inten`=Function which returns intensities of point process(es) given `param`
## `pproc`=Array containing state changes caused by point process(es) [same order as `inten`]]
## `hpp`=Function for sampling inter-event times
##
gillespie <- function(init, times, param, inten, pproc, hpp){
    if(length(times) == 0){
        stop("No time points provided in `times`")
    } else if(times[1] != 0){
        stop("First time point is not 0")
    } else {
        tottime <- times[1]
        tinc <- length(times)
        N <- init
        results <- matrix(nrow = tinc, ncol = length(init))
        results[1, ] <- N
        for(i in 2:tinc){
            results[i, ] <- results[i - 1, ]
            while(tottime <= times[i]){
                intentemp <- inten(tottime, N, param)
                if(all(intentemp == 0)){
                    results[i:tinc, ] <- N
                    i <- tinc
                    warning("Exiting with all intensities equal to 0")
                    break
                } else if(min(intentemp) < 0) {
                    results[i:tinc, ] <- NA
                    i <- tinc
                    warning("Exiting with intensity less than 0")
                    break
                } else {
                    tau <- hpp(intentemp)
                    tottime <- tottime + tau
                    which.pproc <- sample(1:nrow(pproc),
                                          size = 1,
                                          prob = intentemp)
                    if(tottime > times[i]){
                        results[i, ] <- N
                        N <- N + pproc[which.pproc, ]
                        break
                    } else {
                        N <- N + pproc[which.pproc, ]
                    }
                }
            }
            if(i == tinc) break
        }
        cbind(times, results)
    }
}

## Function hpp: inter-event time sampling function (homogenous Poisson process)
## Argument `intentemp` is the output of the intensity function at a given time
##
hpp <- function(intentemp){
    -(1 / (sum(intentemp))) * log(1 - runif(1))
}

## Function nhpp: inter-event time sampling function (nonhomogenous Poisson process)
## Generates inter-event times for a continuous nonhomogenous process using inverse transform.
##
## Function `nhpp` requires the following arguments (provided by `gillespie_plus` below):
## `tottime`=Current time of simulation
## `N`=Current system state
## `param`=Array containing parameters and functions required for calculating intensities
## `inten`=Function which returns intensities of point process(es) given `param`
## `timeleft`=Time left in the simulation
## Subdivisions and tolerances (in integrate) may be altered to increase speed or precision.
nhpp <- function(tottime, N, param, inten, timeleft){
    tryCatch(uniroot(function(X, Y) {
        1 - exp(-integrate(Vectorize(function(X){
            sum(inten(tottime + X, N, param))}), 0, X, subdivisions = 200)$value) - Y},
                     lower = 0, upper = timeleft, tol = 1e-5, Y = runif(1))$root,
             error = function(c) timeleft + 1)
}
##
## Function gillespie_plus: SSA+ function (grid version)
## Performs a run of the nonhomogenous SSA, outputting population sizes at different times.
## The only major difference between this and `gillespie` is the use of `nhpp` over `hpp`.
##
gillespie_plus <- function(init, times, param, inten, pproc, nhpp){
    if(length(times) == 0){
        stop("No time points provided in `times`")
    } else if(times[1] != 0){
        stop("First time point is not 0")
    } else{
        tottime <- times[1]
        tinc <- length(times)
        N <- init                           
        results <- matrix(nrow = tinc, ncol = length(init))
        results[1, ] <- N
        for(i in 2:tinc){
            results[i, ] <- results[i - 1, ]
            while(tottime <= times[i]){
                intentemp <- inten(tottime, N, param)
                if(all(intentemp == 0)){
                    results[i:tinc, ] <- N
                    i <- tinc
                    warning("Exiting with all intensities equal to 0")
                    break
                } else if(min(intentemp) < 0) {
                    results[i:tinc, ] <- NA
                    i <- tinc
                    warning("Exiting with intensity less than 0")
                    break
                } else {
                    tau <- nhpp(tottime, N, param, inten, (times[tinc] - tottime))
                    tottime <- tottime + tau
                    intentemp <- inten(tottime, N, param) # recalculate for new tottime
                    which.pproc <- sample(1:nrow(pproc),
                                      size = 1,
                                      prob = intentemp)
                    if(tottime > times[i]){
                        results[i, ] <- N
                        N <- N + pproc[which.pproc, ]
                        break
                    } else {
                        N <- N + pproc[which.pproc, ]
                    }
                }
            }
            if(i == tinc) break
        }
        cbind(times, results)
    }
}

####### Example application of `gillespie` and `gillespie_plus` [exponential growth]
## Initial population size
init <- 100
## Times
times <- seq(0, 10, 1)
## Terms for `param` list
### Environment function (constant)
constant <- function(t) 1
### Demographic functions
exponentialfun <- list(
    births = function(b, envres) b * envres,
    deaths = function(d, envres) d
)
### Combine parameters and demographic functions into list
param <- list(b = .03, d = .027, exponentialfun, env = constant)
param <- unlist(param) # collapse to one-dimensional list
## Intensity functions which use the demographic functions (above)
inten <- function(t, X, param){
    with(as.list(c(param)),{
        bint <- X * births(b, env(t))
        dint <- X * deaths(d, env(t))
        c(bint, dint)})
}
## State changes caused by the point processes (birth = N + 1, death = N - 1)
pproc <- matrix(c(1,-1), nrow = 2)
## Run 1 simulation of gillespie and gillespie_plus
set.seed(20170915)
gillespie(init, times, param, inten, pproc, hpp)
set.seed(20170915)
gillespie_plus(init, times, param, inten, pproc, nhpp)
