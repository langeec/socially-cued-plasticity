##########################################################################
### The code below loads the files to run the model. It also sets a      #
### base set of parameters and creates an intitial population.           #
##########################################################################



## load source functions
library(parallel)
library(ggplot2)
source('src/demography.R')
source('src/simulate_functions.R')

## basic parameter set (for quick loading and consistency)
base.prms <- function(K=1e3, #population size
                      num.die=5,  #number of individuals that die each day
                      ## life history
                      init.s=0, #initial value of S-trait
                      init.t=0.5, #fixed value of threshold
                      init.g=0.5, #initial value of the G-trait
                      age.mature=1e2, #age at maturation
                      alt.strategy='RP', #alternative strategy (RP for random plasticity, G for genetic, BH for bet hedging)
                      ## selection model and strength
                      selection='S', #type of selection (S for survival or F for fecundity)
                      str.sel=0.5, #strength of selection (s)
                      period=c(1e5,1e5), #period of environmental change
                      asym.days=1, #stength of asymmetry in period, 1 means periods are symmetrical
                      asym.sel=0,#strength of asymmetry in selection
                      scap.cost=0.1, #cost to plasticity
                      ## mutation
                      mut.s=TRUE, #TRUE means that the S-trait mutates
                      mut.sd=0.001, #mutational variance for the S-trait
                      s.binary=FALSE, #TRUE means that the S-trait is binary (e.g., can only be 0 or 1)
                      mut.g=TRUE, #TRUE means that the G-trait mutates
                      mut.sd.g=0.001, #mutational variance for the G-trait
                      ## simulation parameters
                      regime='step', #selection between phenotypes changes in stepwise fashion
                      num.days=1e7, #length of model run
                      to.save=list(type='day', 1e4, 1e3), 
                      print.every=1) {

  inputs <- lapply(as.list(match.call())[-1], eval.parent)
  inputs.f <- lapply(formals(), eval)
  for(v in names(inputs.f))
    if(!(v %in% names(inputs)))
      inputs <- append(inputs, inputs.f[v]) ## add formal args
  inputs[order(names(inputs))] ## order for consistency
}

## create an initial population
make.pop <- function(prms) {
  
  ## age
  age <- sample(1:(2*prms$age.mature), prms$K, replace=T)
  ## mature
  mature <- rep(0, prms$K)
  mature[age>=prms$age.mature] <- 1
  ## phenotype
  phenotype <- sample(x=c(1,2), size=prms$K, replace=TRUE) #1=a, 2=b phenotype

  as.matrix(data.frame(age=age,
                       mature=mature,
                       phenotype=phenotype,
                       loc.s=rep(prms$init.s,prms$K),  
                       loc.g=rep(prms$init.g,prms$K),
                       mom.pheno=rep(NA,prms$K)))
}
