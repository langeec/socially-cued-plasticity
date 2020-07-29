##########################################################################
### This file includes code that determines major demographic events     #
### for individuals including ageing, maturation with strategy and       #
### phenotype expression, death, birth, and mutation.                    #
##########################################################################


## mutation for when s.binary=FALSE (called from within birth)
mutate <- function(x, sigma) {
  ## add mutational effect
  x.new <- x + rnorm(n=length(x),mean=0,sd=sigma)
  ## constrain traits to be between 0 and 1
  if(any(x.new>1 | x.new<0)) {
    x.new <- pmin(x.new, 1)
    x.new <- pmax(x.new, 0)
  }
  x.new
}


## mutation for s.binary=TRUE (called from within birth)
mutate.bin <- function(x, sigma) {
  ## figure out who mutates
  mutants <- which(rbinom(n=length(x),1,prob=sigma)==1)
  if(length(mutants)==0) return(x) ## no mutants
  
  x.new <- x ## create new vector of values
  x.new[mutants] <- (x[mutants]+1)%%2 ## change 0 to 1 and 1 to 0 for mutants
  x.new
}

## iterate one day
next.day <- function(prms, pop) {

  ## individuals get older
  pop[,'age'] <- pop[,'age']+1

  ## maturation
  pop[pop[,'age']==prms$age.mature,'mature'] <- 1
  
  ## 1. calculate threshold (frequency of phenotype a)
  mature <- which(pop[,'mature']==1)
  threshold <- sum(pop[mature,'phenotype']==1) / length(mature)

  ## 2. calculate phenotype-dependent fitness
  s.phenotype <- prms$s[pop[mature,'phenotype']]
  
  ## 3. select individuals to die (only mature individuals die)
  if(prms$selection=='S') {
    ## calculate cost under survival model
    cost <- 1/(1-prms$scap.cost*pop[mature,'loc.s'])
    dies <- sample(mature, size=prms$num.die, prob=cost*s.phenotype)
  } else {
    dies <- sample(mature, size=prms$num.die)
  }
  ## set fitness to zero for dead individuals, so that they cannot be
  ## selected as parents
  s.phenotype[match(dies, mature)] <- 0
  
  ## 4. select parents
  if(prms$selection=='F') {
    ## calculate cost under fecundity model
    cost <- 1-prms$scap.cost*pop[mature,'loc.s']
    parent <- sample(mature, size=prms$num.die, replace=TRUE,
                     prob=cost*s.phenotype)
  } else {
    parent <- sample(mature, size=prms$num.die, replace=TRUE)
  }
  
  ## 5. make offspring
  offspring <- pop[parent,,drop=FALSE]
  offspring[,c('age','mature')] <- 0
  offspring[,'mom.pheno'] <- pop[parent,,drop=FALSE][,'phenotype']
  

  ## mutate traits
  ##mutate s trait
  if(prms$mut.s) {
    
    if(prms$s.binary==FALSE) 
      offspring[,'loc.s'] <- mutate(x=offspring[,'loc.s'],
                                    sigma=prms$mut.sd)
    
    if(prms$s.binary==TRUE)
      offspring[,'loc.s'] <- mutate.bin(x=offspring[,'loc.s'],
                                        sigma=prms$mut.sd)
  }
  
  #mutate g trait
  if(prms$mut.g)
    offspring[,'loc.g'] <-
    mutate(x=offspring[,'loc.g'], sigma=prms$mut.sd.g)
  
  
  ## figure out who matures via SCAP
  ## non-binary s trait
  if(prms$s.binary==FALSE)
    AP <- offspring[,'loc.s'] > runif(nrow(offspring))  ##individuals using social cues for phenotype expression
  
  ## binary s trait
  if(prms$s.binary==TRUE)
    AP <- offspring[,'loc.s']==1  ##individuals using social cues for phenotype expression
  
  num.AP <- sum(AP)
  num.not.AP <- sum(!AP)

  ## determine phenotypes to take upon maturation
  ##individuals using social cues for phenotype expression
  if(num.AP>0)
    offspring[AP, 'phenotype'] <- ifelse(0.5 <= threshold, 1, 2) 

  ##individuals using the alternative strategy
  if(num.not.AP>0) {
    
    ## a. random plasticity
    if(prms$alt.strategy=='RP')
      offspring[!AP,'phenotype'] <- sample(1:2, sum(!AP), replace=TRUE)
    
    ## b. genetic strategy
    ## Probaility of phenotype a (#1) is equal to loc.g, probability of phenotype b (#2) is equal to 1-loc.g 
    if(prms$alt.strategy=='G')
      offspring[!AP,'phenotype'] <- (offspring[!AP,'loc.g']<runif(num.not.AP))+1
    
    ## c. individuals using bet hedging
    if(prms$alt.strategy=='BH'){
      ## figure out which offspring are keeping mom's phenotype
      same <- offspring[!AP,'loc.g'] > runif(num.not.AP)
      offspring[!AP,'phenotype'][same] <- offspring[!AP,'mom.pheno'][same]
      
      ## and which offspring are not keeping mom's phenotype
      offspring[!AP,'phenotype'][!same] <- (offspring[!AP,'mom.pheno'][!same]==1)+1
    }
    }
  
  
    ## replace dead individuals with offspring and return population
  pop[dies,] <- offspring
  pop
}
