##########################################################################
### The code below includes a variety of functions to run the model.     #
### This includes priniting stats of the model as it runs, creating and  #
### updating an object to save model output and calculating the strength #
### of selection on every day. It also includes run.sim which is the     #
### function that runs the simulation.                                   #
##########################################################################


## print pop stats
print.stats <- function(prms, pop, day, time.start) {

  cat('\r')
  time.diff <- round(difftime(Sys.time(), time.start, units='mins'))

  case <- ifelse(prms$selection=='F', 'F', 'S')

  good <- ifelse(prms$selection=='S',
                 names(prms$s)[which.min(prms$s)],
                 names(prms$s)[which.max(prms$s)])
  bad  <- ifelse(prms$selection=='S',
                 names(prms$s)[which.max(prms$s)],
                 names(prms$s)[which.min(prms$s)])

  cat(sprintf('%d min, %2.0f%%: %s, %s=%d, %s=%2.3f, %s=%2.2f, %s>%s',
              time.diff,
              day/prms$num.days*100,
              case,
              'Mature', sum(pop[,'mature']),
              'mean(p)', mean(pop[,'loc.s']),
              'freq(a)', length(which(pop[,'mature']==1 &
                                      pop[,'phenotype']==1)) /
                         max(1, sum(pop[,'mature'])),
              good, bad))
}

## create structure to save model summary info
make.dd.save <- function(prms, period) {
  if(prms$to.save[['type']]=='day') {
    from <- max(1,prms$num.days-prms$to.save[['num']])
  }
  if(prms$to.save[['type']]=='period') {
    if(prms$num.days<(prms$to.save[['num']]*sum(period)))
      cat('num.days too small to store desired number of periods.\n')

    from <- prms$num.days - prms$to.save[['num']]*sum(period) + 1
  }
  to <- prms$num.days
  ll <- min(prms$to.save[['pts']], to-from+1)
  save.days <- round(seq(from=from, to=to, length=ll))
  
  cols <- c('day',
            'favoured',
            'mean.age',
            'num.mature',
            'num.mature.a',
            'num.mature.b',
            'freq.a',
            'freq.b',
            'freq.adapted',
            'mean.loc.s',
            'var.loc.s',
            'mean.loc.g',
            'var.loc.g')
              
  res <- matrix(NA, nrow=length(save.days), ncol=length(cols))
  colnames(res) <- cols
  res[,'day'] <- save.days
  res
}

## update above data-structure for a single day
get.pop.stats <- function(x, prms, pop) {
  x['favoured'] <- ifelse(prms$selection=='F',
                          which.max(prms$s),
                          which.min(prms$s))
  x['mean.age'] <- mean(pop[,'age'])
  x['num.mature']   <- sum(pop[,'mature'])
  x['num.mature.a'] <- sum(pop[,'mature']==1 & pop[,'phenotype']==1)
  x['num.mature.b'] <- sum(pop[,'mature']==1 & pop[,'phenotype']==2)
  x[c('freq.a', 'freq.b')] <-
    x[c('num.mature.a', 'num.mature.b')]/x['num.mature']
  x['freq.adapted'] <- x[c('freq.a','freq.b')][x['favoured']]
  x['mean.loc.s'] <- mean(pop[,'loc.s'])
  x['var.loc.s']  <- var(pop[,'loc.s'])
  x['mean.loc.g'] <- mean(pop[,'loc.g'])
  x['var.loc.g']  <- var(pop[,'loc.g'])
  x
}

## calculate selection strength
calc.s <- function(prms, day, pop) {

  ## step-function selection regime
  if(prms$regime=='step') { 
    favoured <-
      ((day-1) %% (prms$period[1]+prms$period[2]) < prms$period[1])*1
    return(prms$str.sel * (-1)^(favoured)) 
  }


}

## run the model and vary selection type between selection and fecundity
run.sim <- function(prms, print=TRUE) {

  ## start time
  time.start <- Sys.time()
  ## array to save model output in
  dd.save <- make.dd.save(prms, period)
  ## initial population
  pop <- make.pop(prms)
  
  for(day in 1:prms$num.days) {
    ## when sel.coeff>0, a is favoured and
    ## when sel.coeff<0, b is favoured
    sel.coeff <- calc.s(prms, day, pop) 
    pos <- sel.coeff>0
    
    if(prms$selection=='F')
      prms$s <- c(a=1+sel.coeff*pos + prms$asym.sel , b=1-sel.coeff*(!pos)) #a favored in asym
    
    if(prms$selection=='S')
      prms$s <- c(a=1-sel.coeff*(!pos), b=1+sel.coeff*pos + prms$asym.sel) #a favored in asym
    
    ## next generation
    pop <- next.day(prms, pop)
    
    ## print progress to console
    if(print & day%%prms$print.every==0)
      print.stats(prms, pop, day, time.start)

    ## save output
    if(day %in% dd.save[,'day']) {
      ind <- match(day,dd.save[,'day'])
      dd.save[ind,] <- get.pop.stats(dd.save[ind,], prms, pop)
    }
  }
  dd.save
}
