##########################################################################
### The code below runs the model. We changed the base.prms to explore   #
### variety of parameters including period of environmental change,      #
### initial p, strength of selection, cost to plasticity, selection type #
### asymmetries in selection strength, asymmetries in selection period.  #
### We did so by changing the prms statement to update the base.prms     #
##########################################################################

rm(list=ls())
source('src/initialize.R')  ##sources files necesary to run the model. 
prms <- base.prms()         ##allows explorated of a variety of parameters by changing the parameter values in initialize.R
res <- run.sim(prms=prms, print=T) ##runs the model and saves output into res object
