A model of socially cued plasticity.

To run the model open run_model.R.  To explore conditions where socially cued platicity might evolve, you can change the base.prms() to explore a variety of parameters including period of environmental change, initial s and g locus values, strength of selection, cost to plasticity, selection type asymmetries in selection strength, asymmetries in selection period.   

The initialize file loads the files to run the model. It also sets a base set of parameters and creates an intitial population.  

The simulate_functions file includes functions used to run the model. This includes printing stats of the model as it runs, creating and updating an object to save model output and calculating the strength of selection on every day. It also includes run.sim which is the function that runs the simulation.   

The demography file has code that determines major demographic events for individuals including ageing, maturation with strategy and phenotype expression, death, birth, and mutation. 

For more deatils see our forthcoming publication:

Elizabeth C. Lange, Joseph Travis, Kimberly A. Hughes & Leithen M’Gonigle (Accepted). Can you trust what you see? The evolution of socially-cued anticipatory plasticity. The American Naturalist.

