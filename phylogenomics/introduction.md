# Preampule
We will present here different methods of making **"not wrong trees"**.

We will focus on tree reconstruction based on characters
(as many of you already have some experience with making trees based on distance)

Some methods and steps are shared between the different models used to reconstruct trees (ie. matrices models of
character changes, robustness evaluation).

I choose to embed the description of those shared method in each document describing the tree building method
- so you can have an overview of all steps for each tree building method independently.
This means that you will have the complete information for each method (without the need to look at the other method).

Those common parts will be delimited with [?? icone]

But if you are curious and look at different methods: then some parts will be totally identical.

I will heavily rely on worflows, so we can recover the logical and chronological steps of what we have to pay attention to.

> *"Phylogenetic analysis is the most widely used tool to estimate biological relationships, but it is often not the most appropriate tool to estimate relationships among bacterial genomes because, as the result of frequennt genetic exchange, different parts of the genome can have different evolutionary histories" [ksnp2 paper citat] * - well formulated - Assumes: character chared by pairs of individuals by descent - vertical not horizontal -> deviations = homoplasies (convergence, hztal transfer osv)
>   - check your assumptions: reasonable? - identify and remove those htzal parts
> - correction methds ? ie include different histories -> check ClonalFrameML - account for recombination , others ?
> - use a different method if not adapted or complementary ...w

# Warnings
- production of an optimal tree: no deterministic -> NP = polynomical time problem
(not exact solution in a reasonable amount of time)

# Tree building based on characters

![Character_tree_workflow]

common methods
Robustness

> check -[ ]
  >ML -> rate variation among branches
  > parsimony: rates sites - different rates OVER time (not branches)
  > Kuhner felsenstien 1994 - comparison phylogeny algorithms under equal and unequal evolutionary rates - mol biology and evolution
  > Kalaczkowski 2009 - Long-branch attraction bias and inconsistency in Bayesian phylogenetics plos one e7891
  > Kalaczkowski 2004 performance of maximum parsimony and Likelihood phylogenetics when evolution is heterogeneous nature
  > site specific changes in evolutionary rates over time = heterotachy (non-parametric estimation of trees : parsimony more accurate than ML) -> bacteria likely fast influence -> prefer ML (appears parsinony better for virus - explanation ksnp2)
  > minimum spanning trees -> other representations (alternatives - going further, networks osv..) -> assumes id by state (not descent - > only similarity) -
  > phylogenetic tree comparison (ie Compare Trees)
  > ? correlation between phylog trees -> and g. data after more...

## Parsimony

![Parsimony]

## Maximum Likelihood

![ML_workflow]

## Bayesian

![Bayesian_workflow]
