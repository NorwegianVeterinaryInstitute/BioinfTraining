# 1. What is [PopPUNK](https://genome.cshlp.org/content/29/2/304) : Population Partitioning Using Nucleotide K-mers
## 1.1 Summary - Overview

[Manual](https://poppunk.readthedocs.io/)
[PopPUNK Documentation](https://poppunk.readthedocs.io/en/latest/index.html)
[Presentation](https://docs.google.com/presentation/d/1StmmM02lSpFPdevQT3iDB3BAKRoMC8Q4NtJ4nzu7MdY/edit?usp=sharing)

- whole genome (core + assessory) population analysis/clustering
- distinction between isolates : uses k-mer of different length: `mash` to find core and accessory distances between isolates (**pairwise**)
- the distribution of those distances is used to discriminate clusters (defined as strains) of closely related isolates (similarity:both core and accessory)
- clustering of newly added isolates: EXTENDABLE4 - without the need of reanalizing all samples
- maintenance free and auto-reduce database
> aimed for consistent naming clusters between studies, outbreak detection in minutes
- analyse up to 10^4 samples in single step - possible to add new samples
-
# 2. How to use PopPUNK

## 2.1 Create datadase

- create a list of your assemblies and other sequences you want to include: you need to have the path of the file included, so you need to create the list from the folder you will run PopPUNK: `ls <path>/*fasta > reference_list.txt`
-
-


poppunk --fit-model --distances <folder/folder.dists> --ref-db <folder_to_db/model_db> --output <output_folder> --full-db --K <nb_blubs>

## 2.2 Fitting a model: what to look for



<p align="center">
<img width="500" src="./figures/refine_poppunk.png">
</p>

> trick for seing plot : test! - [ ] use move with bondary with 0

> kernet density estimate - > can be used to identify outliers and contamination (program to remove those isolates from DB) > - [ ] where?



- how to do: accessory, core, and combined clusters (see sup material) - [ ]

### 2.2.1 2D GMM
> equal likelihood contours and decision boundary (within and between cluster-assignments )

### 2.2.2 HDBSCAN -
`--dbscan`

### 2.2.3 Evaluate the chosen model
- silhouette index = measure of how similar an object is to its own cluster (cohesion) compared to other clusters (separation) -> must be close to 1
- network density must be low (mean good separation between clusters: fewer within than between strains links)

> **important to look at the network**

#### 2.2.4 Only if the model is really bad: (otherwise -> go to refine model)

- possibility to optimize using core distance only (vertical boundary) or accessory only (horizontal boundary)
> when?  if core and accessory genomes have independent evolution histories - ex. lots of recombination (blurs in blobs), prophages...
> NB: depending on level problem - can be sufficient to only refine model

- **modify sketch size** (increase sensibility detection SNPs also increases running time) -> finer differences in %pi (down to single SNPs) --option -[ ]

## 2.3. refining the model `--refine-model`

`poppunk --refine-model --distances <folder/folder.dists> --ref-db <folder_to_db/model_db> --output <output_folder> --full-db` +  For 2D GMM `--K <nb_blubs>` OR For HDBSCAN `--dbscan`

Common optional options: HDBSCAN & 2D GMM
`--pos-shift <POS_SHIFT>` (away from origin)
`--neg-shift <NEG_SHIFT>` (towards origin)
`--manual-start <filename>`
`--indiv-refine` allow these boundaries to be placed independently on core/accessory

<p align="center">
> need to create a triangular⁵ boundary - move forward and backward FROM starting point (range)
```
mean0 (x,y) #for within strains blob
mean1 (x,y) #for between strains blob
start (x,y) #starting point to move boundary to
```
</p>


example: `poppunk --refine-model --distances lm_example/lm_example.dists --ref-db lm_example_2DGMM --output refine_gmm --full-db --neg-shift 0.1`
to obtain the same refinement as published in their article



NB: **Better to redo all the steps in the same folder - now that we now the best fit - to facilitate further use**

## 2.4. Simplification database (optional)
When model is good -> then we can stop using `--full-db` option (but not compulsory)

## 2.5 Visualisation resuts
> output can be made either from --model-fit or --refine-model
`--microreact`	Generate output files for microreact visualisation
`--cytoscape	`Generate network output files for Cytoscape
`--phandango`	Generate phylogeny and TSV for Phandango visualisation
`--grapetree`	Generate phylogeny and CSV for grapetree visualisation
`--rapidnj` RAPIDNJ (Path to rapidNJ binary to build NJ tree for Microreact)
`--perplexity `PERPLEXITY (Perplexity used to calculate t-SNE projection (with –microreact) [default=20.0])
`--info-csv` INFO_CSV (additional metadata: Epidemiological information CSV formatted for microreact (can be used with other outputs))


### 2.5.1 Microreact

### 2.5.2.GrapeTree

### 2.5.3. Phandango

### 2.5.4 Cytoscape


## 2.6 Adding new sequences = assigning queries
> addition to the reference network
1. pairwise distances are calculated
2. added as nodes in network to clusters
3. clusters name **DO NOT CHANGE** - unless merged -> then both labels displayed on merged cluster

`poppunk --assign-query --ref-db <database> --q-files <query_list.txt> --output <strain_query> --threads <3> --update-db`

> optional `--model-dir <directory>`  if fitted model is in a separated directory
> previsous clustering/network: `--previous-clustering`

> **OOPS!** for further adding queries: need to use the new database: stored in strain-query `--ref-db <strain_query>`

### Queries using core or accessory only;

--ref-db <refine-model-db> `--core-only`
--ref-db <refine-model-db> `--accessory-only `


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX





# Poppunk commands


# Trics for difficult species:
- Low diversity -> a and %pi not necessary correlated - analyse independently (see supplement)
  - genome wide association study: genes responsible differences in a and %pi clusters (pyseer) -> used a-cluster as phenotype


# My Points
- [ ] as use random sampling of distances to create network -> maybe important not to unbalanced sampling to create model -> otherwise model risicate to be pushed toward most common lineage (specially if expect new samples different) - same pb as usual ...
- [ ] how to extracti within and between-cluster SNPs distances -> %pi??
- [ ] understand this sketch size -> per scaffold or over whole genome - must be high!! -> the scetch size -> determines what will be most computationnaly intensive

- within cluster distance: appear under-estimate - but between comparable to RhierBAPS (maybe depend on model used)

# ## 2.opening PopPUNK (how does it functions)

- p_random : maximum probability of random matches kmerse (set to 5% ?? check) -> determines k_min
- kmax set at 29 by default (can increase to max mash is able to support = - [ ] find )

## 2.1 Models: classify which distances pairs (%pi and a) are within the same cluster

### 2.1.1 Model to fit -> Two-dimensional Gaussian mixture model (2D GMM)
Fits random subsample of up to 10⁵ distance pairs
scikit-learn 0.19
Dirichlet Process prior on weights (nearly 0 to define within-cluster distances)
best final likelihood from five k-means initial starts
K: maximum allowed of mixture components (default = 2)
Distances classified with fitted model

### 2.1.2 Model to fit -> HDBSCAN (option dbscan)
Classify a subsample of 10⁵ points (distance pairs) using:
Boruvka ball tree Algorithm
Iterations - progressive reductions -> minimum of samples required to initiate the seqrch for a cluster
 - -> how conservative clustering is
 - -> minimum cluster size (threshold nb points(distances a cluster must contain)
 - -> extends points in cluster closest to origin (represent within strains) - should not overlap with between strains clusters

### 2.2 Networks constructed:
networkx v2.1 -> undirected graph with unweighted egdes -> population clusters
extraction:
- connected components (nodes = samples, connections = pairwise distances)
- ordered by nb isolates (largest to smallest)
Evaluation netwok structure: ns = transitivity (1 - density)
- [ ] transitivity?
- (1-density) -> to subdivise the population

ns  > 0.8 - 1 good fit

database simplification: randomly select one isolate - per cluster (contected to all other members clusters)

### 2.3 Refinement distance classification
Models -> treat clusters symmetrically (force clusters same distance/structure)
Refinement -> delimit precisely range %pi - a distances trated as within strain links - to maximize ns
  - line between means of within and between strain clusters
  - if not correct fit -> provide within and between stran cluster means manually - move boundary starting point (range of move) -> maximize ns (40 equaly spaced points over the allowed range)

# Summary - Files produced

| File extension                   | what                                                                                                                     | when                        |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | --------------------------- |
| *.search.out                     | pairwise core & accessory distances                                                                                      |                             |
| *graph.gpickle                   | network used to predict clusters                                                                                         |                             |
| *DPGMM_fit.png                   | scatter plot of all distances, and mixture model fit and assignment                                                      |                             |
| *DPGMM_fit_contours.png          | contours of likelihood function fitted to data                                                                           |                             |
| *distanceDistribution.png        | scatter plot of the distance distribution fitted by the model + kernel-density estimate                                  |                             |
| *.csv                            | isolate names and the cluster assigned                                                                                   |                             |
| *(db).png                        | unclustered distribution of distances used in the fit (subsampled from total)                                            |                             |
| *.npz                            | save fit parameters                                                                                                      |                             |
| *refs                            | representative references in the new database                                                                            | (unless --full-db was used) |
| *dbscan.png                      | scatter plot of all distances, and DBSCAN assignment.                                                                    | --dbscan                    |
| *external_clusters.csv           | CSV file relating the samples to previous clusters provided in the input CSV.                                            | --external-clustering       |
| *core_dists.csv                  | matrix of pairwise core distances                                                                                        | --microreact                |
| *acc_dists.csv                   | matrix of pairwise accessory distances                                                                                   | --microreact                |
| *core_NJ_microreact.nwk          | neighbour joining tree using core distances (for microreact)                                                             | --microreact                |
| *perplexity5.0_accessory_tsne.dot | t-SNE embedding of accessory distances at given perplexity                                                               | --microreact                |
| *microreact_clusters.csv         | cluster assignments plus any epi data added with the --info-csv option (for microreact)                                  | --microreact                |
| *cytoscape.csv                   | cluster assignments plus any epi data added with the --info-csv option (for cytoscape)                                   | --cytoscape                 |
| *cytoscape.graphml               | XML representation of resulting network (for cytoscape)                                                                  | --cytoscape                 |
| *refined_fit.png                 | plot of the new linear boundary, and core and accessory distances coloured by assignment to either side of this boundary | --fit-model                 |
| *refined_fit.npz                | **The saved parameters of the refined fit.**                                                                             | --fit-model                 |
| *clusters.csv + .gpickle       |  for core and accessory                                                                      |--fit-model  --indiv-refine |
