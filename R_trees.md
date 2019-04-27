start with data prep
send a link

# R and RStudio
R is a free software environment for statistical computing and graphics, and can be downloaded [here](https://cran.uib.no/). RStudio is an integrated development environment for R, which in short makes R easier to use. Rstudio can be downloaded [here](https://www.rstudio.com/products/rstudio/download/).

For this session, you will need both R and RStudio.

# Packages
In R, the fundamental unit of shareable code is a package. There is a myriad of packages available for download, and the ones published on the Comprehensive R Archive Network [CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html) are trustworthy and of high quality. Packages can be created by anyone, and can also be hosted on GitHub. 

To install packages from CRAN, use the following command:
```{R}
install.packages("package_name", dependencies = TRUE)
```

For this session, you will need the following packages:
- `ape` : [webpage](http://ape-package.ird.fr/), [Manual](https://cran.r-project.org/web/packages/ape/ape.pdf)
- `ggtree` : [vignettes](http://www.bioconductor.org/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html)
- `cluster` : [documentation](https://cran.r-project.org/web/packages/cluster/cluster.pdf)
- `tibble` : [vignette](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html)
- `dplyr` : [website](https://dplyr.tidyverse.org/)

To install packages from GitHub, you will need the following package:
- `devtools`: [github](https://github.com/r-lib/devtools)

To install from github, you need to use the following functions:
```{R}
library(devtools)
install_github("hkaspersen/distanceR")
```

# Training data
The data files used in this session can be downloaded below (left click "download as..." on the link).

[cgMLST data](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/cgMLST.tsv)

[Tree metadata](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/tree_metadata.txt)

[Tree heatmap data](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/tree_heatmap_data.txt)

# Calculating distances from cgMLST data
Based on the allele data in the cgMLST file, one can calculate the percentage of similar "labels" in a pairwise fashion. In other words, the first sample in the data is compared to the second, and the percentage of similar labels is calculated (excluding the loci if one label is missing). The first sample is then compared to the next sample, until all samples have been compared to all. The result of this calculation is a dissimilarity matrix, which is N x N big. This data can then be clustered, and a tree object can be created from this data.

To do this in R, one has to import the cgMLST data first:
```{R}
library(tibble)
library(dplyr)

cgMLST_data <- read.table("cgMLST.tsv", sep = "\t", colclasses = "factor") %>%
  na_if("0") %>%
  column_to_rownames("FILE")
```

To calculate distances from the data:
```{R}
library(cluster)
library(ape)

distances <- as.phylo(hclust(daisy(cgMLST_data, metric = "gower"), method = "average"))
```

**Important: Labeling**:
- Labels of isolates in pairwise distance matrix and labels for your isolates-labels in your metadata file
have to be **perfectly identical** . This identity allow us to create a link (as linking database by an identical key)
between the data representing the tree and the metadata. If not identical, you won't be able to annotate the tree graphic.




A good introduction for learning how to use trees data in R: [data Integration, Manipulation and Visualisation of phylogenetic trees]
(We used it a lot to make this page.)

# Create a tree - and simple plot (just labels)

```{R}
#load the remaining necessary packages
library("packagename") #you need to load: ape, ggtree, distanceR
#
```


# Importing an existing tree -> Eve
It is possible to import trees that were created which clustering/phylogenetic softwares.
The [`treeio` package](https://bioconductor.org/packages/release/bioc/html/treeio.html), that has been developped from parsing different trees format into R.
And have been designed to import both tree and metadata associated with the tree.

- [ ] can we check if treeio is automatically imported with ggtree?

As written in [data Integration, Manipulation and Visualisation of phylogenetic trees]
Treeio supports importation of:
> - `Newick`, `Nexus`, `New Hampshire eXtended format (NHX)`, `jplace` and `Phylip`
  - outputs from programs: `BEAST`, `EPA`, `HyPhy`, `MrBayes`, `PAML`, `PHYLDOG`, `pplacer`, `r8s`, `RAxML` and `RevBayes`.

Example import:
```
mytree <- read.newick("file_name") #import the tree in R and assign it to an object of class phylo, with name mytree
plot(mytree) # an very simple way to check that your import succeded (but do not expect it to look nice)
```

NB: Treeio also offers usefull functions such as merging trees from different sources and exporting trees
Example export:

```{R}
write.beast(tree_object_in_R, file ="export_file.nex")
```


# Adding Annotations (decorations) and linking the necessary metadata


```{R code}

```

collors annotations heatmap?

# Going further
### other tree related packages but not used R

`tidytree` , `treeio` -> some functions described in [data Integration, Manipulation and Visualisation of phylogenetic trees]
phangorm


gheatmap -> eve check
> if use gheat (heatmap) -> note row.names have to be exacly the same as isolate names in tree/metadata (in col1) - in heatmap the order of rows does not need to be ordered

### links

[Analysis of Phylogenetics and Evolution with R](http://ape-package.ird.fr/APER.html)

Many R packages can be found here ...

[Bioconductor package portal](https://bioconductor.org/packages/release/BiocViews.html#___Software)

[CRAN Task view](https://cran.r-project.org/) R packages sorted by themes

There is no overview of all the R-packages hosted on GitHub, so we need some `googling` to find them.


[data Integration, Manipulation and Visualisation of phylogenetic trees]:https://yulab-smu.github.io/treedata-book/index.html
[Haukon's github]:https://github.com/hkaspersen/distanceR.git
