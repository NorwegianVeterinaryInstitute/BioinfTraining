oA good introduction for learning how to use trees data in R: [data Integration, Manipulation and Visualisation of phylogenetic trees]. (We used it a lot to make this page.)

# R and RStudio
R is a free software environment for statistical computing and graphics, and can be downloaded [here](https://cran.uib.no/). RStudio is an integrated development environment for R, which in short makes R easier to use. Rstudio can be downloaded [here](https://www.rstudio.com/products/rstudio/download/).

For this session, you will need both R and RStudio.

R is an object orientated programming language. What you need to understand for this lesson is that you can assign  with the symbol `<-` objects (dataset in form of tables, results of calculations, trees ...) into an object that you name. Example: `my_color <- "red"` and you can recall objects from memory by using their names: example type: `my_color`

# Packages
In R, the fundamental unit of shareable code is a package. There is a myriad of packages available for download, and the ones published on the Comprehensive R Archive Network [CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html) are trustworthy and of high quality. Packages can be created by anyone, and can also be hosted on GitHub. Some packages, specialized for biological analyses, are hosted by [Bioconductor](https://www.bioconductor.org/), and have their own installation method.

To install packages from CRAN, use the following command:
```{R}
install.packages("package_name", dependencies = TRUE)
```

For this session, you will need the following packages
- `ape` : [website](http://ape-package.ird.fr/), [manual](https://cran.r-project.org/web/packages/ape/ape.pdf)
- `cluster` : [documentation](https://cran.r-project.org/web/packages/cluster/cluster.pdf)
- `tibble` : [vignettes](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html) *(vignettes are nearly equivalent to manual pages.)*
- `dplyr` : [website](https://dplyr.tidyverse.org/)
- `ggtree` : [vignettes](http://www.bioconductor.org/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html). We have to install this package from GitHub (see below)
- `distanceR` Also need to be installed from GitHub: [Haukons' GitHub](https://github.com/hkaspersen/distanceR)

To install packages from GitHub, you will need the following package:
- `devtools`: [github manual](https://github.com/r-lib/devtools)

To install from github, you need to use the following functions:

```{R}
# loading the packages we will need in R environment
# NB: loading packages allows to import functions defined in to the packages into R memory.
# This makes the functions defined in those pacakges are available for use in R.

library(devtools)
install_github("hkaspersen/distanceR")
install_github("GuangchuangYu/ggtree")
```

# Training data
The data files used in this session can be downloaded below (left click "download as..." on the link).

[cgMLST data](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/cgMLST.tsv)

[Tree metadata](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/tree_metadata.txt)

[Tree heatmap data](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/tree_heatmap_data.txt)

NB: if you are working on a Linux machine, you can download those files from the command line with: `wget <url_file>`

# Calculating distances from cgMLST data
Based on the allele data in the cgMLST file, one can calculate the percentage of similar "labels" in a pairwise fashion. The labels in the cgMLST data represent the loci variant, based on the nucleotide sequence (same as regular 7 gene MLST).

In other words, the first sample in the data is compared to the second, and the percentage of similar labels is calculated (excluding the loci if one label is missing). The first sample is then compared to the next sample, until all samples have been compared to all.

The result of this calculation is represented in a dissimilarity matrix, which is N x N big. Diagonals of the matrix represent each isolate compared to iself (0 % because totally identical). This matrix is symetrical (think why...). The most dissimilar the pair will be the closest to 1 % will be the value of comparison for the pair will be.

This data can then be clustered hierachically (groupped based on dissimilarity), and a tree object can be created from this hierarchical clustering.

To do this in R, one has to import and clean the cgMLST data first:
```{R}
# loading the packages we will need in R environment

library(tibble)
library(dplyr)

# Bellow we import the data (tablular data, separator = tabulation, each column is a factor: categorical data) - %>% is a pipe simbol in R ...

cgMLST_data <- read.table("cgMLST.tsv", sep = "\t", colClasses = "factor", header = TRUE) %>%

  # change all "0" to NA
  na_if("0") %>%
  # create rownames from the values in the column "FILE"
  column_to_rownames("FILE")
```

To calculate distances from the data and create a tree:

```{R}
library(cluster)
library(ape)

tree <- as.phylo(hclust(daisy(cgMLST_data, metric = "gower"), method = "average"))
```
Check out the function information with `?daisy` and `?hclust` to see details on the metric and method arguments

To then visualize this tree, we can either use the `plot()` function from base R, or use ggtree function from ggtree package.

*NB: base R is a minimal set of packages that are loaded automatically when you install R (functions associated with the language functionning of R, basic statistical and graphical functions)*
> What are the advantages of using ggtree?
> - using base plot() function allows us to have a fast view at our tree object (check that it worked)
> - ggtree has advanced anotations functions that allow us to add layers and follows a specifig synthax that allows to do advanced graphics "automatically and easily". This special synthax is called ["grammar of graphics"](https://en.wikipedia.org/wiki/Leland_Wilkinson) which give the gg from the ggtree package.

> What are layers?:
  > -  Imagine that me made a basic graphic on a white A4 paper. We can add labels, annotations, colors, title whith layers: think of a layer as a transparent sheet where you only draw one type of information (ex: leaves labels) and you put this transparent sheet on top of your basic tree. You can now see your tree annotated with the leaves labels.

```{R}
library(ggtree)

plot(tree)
ggtree(tree)
```
In base R, the tree will have sample name labels on the tips. In ggtree, nothing was specified, so only the tree was plotted. Ggtree have multiple layouts to choose from, i.e. "circular", "rectangular", "fan", "equal_angle" etc. See the vignette on ggtree above for more information. To add more information to the tree with ggtree, we can use additional layers of information, like this:

```{R}
library(ggtree)

ggtree(tree,
       layout = "rectangular") +
     geom_tiplab()
```
Now, the isolate names have been added to the tree tip points (leaves labels). However, if you want to combine metadata such as sequence types, animal species etc. to the tree and visualize it, one first has to connect the tree and the metadata:

# Linking metadata to your tree data allows you to add annotations (decorations)

```{R}
library(ggtree)

# Import metadata
metadata <- read.table("tree_metadata.txt", sep = "\t", header = TRUE)

# Plot tree
ggtree(tree,
       layout = "rectangular") %<+% metadata +
     geom_tiplab(aes(label = ST))
```

Now, the ST column in the metadata is plotted in the tree instead of the isolate names. **Important**: For this to work, one has to make sure that the first column in the metadata file has exact matches to the tip labels in the tree object. In other words, if one sample is named "2014-01-2355" in the tree, and the same isolate is named "68-2014-01-2355" in the metadata, the information will not be plotted for that isolate. Therefore, to check the tip labels on the tree, one can type `tree$tip.labels` to see what they look like. Note that the order of the samples doesn't matter in the metadata file, just the names of the samples.
Notice that the `aes()` function was added to the `geom_tiplab()` function. This links the column "ST" in the metadata to the corresponding samples in the tree, plotting the ST where the specific sample is placed in the tree. Anything written outside the `aes()` function isn't linked to the metadata in any way. Thus, to adjust the size of the labels, one can type in `geom_tiplab(aes(label = ST), size = 3)`.

Now, with the metadata file in hand, one can add more information to the tree, for example colored tip nodes:

```{R}
# Plot tree
ggtree(tree,
       layout = "rectangular") %<+% metadata +
     geom_tiplab(aes(label = ST)) +
     geom_tippoint(aes(color = species))
```
Ggtree will here plot the animal species in the metadata as colored nodes on the tree. Since no specific palette is specified, it will use default ggplot2 coloring. If you want to specify colors, a specific palette need to be created:

```{R}
# create palette
palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43",
             "Wild bird" = "#fdae61")

# Plot tree
ggtree(tree,
       layout = "rectangular") %<+% metadata +
     geom_tiplab(aes(label = ST)) +
     geom_tippoint(aes(color = species)) +
     scale_color_manual(values = palette)
```
- [ ] scale_color_manual is a ggplot2 function -> so we need to import ggplot2

# Importing an existing tree
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
mytree <- read.newick("file_name") # import the tree in R and assign it to an object of class phylo, with name mytree
plot(mytree) # an very simple way to check that your import succeded (but do not expect it to look nice)
```

NB: Treeio also offers usefull functions such as merging trees from different sources and exporting trees
Example export:

```{R}
write.beast(tree_object_in_R, file ="export_file.nex")
```

**collors annotations heatmap?**

# Going further
### other tree related packages but not used in here

`tidytree` , `treeio` -> some functions described in [data Integration, Manipulation and Visualisation of phylogenetic trees]

`ggplot2` has some additional functions that can be used to annotate trees (ggtree package is based on ggplo2)


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
