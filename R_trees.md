start with data prep
send a link

# What is R short introduction
R and R studio
console

# Training data :
E.coli cgMLST data obtained with
The data can be downloaded from Abel: `<path>`

Data information:
- pairwise distance matrix : Distance Metrics = "gowermetric"

**Important: Labeling**:
- Labels of isolates in pairwise distance matrix and labels for your isolates-labels in your metadata file
have to be **perfectly identical** . This identity allow us to create a link (as linking database by an identical key)
between the data representing the tree and the metadata. If not identical, you won't be able to annotate the tree graphic.

# Required R packages:
- `distanceR` -> wrapper from the "cluster package"
- `ape` : [webpage](http://ape-package.ird.fr/), [Manual](https://cran.r-project.org/web/packages/ape/ape.pdf)
- `ggtree` : [vignettes](http://www.bioconductor.org/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html)
- `devtools`: [github](https://github.com/r-lib/devtools) we will use it to install `distanceR` from github
- Haukon's `distanceR` [Haukon's github]

`
to install packages in R type:
```{R}
install.packages("package_name", dependencies = TRUE)
```

If you want to use the functions availables from packages hosted on GitHub like Haukon's package: [Haukon's github]
you need to install as such:

```{R}
library(devtools) # we need to load devtools package into the enviromnment before we can use it
install_github("hkaspersen/distanceR")
```
- [ ] This one I will have to test,

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
