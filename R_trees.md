# Introduction to tree generation and visualization in R

# R and RStudio
R is a free software environment for statistical computing and graphics, and can be downloaded [here](https://cran.uib.no/). RStudio is an integrated development environment for R, which in short makes R easier to use. Rstudio can be downloaded [here](https://www.rstudio.com/products/rstudio/download/).

For this session, you will need both R and RStudio.

R is an object orientated programming language. You can assign values to an object with the name of your choosing, for example `my_color <- "red"`. When doing so, the object is stored in the computers memory. The object `my_color` will then hold the value `"red"`, and can be called by typing: `my_color`.

# Projects
In RStudio, one can create a project to store all scripts and notes from one project in one folder. After starting RStudio, go to `file` and `New project`. In the window, click `New directory`, `New project`, and type in the name of the project folder, for example `Ggtree tutorial`. Choose which directory to create the project in, and click `create project`.

Now, RStudio will restart, and the working directory of R will be in the folder you specified above. The advantage of dividing your work into projects is that all your scripts and files can be easily accessed in the project, or subfolders of the project, making it easier to have control of your code. The project folders may also be under version control with git.

Now we are ready to install some packages!

> NOTE: In the code boxes below you will see comments starting with a `#`. This is not code, but helpful comments that explain the arguments and functions found in the box.

# Packages
In R, the fundamental unit of shareable code is a package. There is a myriad of packages available for download, and the ones published on the Comprehensive R Archive Network [CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html) are trustworthy and of high quality. Packages can be created by anyone, and can also be hosted on GitHub. Some packages, specialized for biological analyses, are hosted by [Bioconductor](https://www.bioconductor.org/), and have their own installation method. Note that all installed libraries are available for activation regardless of which project folder you are in.

To install packages from CRAN, use the following command:
```{R}
install.packages("package_name")
```

For this session, you will need the following packages
- `ape` : [website](http://ape-package.ird.fr/), [manual](https://cran.r-project.org/web/packages/ape/ape.pdf)
- `cluster` : [documentation](https://cran.r-project.org/web/packages/cluster/cluster.pdf)
- `tibble` : [vignettes](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html) *(vignettes are nearly equivalent to manual pages.)*
- `dplyr` : [website](https://dplyr.tidyverse.org/)
- `ggtree` : [vignettes](http://www.bioconductor.org/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html). We have to install this package from GitHub (see below)
- `distanceR` Also need to be installed from GitHub: [Haukons' GitHub](https://github.com/hkaspersen/distanceR)

To be able to install packages from GitHub, you will need the following package:
- `devtools`: [manual hosted on GitHub](https://github.com/r-lib/devtools)

To install packages from github, you need to use the following functions:

```{R}
# Loading the packages we will need in the R environment with the "library" function
# NB: loading a package makes the functions defined in the package available in R
# This only has to be done once per session

library(devtools)
install_github("GuangchuangYu/ggtree")
install_github("hkaspersen/distanceR")
```

# Training data
The data files used in this session can be downloaded below (left click "download as..." on the link). Copy or download directly into the project folder you created earlier.

[cgMLST data](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/cgMLST.tsv)

[Tree metadata](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/tree_metadata.txt)

[Tree heatmap data](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/hk_ef_R_trees/training_files/tree_heatmap_data.txt)

NB: if you are working on a Linux machine, you can download those files from the command line with: `wget <url_file>`

# Calculating distances from cgMLST data
Based on the allele data in the cgMLST file, one can calculate the percentage of similar "labels" in a pairwise fashion. The labels in the cgMLST data represent the loci variant, based on the nucleotide sequence (same as regular 7 gene MLST).

In other words, the first sample in the data is compared to the second, and the percentage of similar labels is calculated (excluding the loci if one label is missing). The first sample is then compared to the next sample, until all sample pairs are compared.

The result of this calculation is represented in a symmetrical dissimilarity matrix, which is N x N big. Diagonals of the matrix represent each isolate compared to iself, and is therefore set as 0, as there are no dissimilarities in the comparison. The rest of the comparisons result in values between 1 and 0, where values closer to zero represent a more similar comparison.

This data can then be clustered hierachically (here by the UPGMA method), and a tree object can be created from the resulting clusters.

To do this in R, one has to import and clean the cgMLST data first:
```{R}
# loading the packages we will need in the R environment

library(tibble)
library(dplyr)

cgMLST_data <- read.table("cgMLST.tsv",             # The cgMLST data file name
                          sep = "\t",               # The separator in the file, here tab
                          colClasses = "factor",    # Set all columns in the data to factor
                          header = TRUE) %>%        # Does the file contain column names?
  na_if("0") %>%                                    # change all "0" to NA
  column_to_rownames("FILE")                        # Change rownames to the values in the "FILE" column

```
Here you see the `%>%` operator. This operator is a pipe, and sends the output from one function to the next, similar to `|` in bash.
Note that a specific data structure is needed to calculate distances correctly. If an allele is missing, it should be represented by `NA`, which is why we change all `0` to `NA`. Also note that the row names need to have the sample ID's, which can be found in the "FILE" column. A factor is a categorical variable.

To calculate distances from the data and create a tree:

```{R}
library(cluster)
library(ape)

# calculate distances
distances <- daisy(cgMLST_data, metric = "gower")

# cluster the samples from the dissimilarity matrix
clust_dist <- hclust(distances, method = "average")

# create tree object
tree <- as.phylo(clust_dist)

# to do all three calculations at once you can nest your commands within each other:
tree <- as.phylo(hclust(daisy(cgMLST_data, metric = "gower"), method = "average"))
```
Check out the function information with `?daisy` and `?hclust` to see details on the metric and method arguments

To then visualize this tree, we can either use the `plot()` function from base R, or use `ggtree()` function from ggtree package.

*NB: base R is a set of packages that are loaded automatically when you start R*
> What are the advantages of using ggtree?
> - using base plot() function allows us to have a fast view at our tree object
> - ggtree has advanced anotation functions that allow us to add layers, and follows a specific syntax that allows to do advanced graphics "automatically and easy". This special syntax is called ["grammar of graphics"](https://en.wikipedia.org/wiki/Leland_Wilkinson) which is the gg from the ggtree or ggplot2 package.

> What are layers?:
  > -  Imagine that I made a basic graphic on a white A4 paper. We can add labels, annotations, colors, title with layers: think of a layer as a transparent sheet where you only draw one type of information (ex: leaves labels) and you put this transparent sheet on top of your basic tree. You can now see your tree annotated with the leaves labels.

```{R}
library(ggtree)

plot(tree)
ggtree(tree)
```
In base R, the tree will have sample name labels on the tips. In ggtree, nothing was specified, so only the tree was plotted. Ggtree have multiple layouts to choose from, i.e. "circular", "rectangular", "fan", "equal_angle" etc. See the vignette on ggtree above for more information. To add more information to the tree with ggtree, we can use additional layers of information, like this:

```{R}
ggtree(tree,
       layout = "rectangular") +
     geom_tiplab()
```
Now, the isolate names have been added to the tree tip points (leaves labels). However, if you want to combine metadata such as sequence types, animal species etc. to the tree and visualize it, one first has to connect the tree and the metadata:

# Linking metadata to your tree data allows you to add specific annotations

```{R}
# Import metadata
metadata <- read.table("tree_metadata.txt", sep = "\t", header = TRUE)

# Plot tree
# The '%<+%' operator links the data to the tree
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
       layout = "rectangular") %<+% metadata +    # links your metadata to the graph
     geom_tiplab(aes(label = ST)) +               # add a labels layer (for the leaves), labels are taken from ST column
     geom_tippoint(aes(color = species))          # add a points layer, each species represented by its own color
```

`ggtree()` will here plot the animal species in the metadata as colored nodes on the tree. Since no specific palette is specified, it will use default ggplot2 coloring. If you want to specify colors, a specific palette need to be created:

```{R}
library(ggplot2)

# create palette
palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43")

# Plot tree and assign to object
my_tree <- ggtree(tree,
       layout = "rectangular") %<+% metadata +
     geom_tiplab(aes(label = ST)) +
     geom_tippoint(aes(color = species)) +
     scale_color_manual(values = palette)   # here we specify the palette
     
# To look at the tree, type the name of the object:
my_tree
```
This will then plot the tree with the specified colors. Further adjustments to the look of the tree can be made by adding arguments like 'size', 'alpha' (transparency), 'offset' (distance of labels from the tree) and others, on each layer. Make sure to put these arguments outside the `aes()` function (unless you want it to represent something in your metadata).

Adding a legend to the tree if it does not show-up automatically: look at `?theme` for all options.
Ex: add the following line to commands specifying `my_tree`: `+ theme(legend.position = "bottom")`

> *NB: If you did not save your script, you can run [R_tree_course_script](https://raw.githubusercontent.com/NorwegianVeterinaryInstitute/BioinfTraining/master/training_files/R_tree_course.R)
to catch up fast up to this point*

# Adding heatmaps
We have now created an annotated tree connected to our metadata. However, one can also add a heatmap to the outside of the tree, representing for example presence/absence of genes. To do this, we will use the function `gheatmap` from ggtree. NOTE: The heatmap data need to have the rownames as the sample names, unlike the metadata file, which had the first column as the sample names. We fix this when importing the table into R below.

> NB: a heatmap is a plot where colors are associated to values in your data (ex: white will be high value, red will be low value). If you plot a heatmap of a matrix, you will see low and high values as a gradient of colors.  

```{R}
# Import heatmap data
heatmap_data <- read.table("tree_heatmap_data.txt",
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE) %>%
  mutate_at(vars(-id),
            funs(as.character)) %>% # change columns to character (-id from dataframe representation, generates a named list)
  column_to_rownames("id") # set row names as the values in the column "id" (passes the named list to rownames)

gheatmap(my_tree,               # your tree
         heatmap_data,          # the heatmap data
         offset = 0.07,         # distance of heatmap from tree tips
         width = 0.5,           # the width of the heatmap
         font.size = 3)         # the font size of the column headers
         
```
Take a look at the tree, and see that the column headers looks like they are all jumbled up. We can fix this by changing the direction of the column headers and placing them on top of the plot:

```{R}
gheatmap(my_tree,               
         heatmap_data,         
         offset = 0.07,
         width = 0.5,
         font.size = 3,
         colnames_position = "top",
         colnames_angle = 90)
         
```
Similar to the tree annotation above, one can use a specific color palette of your choosing to represent the data.

NB: type `colors()` to see an example of colors you can use 

```{R}
palette2 <- c("0" = "grey95",
              "1" = "steelblue")

complete_tree <- gheatmap(my_tree,               
         heatmap_data,         
         offset = 0.07,
         width = 0.5,
         font.size = 3,
         colnames_position = "top",
         colnames_angle = 90) +
  scale_fill_manual(values = palette2)

complete_tree
```

Now you can save the tree with high resolution with the `ggsave` function.

```{R}
ggsave("complete_tree.tiff",  # filename of your choosing
       complete_tree,         # your tree
       device = "tiff",       # what kind of file should it be saved as?
       dpi = 300,             # resolution of the tree
       units = "cm",          # units to use in the height and width arguments
       height = 20,           # height of the image
       width = 20)            # width of the image
```
The image file can be found in your project folder.

Available file formats you can be accessed with `?devices`, ex: bmp, jpg, pdf ...

Finally, there is a lot of functionality in ggtree and other functions that isn't discussed here. use `?ggtree` and `?gheatmap` in R to see additional arguments and possibilities for different visualization of trees.

# The distanceR package
The distanceR package make it easier to calculate trees, annotate and add heatmaps without the use of all the functions above. The package is available on Håkon's github page listed above.
In short, the functions can be used as follows:

```{R}
library(distanceR)

# calculate tree
tree <- calc_tree("cgMLST.tsv",
                  metric = "gower",
                  method = "average")
                  
# note that in the function above, the file in the 
# folder is specified, not an object in R. That is 
# because the function imports the data in the correct
# format and then runs distance calculation.

# Annotate tree with metadata
annotated_tree <- annotate_tree(tree,
                                "tree_metadata.txt",
                                layout = "rectangular",
                                label_variable = "ST",      # The labels on the tips ("ST": column name in the metadata)
                                color_variable = "species")  # The colored nodes

# Similar to above, the file in the folder is specified as the metadata. 
# Additional settings can be added, such as point size, font size etc. 
# See `?annotate_tree` for details.

# Add heatmap to the tree
add_heatmap(annot_tree,
            "testdata/heatmap_data.txt",
            layout = "rectangular",
            colnames_position = "top",
            colnames_angle = 90,
            colnames_offset = 0.5,
            font_size = 2,
            heatmap_offset = 0.1)

```

# Going further
## Importing an existing tree
It is possible to import trees that were created which clustering/phylogenetic softwares.
The [`treeio` package](https://bioconductor.org/packages/release/bioc/html/treeio.html), that has been developped from parsing different trees format into R.
And have been designed to import both tree and metadata associated with the tree.

As written in [data Integration, Manipulation and Visualisation of phylogenetic trees] in chapter: [linking ananotation data to tree with tidytree](https://yulab-smu.github.io/treedata-book/chapter7.html)

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

## There are several other packages and other functions in the packages we presented 

For `[tidytree]` , `[treeio]` and `[ggtree]` -> many more functions are described in [data Integration, Manipulation and Visualisation of phylogenetic trees]

`[ggplot2]` has some additional functions that can be used to annotate trees (note that: ggtree package is based on ggplo2)

`ape` has also many more functions to manipulate trees. See: [Analysis of Phylogenetics and Evolution with R](http://ape-package.ird.fr/APER.html).

... and there are many more ...

[tidytree]:https://cran.r-project.org/web/packages/tidytree/vignettes/tidytree.html
[ggplot2]:https://www.rdocumentation.org/packages/ggplot2/versions/3.1.1 
[ggtree]:https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html
[treeio]:https://github.com/GuangchuangYu/treeio

## How to find R packages: 

[Bioconductor package portal](https://bioconductor.org/packages/release/BiocViews.html#___Software)

[CRAN Task view](https://cran.r-project.org/) R packages sorted by themes

There is no overview of all the R-packages hosted on GitHub, so we need some `googling` to find them.

[data Integration, Manipulation and Visualisation of phylogenetic trees]:https://yulab-smu.github.io/treedata-book/index.html
[Håkon's github]:https://github.com/hkaspersen/distanceR.git

## Other resources
A good introduction for learning how to use trees data in R: [data Integration, Manipulation and Visualisation of phylogenetic trees]. (We used it a lot to make this page.)

