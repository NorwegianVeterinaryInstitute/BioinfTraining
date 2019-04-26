start with data prep
send a link

# What is R short introduction
R and R studio
console

# Data:
cgMLST data
pairwise distance matrix : Metrics - "gowermetric"

Labeling - distance labels and metadata label for isolates has to be **perfectly identical** (because it allows to create a link between those 2 datasets)

# Required R packages:
`distanceR` -> wrapper from the "cluster package"
`ape`
`ggtree`
`haukons package` !!!
`devtools`
`
to install packages in R type
install.packages("package_name", dependencies = TRUE)

if you want to use haukons skrip: as it is a git hub repository : <https://github.com/hkaspersen/distanceR.git>
you need to install as such:




# Create a tree - and simple plot (just labels)
- E.coli data from haukon download from Abel : they are found here ``



export - Eve ...

# Importing an existing tree -> Eve
you can important trees directly from other package
- from specific programs ->
- https://yulab-smu.github.io/treedata-book/index.html can look here for a short introduction to ggtree


# Add decorations - link the metadata
collors annotations heatmap

each step has to build on the previous one
but you need to get something on screen, very fast
```{R code}

```

# Going further
### other tree related packages but not used R

`phylo`
treeio
phangorm

### links
gheatmap -> eve check
> if use gheatmap -> note row.names have to be exacly the same as isolate names in tree/metadata (in col1) - in heatmap the order of rows does not need to be ordered
