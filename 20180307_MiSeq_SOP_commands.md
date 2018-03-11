Introduce the commands for doing beta-Diversity

plotting a heatmap with the abundances of the first 50 OTUs in all samples
```
heatmap.bin(shared=final.0.03.subsample.shared, scale=log2, numotu=50)
```

Next we compare using similarity of membership and structure the samples

```
dist.shared(shared=final.shared, calc=jclass-braycurtis-thetayc,
 subsample=2401)
```
then vizualize the distance matrices using `heatmap.sim`

```
heatmap.sim(phylip=final.jclass.0.03.lt.ave.dist)

heatmap.sim(phylip=final.braycurtis.0.03.lt.ave.dist)

heatmap.sim(phylip=final.thetayc.0.03.lt.ave.dist)
```

Next we create venn diagrams of the first four samples.

```
venn(shared=final.0.03.subsample.shared, groups=F3D0-F3D1-F3D2-F3D3)
```

Creating a dendrogram showing the relationships between the samples based
on their similarity.
```
tree.shared(phylip=final.thetayc.0.03.lt.ave.dist)
```
View the `*.tre` file with the program:
 `treeviewX`(available on your biolinux machine)

Next we test using the "parsimony" command if there is significant clustering

```
parsimony(tree=final.thetayc.0.03.lt.ave.tre, group=mouse.time.design,
   groups=all)
```

Next we can vizualize the sample clustering using Principal Coordinate Analysis

```
pcoa(phylip=final.thetayc.0.03.lt.ave.dist)
```
This produces output on screen.
it shows the R-square or correlation between the original distance matrix and the dimensions in pcoa space. so the R-square in 2D space is 0.88 and in 3D space almost 0.98.



Next up is non-metric multidimensional scaling

```
nmds(phylip=final.thetayc.0.03.lt.ave.dist)
```
and with three dimensions

```
nmds(phylip=final.thetayc.0.03.lt.ave.dist, mindim=3, maxdim=3)
```
Here we see that the stress level = 0.04.
stress should preferably be below 0.20 to trust the results from the ordination. The lower the better.

Then we do some testing to see if these ordinations are different

First we try analysis of molecular variance

```
amova(phylip=final.thetayc.0.03.lt.ave.dist, design=mouse.time.design)
```

This shows that the separation is significant.
Next we test the variance with homova

```
homova(phylip=final.thetayc.0.03.lt.ave.dist, design=mouse.time.design)
```
This show that there is also a significant difference in the variation between the samples of early and late. That we also saw with the early beta.diversity estimators.


Testing which OTUs are responsible for the ordination patterns

Here we use the command: metastats
```
metastats(shared=final.0.03.subsample.shared, design=mouse.time.design)
```

check the output file to identify which otu's have a q-value below 0.05.
