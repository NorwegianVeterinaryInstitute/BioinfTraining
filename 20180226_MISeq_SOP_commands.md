# Mothur commands used 20180226

## Preparing for clustering and OTU picking

### Removing the mock communtiy from the dataset
```
remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.count_table,
 fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, 
 taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy,
 groups=Mock)
```

### Build a distance matrix with the maximum distance set
```
dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta,
 cutoff=0.20, processors=2)
```

Why is the cutoff used here? and how does mothur use that cutoff to build the distance matrix

### Cluster sequences using average linkage clustering
```
cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist,
  count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table,
  method=average, cutoff=0.20)
```

Note: in the newest MiSeq SOP this method is replaced by the `optiClust` algorithm.

### Create an otu table, as well as relative abundance files for each sample
```
make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.list,
 count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table,
 label=0.03)
```

### Classifying the OTU sequences see in what taxon they fall.
```
classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.list,
 count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table,
 taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy,
 label=0.03)
```

### renaming output files so we can use them easily in the diversity analysis

```
system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist final.dist)

system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta final.fasta)

system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an_unique_list final.list)

system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an_unique_list.0.03.cons.taxonomy final.0.03.taxonomy)

system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an_unique_list.0.03.cons.tax.summary final.0.03.tax.summary)

system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an_unique_list.shared final.shared

system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table final.count_table)

```
**Now we have the final datasets, that can be used for diversity analysis using otu**


##Preparing for diversity analysis

### identifying the sample with the least sequences
```
count.groups(shared=final.shared)
```

### subsampling all datasets to downsample to the smallest sample size
It is important that you check if your samples match the **_size_**(2401) indicated here, if not that use the smallest sample size in your analysis.
```
subsample(shared=final.shared, size=2401)
```

### calculating data to create rarefaction curves
```
rarefaction.single(shared=final.0.03.subsample.shared,
 calc=sobs, freq=100)
```

### calculating alpha-diversity estimators using subsampling
```
summary.single(shared=final.shared,
 calc=nseqs-coverage-sobs-invsimpson,
 subsample=2401)
```


