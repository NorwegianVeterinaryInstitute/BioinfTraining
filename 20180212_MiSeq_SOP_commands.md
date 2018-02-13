# Mothur commands used 20180205
This section will cover

1. Removing chimeric sequences
2. Removing unwanted taxa from the dataset
3. Calculate the error rate using the Mock community
4. Determine OTUs from the tutorial data at thr 97% sequence similarity level.

### Searching for chimeric sequences with chimera.uchime
Note that we use uchime here, but the current mothur standard is Vsearch and is recommended for chimera checking. We use uchime since the biolinux version does not have Vsearch.

`chimera.uchime(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table)`

### Removing chimeric sequences, the INCORRECT way
`remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos)`

`summary.seqs(fasta=current, count=current)`

### Removing chimeric sequences, the CORRECT way
`remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table)`

Why is this the correct way?

`summary.seqs(fasta=current, count=current)`

### Classifying seqs to detect undesired lineages
`classify.seqs(fasta=current, count=current, taxonomy=trainset9_032012.pds.tax, reference=trainset9_032012.pds.fasta, cutoff=80)`

### Removing undesired lineages
`remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloropast-Mitochondria-unknown-Archaea-Eukaryota)`

### Summary.seq file to create an overview of the clean dataset
`summary.tax(taxonomy=current, count=current)`

### Summary.tax to create an overview of the taxonomic classifications present in the clean dataset
`summary.tax(taxonomy=current, count=current)`

**At this point we have a clean dataset which can be used for generating OTUs and doing diversity studies. No matter what toolor  pipeline you have used al of the above filtering/removal /screening steps have to be performed. Of course this is debatable.**


##Assessing error rates
In this section we use a mock community consisting of 32 unaligned fasta sequences, to analyze the error rate that we have based on the clean-up that we have performed.

### extracting the mock community dataset from the complete dataset
`get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, groups=Mock)`

### Calculating the seq error rate using the mock community fasta sequences

`seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table, reference=HMP_MOCK.v35.fasta, aligned=F)`


### Calculating distances in the mock community sequences
`dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, cutoff=0.03)`

### Clustering sequences based on distances to generate OTUs
`cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table, cutoff=0.03)`

### Make shared command
`make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table, label=0.03)`

### Rarefaction command
`rarefaction.single(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared)`
This generates a file which you can us to check how the diversity is distributed in the dataset.



## Preparing for clustering and OTU picking

### Removing the mock communtiy from the dataset
`remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, groups=Mock)`

### Build a distance matrix with the maximum distance set
`dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, cutoff=0.20)`

Why is the cutoff used here?

### Cluster sequences using average linkage clustering
`cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table, method=average, cutoff=0.20)`
Note: in the newest MiSeq SOP this method is replaced by the `optiClust` algorithm.

### Create an otu table, as well as relative abundance files for each sample
`make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table, label=0.03)`

### Classifying the OTU sequences see in what taxon they fall.
`classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy, label=0.03)`
