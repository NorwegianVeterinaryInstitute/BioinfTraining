# Mothur commands used 20180205

### Combining forward and reverse reads
`make.contigs(file=stability.files, processors=2)`

`summary.seqs(fasta=stability.trim.contigs.fasta)`

### Remove sequence that are too long, or have ambigous bases
`screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=275)`

### Repeat of above, but using the summary output as well
`screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=275, summary=stability.trim.contigs.summary)`

### Obtaining what is in memory
`get.current()`

`summary.seqs(fasta=current)`

### Creating unique sequences by clustering identical sequences, creates a names file
`unique.seqs(fasta=stability.trim.contigs.good.fasta)`

### Counting sequences in the names file to create a count_table
`count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)`
Check that the output file has all the groups in various columns.

### Summary of the dataset sofar, note this uses the fasta, groups and names files as well
`summary.seqs(count=stability.trim.contigs.good.count_table)`

### Extracting the 16S rRNA V4 region from the silva.bacteria.fasta database
`pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=2)`

### Renaming the seqs.pcr output file
NOTE this command does not work in mothur on the biolinux machine. This needs to be done manually!!!

Nonetheless here is the command:
`rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)`

### Summary of the silva.v4.fasta file
`summary.seqs(fasta=silva.v4.fasta)`

### Aligning the contigs against the silva v4 dataset
`align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)`

### Summary of the aligned sequences
`summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)`

### Removing sequences that do not align correctly to the V4 region, or that have too many ambigous bases
this duplicated the last sequence in the count_table file.
the duplication was manually removed

`screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8)`

### summary after screening
`summary.seqs(fasta=current, count=current)` 

**or**

 `summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table, processors=2)`

### Filtering the alignment from empty columns and removing overhangs from the alignment
`filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)`

### Rerunning the unique.seqs command to adjust for the sequences being indentical
`unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)`

### Running precluster to account for sequencing and pcr errors
`pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)`

