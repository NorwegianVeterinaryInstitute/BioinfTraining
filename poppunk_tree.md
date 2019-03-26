# What is [PopPUNK](https://genome.cshlp.org/content/29/2/304) : Population Partitioning Using Nucleotide K-mers

[Manual](https://poppunk.readthedocs.io/)
[PopPUNK Documentation](https://poppunk.readthedocs.io/en/latest/index.html)
[Presentation](https://docs.google.com/presentation/d/1StmmM02lSpFPdevQT3iDB3BAKRoMC8Q4NtJ4nzu7MdY/edit?usp=sharing)

- whole genome (core + assessory) population analysis/clustering
- distinction between isolates : uses k-mer of different length: `mash` to find core and accessory distances between isolates (**pairwise**)
- the distribution of those distances is used to discriminate clusters (defined as strains) of closely related isolates (similarity:both core and accessory)
- clustering of newly added isolates: EXTENDABLE4 - without the need of reanalizing all samples
- maintenance free and auto-reduce database
> aimed for consistent naming clusters between studies, outbreak detection in minutes


# How to use PopPUNK




# What poppunk does and questions: to prepare presentation


# Poppunk commands


poppunk --fit-model --distances <folder/folder.dists> --ref-db <folder_to_db/model_db> --output <output_folder> --full-db --K <nb_blubs>
