# visualization networks with Cytoscape
> from PopPUNK output

- from ubuntu (native install)
- from R - with bioconductor

## For Ubuntu: Install [Cytoscape](https://cytoscape.org/)

> for ubuntu
- you will need java 8 installed (might work with higher):
- you will need to download cytoscape and run the install file

```
# install java 8
sudo apt update
sudo apt install openjdk-8-jdk

# go to the direcory where you saved Cytoscape
chmod u+x  Cytoscape_3_7_1_unix.sh
sudo sh Cytoscape_3_7_1_unix.sh

#if error:
sudo apt install libcanberra-gtk-module libcanberra-gtk3-module

# run cytoscape
Cyctoscape &

```

### Import network
> if you do not want all the path in sample names you can do a little cleaning

> little script to embelish: remove path and extension from sample id
usage: python 3 environment
`python path/cytoscape_labels.py <input.csv> <path_pattern_to_remove> <extension_to_remove>`
ex: `python /home/evezeyl/Documents/gits/Listeria_listadapt/Evfi_test/scripts/cytoscape_labels.py all_in_one3_cytoscape.csv sequences/ .fas`

- import network from file > select *.graphml
- import table from file > select * cytoscape_(modified).csv
  - select id (with the path as key : as it is this key that is used in the netwok) -> this will bind the graph to the isolates clusters numbers and names

- go to styles:
  - can change the shape
  - fill color : column: combined_cluster, mapping type: discret -> then on right size can choose color for each cluster

- export network to image file - to export the vizualization

## R cytoscape - in bioconductor
### install bioconductor on your laptop
[R studio in different conda environmnets](https://support.rstudio.com/hc/en-us/articles/200486138-Changing-R-versions-for-RStudio-desktop)
How to in [techstuff](./techstuff.md#One-Rstudio-for-all-condas)

```
conda create --name rcy3 bioconductor-rcy3
conda activate rcy3
conda install -c conda-forge r-biocmanager
-c bioconda bioconductor-biocinstaller
rstudio &

#in R studio: to view documentation package:
browseVignettes("RCy3")
```
NB: rcytoscape for cytoscape 2.8.1
NB: RCy3 for cytoscape 3.x
