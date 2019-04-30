
# Visualization networks with Cytoscape


## Install [Cytoscape](https://cytoscape.org/)

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

> for windows: follow instructions of cytoscape web page

NB: for R enthusiasts there is also a way to control cytoscape from R with `rcy3` packages.
Or you can also with with specific R "network" packages

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

# Tree vizualisation: ressources

Please refer to our [R_trees course](./R_trees.md) for importing, visualizing and annotating trees with R
