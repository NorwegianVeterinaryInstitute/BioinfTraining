
# Visualization networks with Cytoscape

## Install [Cytoscape](https://cytoscape.org/)

#### Ubuntu
- you will need java 8 installed (or higher installed)
- you will need to download cytoscape and run the install file (follow instructions from Cytoscape webpage)

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

#### windows: follow instructions of cytoscape webpage

NB: for R enthusiasts there is also a way to control cytoscape from R with `rcy3` packages.
Or you can also with with specific R "network" packages

## Import network
> if you do not want all the path in sample names you can do a little cleaning (R or in python)

Optional: an exemple little [python script](https://github.com/evezeyl/R_poppunk/blob/master/accessory_scripts_py_R/cytoscape_labels.py) add a column without path and extension accolated to sample id (ask if you need improvements) so you can use this as labels:
```
# usage: python 3 environment
python path/cytoscape_labels.py <input.csv> <path_pattern_to_remove> <extension_to_remove>

# example:
python ~/scripts/cytoscape_labels.py all_cytoscape.csv sequences/ .fas
```
After cleaning you can import your file:
- import network from file > select `*.graphml`
- import table from file > select `*.cytoscape_(modified).csv`
  - select `id` (key : column id that is identical in the network and data you want to associate with the network): -> this binds the graph to the isolates clusters numbers and names

- go to styles:
  - you can change the `shape`
  - `fill color` : column: combined_cluster, mapping type: discret -> then on right size can choose the color for each cluster

- export network to image file

# Tree vizualisation

Please look at our [R_trees course](./R_trees.md) for importing, visualizing and annotating trees with R
