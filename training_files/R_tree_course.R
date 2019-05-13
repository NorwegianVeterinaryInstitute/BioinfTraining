library(tibble)
library(dplyr)
library(cluster)
library(ape)
library(ggtree)
library(ggplot2)


cgMLST_data <- read.table("cgMLST.tsv", sep = "\t", colClasses = "factor", 
                          header = TRUE) %>%
		na_if("0") %>%
		column_to_rownames("FILE")
                    
distances <- daisy(cgMLST_data, metric = "gower")

clust_dist <- hclust(distances, method = "average")

tree <- as.phylo(clust_dist)

metadata <- read.table("tree_metadata.txt", sep = "\t", header = TRUE)

palette <- c("Broiler" = "#4575b4",
             "Pig" = "#74add1",
             "Red fox" = "#f46d43")

my_tree <- ggtree(tree, layout = "rectangular") %<+% metadata +
     geom_tiplab(aes(label = ST)) +
     geom_tippoint(aes(color = species)) +
     scale_color_manual(values = palette) +
     theme(legend.position = "bottom")	# if does not add automatically you can add legend like that. 

?theme

# heatmap data need to have the rownames as the sample names, 
heatmap_data <- read.table("tree_heatmap_data.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
		mutate_at(vars(-id), funs(as.character)) %>% # change columns to character  
  		column_to_rownames("id") # set row names as the values in the column "id" 

gheatmap(my_tree, heatmap_data,          
         offset = 0.07, width = 0.5, font.size = 3) 

gheatmap(my_tree, heatmap_data,          
         offset = 0.07, width = 0.5, font.size = 3, 
		colnames_position = "top",
         	colnames_angle = 90)

palette2 <- c("0" = "grey95",  
              "1" = "steelblue")

complete_tree <- gheatmap(my_tree, heatmap_data,         
         offset = 0.07, width = 0.5, font.size = 3,
         colnames_position = "top",
         colnames_angle = 90) +
  	 scale_fill_manual(values = palette2) # note here the legend is moved on the left

complete_tree

ggsave("complete_tree.tiff",  complete_tree, 
       device = "tiff", dpi = 300,        
       units = "cm",  height = 20, width = 20)    

?devices


