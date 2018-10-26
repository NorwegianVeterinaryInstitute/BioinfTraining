#How to download multiple genomes.

This page is based on the work of Kai Bling: [https://github.com/kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download). We recommend that you visit his github repository for seeing an extended list of examples on how to download microbial genomes from the NCBI web page. Here we describe how you can use it on the abel cluster, when using the Veterinary Institute installation of the NCBI-genome-download scripts.

In order to be able to download bacterial genomes to your directory of choice, you first need to start the conda environment that contains the scripts.

#####Activating the environment
	
	conda activate ncbidown
	
that will create a prompt on your screen that looks something like this:

	(nbcidown) thhaverk@login-0-0

#####Downloading all *Staphylococcus epidermidis* genomes from refseq
	ncbi-genome-download --genus "Staphylococcus epidermidis" bacteria
	
This would download for all 482 genomes the genbank flatfiles (extension: *.gbff.gz). However, many of these genomes are not finished genomes. If you want to only have genomes that are complete and only in FASTA format (e.g. without any annotations), than use the following command:
	
	ncbi-genome-download --format fasta --assembly-level complete --genus "Staphylococcus epidermidis" bacteria
	
##### downloading multiple species
It could be interesting to compare the pangenomes of *Staphylococcus aureus* with *Staphylococcus epidermidis*. Thus we need to download the genomes of both species. We can use the taxon ids for both species, which are found [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1280&lvl=3&lin=f&keep=1&srchmode=1&unlock): 

	ncbi-genome-download --taxid 1280,1282  --assembly-level complete bacteria
	
For more examples please visit the page of the program [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
		

 