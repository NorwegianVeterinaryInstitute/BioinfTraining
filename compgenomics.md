# Comparative genomics

Teacher: Karin Lagesen


## Format for this course

The course takes place on Mondays, at 9am, in Resepten. Remember to bring your
laptops. We will be running analyses on abel, and using tools within the 
virtual machine for evaluation of the results.

Each week we will do the following:
  * Introduction to today's theme
  * Practical exercise on a small dataset
  * Homework, where you will perform the same types of exercise on a different 
  given dataset
  
Most weeks, two students will be asked to review some more material regarding
how something works, and then prepare a short (10-15 minutes) presentation on
the subject.


## Practicals

For the first part of this course, we will follow the assembly exercises
that are part of the INF-BIOx121 course at the University of Oslo.

[Link to the relevant pages](https://github.com/karinlag/INF-BIOx121/tree/2017/Assembly/practicals) 


### 2018-02-05 ###

#### Preparatory work

Make sure that access to abel works, and that the jupyter notebook (anaconda)
is installed on the virtual machine.

#### Today's practical

Program for today:
  * Short lecture on how assembly works
  * Experiment with de Bruijn graphs, using the DeBruijnGraph.ipynb notebook
  * Go through most of the Velvet exercise from the INF-BIO site
  
The dataset that we will use for the in-class exercises can be found here:

`/work/projects/nn9305k/bioinf_course/compgenomics/rawdata/inclass`

Note, for the INF-BIO course, the students were working on a non-abel computer.
To make these exercises work, we have to do some things to make things work.

1. Log onto abel
2. Check which login node you are on (`hostname`)
3. Start the `screen` program
4. Go to your project home directory
5. Create a course directory, call it `2018_compgenome`, and go into it
6. Ask for a qlogin session:  
`qlogin --account=nn9305k --time=03:00:00 --ntasks=2 --mem-per-cpu=4G`
7. To run the velvet program, we have to use the `module` system to load it:  
`module load velvet`

We can now work through the exercises more or less as described. Note,
we only have paired end data, and will only be working with that data.

We will be using a google spreadsheet to record our assembly data:

[Velvet k-mer test](https://docs.google.com/spreadsheets/d/1mvIV0jenKBWGxIVyHTMe2Stb2OkNuLRv0iqpliC3URY/edit?usp=sharing)


When it comes to downloading and running notebooks:
1. On github, click on the filename so that the notebook is displayed.
2. Click on `raw` in the upper corner
3. Copy the URL in the address field
4. Go to a terminal on your local computer, and type in `wget URL`
5. To run the notebook: `jupyter notebook filename.ipynb`


#### Homework
 
There are six genomes in the directory 

`/work/projects/nn9305k/bioinf_course/compgenomics/rawdata/homework`

You will work with a partner. Each team will pick three k-mer sizes and use 
velvet to assemble the six genomes. Each person in the team is responsible
for three genomes each.

[Record your results here](https://docs.google.com/spreadsheets/d/124Eb6IQ44coSKMH0kRLU18AJ5FZ7-ijwsxqf1NsC9Ys/edit?usp=sharing)




