## Installing Dtool 
`pip install dtool`

## In Saga
`conda activate dtool`

## Creating data set 
`dtool create Ecoli-Test`

The above command creates a folder called "Ecoli-Test" and adds "data" folder and README.yml file into "Ecoli-Test"
"Ecoli-Test" is called a dataset. 

## Adding data to dat set
`Dtool add item Test.fastq Ecoli-Test`

## Adding metadata to dataset  
`dtool readme interactive Ecoli-Test`

The above command will ask for the information below

`description [Dataset description]:`
`project [Project name]:`
`confidential [False]:`
`personally_identifiable_information [False]:`
`name [Your Name]:`
`email [you@example.com]:`
`username [jeevka]:`

Updated readme
To edit the readme using your default editor:
dtool readme edit file://Falcon-Pro.local/Users/jeevka/Dtool/Test2

# Checking the details of the dataset
`dtool ls Ecoli-Test`
