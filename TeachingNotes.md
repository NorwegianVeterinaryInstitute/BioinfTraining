# Teaching notes

## 2018-03-05

### Homework

* Transfer the files to your computer, and open them in firefox. Most 
  of the data is untrimmed, and 150 bp or shorter. Only the ESBL is longer.
    
  Re fastqc results: you have to know your data. This time, we went into 
  it blind. We didn't check what we had. We could have had very crappy
  data. As it was, most of it was fine. However, check what you have.
  When we assembled things, we didn't even know the read length, which
  complicated things.
    
    * 151: 14042624, DTU2014 - can use auto
    * 100: F27, F159, F168 - can use auto
    * variable: S19 - need to fix
    
* Open one of the files with `less`. As you see it says HISEQ here.
    If you open the ESBL file, we will see M- which indicates Miseq data
    
* Re SPAdes manual: we need to make sure that we're looking at the manual
  for the right version.
    
  How to figure out what version of spades we have? 
    
  module load spades
  spades.py -v 
  
* QUAST results 
    
    