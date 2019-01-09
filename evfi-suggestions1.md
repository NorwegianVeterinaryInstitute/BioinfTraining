
# course itself - videos introduction + Lectures.pdf

Q: Questions

## Small things... 

**Suggestions:**

> maybe add complete reference for bibliography (end file) only to paste - in order from course. checked all links (in pdfs. to see if still valid). 
> All cited articles (but one - not free - asked for it) -> found and sorted in folder for course: if you want I can put all in your course folder (avoid others to have to search for articles if interested) - (I just love to have a look at articles cited in courses, as they are actually chosen because of quality/reprensentativity). 

-About:  ..the building blocks dna ... at the end of first course. Not sure what to think about it: very simple and clear but THEN - > there is a sudden jump to more complicated concepts that have not been presented. Maybe too limited? (ex: ? Arvind Speaks also of libraries for **transposons ...**) ...
- so maybe short background refresh: (problem increase to more of a genetic course): or evt provide ref for good review papers if needed... 
      - mobile elements /"junk"
      - methylation : -> Thomas: presentation about methylation detection (packBio) - epigenetics
      - RNAs - kinds, and functions ...exome sequencing
      
- biais - method/platform specific: Maybe you say that later in the course (think would be to know what are potential problems, pitfalls, -> also important to understand proprely articles). Mention ok...Arvind speaks a bit about optical duplicates course 3 (new-machines with 2 color dyes) -> how will influence analyses?
      
**For each lecture:**            
- Lecture1.pdf: Karch et. al. EMBO Mol. Med. 20**12** (typo)

- Lecture2.pdf: in case you  want to actuallize your links
  - p 22 Link: <https://www.sequencing.uio.no/forms/guidelines-submission-form-(illumina).pdf> they moved it to https
  - p 23 Link: <http://www.micronautomata.com/bioinformatics> does not seem to respond - retry...
  - p 31: Illumina - video link (see bellow) 
  
- Lecture 2 Video (Thomas) :
   - library: define earlier - maybe set the terms to know p 33 .  earlier in the presentation
   - Illumina ca 50 min on video
     Q: when watching video: was wondering: 
     - how fliped over for paired end reads - when first sequencing pass finished `atomatic flip when DNA fraqument amplification finished? - Did not really get at first the **bridge amplification** (pair ends) - thought: easier with wideo 
     - > Suggestion: show video illumina: explain the bridge add ? <https://www.youtube.com/watch?v=fCd6B5HRaZ8> but I find somehow older video cleared - its also shorter (also explain reads indexes 1 and 2)  <https://www.youtube.com/watch?v=HMyCqWhwB8E>    
    
Q: MinION - can get high quality consesus continuous sequencing -> not sure how see how sequence many times same seq? (so each pore 1 DNA molecule, so many porea ...some same DNA molecule -> consensus -> from many pores-- ? do I understand right?

   - p57 (was just wondering what were the different ending parts - as not yet explained : ex. P5-index2-Rd1SP)

- Lecture3.pdf: p3 link changed <https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-04.pdf>
   - so indexes - short tagging sequences -> representing the 2 complementary directions (old school: Forward and Reverse) of the sequenced DNA fragment.
   
- Lecture4.pdf (Karin): p 43: changed link to <http://www.cbcb.umd.edu/research/assembly_primer>
Q: - BLAST: K-mer World -> optimal choice (why not chose 36/17) - comparing all genome or part genome ? do we have idea about query seq?
  - Outch! *burrow wheeler transform (at the end of day) - comming back to that later on*
  
  - Ok- some clarrification assembly with choice length word....Q: From experience; lets say WGS - how genome structure (several chromosomes pairs - feks also with familly genes , or plasmid DNA -> bacteries ...affect assembly, think Bacteria, if HGT most likely plasmid ...) 

Q: **5 days course** - 5th day lecture - no video? - ? going throught the ressources?

### additional questions - to be sure I understand
- Q? Curious - whole genome amplification ? works ??? tested before making libraries (we did try, some ok, some not) - applicable - sure biais some taxa --- experience? 
- Q: depth and coverage:
   - so if denovo genome -> have no idea of coverage -> estimate -> expectated closely related?  => yes need example reference talking about ... works ok for bacteria - (what about plants - polyploidie)  
   - depth -> homomgeneity ...-> discriminate different OTU (use to help or unreliable?) 
   
- Q: good enough overlapp between reads (Assembly) -> Q: chose length word ~17+ + some mismatches - way to determine optimal?

- Q: 454: do not use (as not presented) - need to read a bit more. Maybe worth short presentation or link - to be able to uderstand if specificity when reading papers. 
--------------------
# ? ideas if not already suggested
- Introduction to Genomic Data Science - <https://courses.edx.org/courses/course-v1:UCSanDiegoX+CSE181.1x+1T2017/course/> - <https://stepik.org/> can help train programming with python: learn about sliding window, some text algorithms, charging machines ...
 
#  REFERENSE LIST: READY TO PASTE

Only one in bold = not found -> asked (research gate)

#### Lecture 1

Did not include the historical papers

Lex Nederbragt blog: https://flxlexblog.wordpress.com/2016/07/08/developments-in-high-throughput-sequencing-july-2016-edition/

Stratton, Michael R., Peter J. Campbell, and P. Andrew Futreal. “The Cancer Genome.” Nature 458, no. 7239 (April 2009): 719–24. https://doi.org/10.1038/nature07943.

Nasir, Arshan, Kyung Mo Kim, and Gustavo Caetano-Anolles. “Giant Viruses Coexisted with the Cellular Ancestors and Represent a Distinct Supergroup along with Superkingdoms Archaea, Bacteria and Eukarya.” BMC Evolutionary Biology 12, no. 1 (August 24, 2012): 156. https://doi.org/10.1186/1471-2148-12-156.

Kujiraoka, Manabu, Makoto Kuroda, Koji Asai, Tsuyoshi Sekizuka, Kengo Kato, Manabu Watanabe, Hiroshi Matsukiyo, et al. “Comprehensive Diagnosis of Bacterial Infection Associated with Acute Cholecystitis Using Metagenomic Approach.” Frontiers in Microbiology 8 (April 20, 2017). https://doi.org/10.3389/fmicb.2017.00685.

Tyson, Gene W., Jarrod Chapman, Philip Hugenholtz, Eric E. Allen, Rachna J. Ram, Paul M. Richardson, Victor V. Solovyev, Edward M. Rubin, Daniel S. Rokhsar, and Jillian F. Banfield. “Community Structure and Metabolism through Reconstruction of Microbial Genomes from the Environment.” Nature 428, no. 6978 (March 2004): 37–43. https://doi.org/10.1038/nature02340.

Karch, Helge, Erick Denamur, Ulrich Dobrindt, B. Brett Finlay, Regine Hengge, Ludgers Johannes, Eliora Z. Ron, Tone Tønjum, Philippe J. Sansonetti, and Miguel Vicente. “The Enemy within Us: Lessons from the 2011 European Escherichia Coli O104:H4 Outbreak.” EMBO Molecular Medicine 4, no. 9 (September 4, 2012): 841–48. https://doi.org/10.1002/emmm.201201662.

Rohde, Holger, Junjie Qin, Yujun Cui, Dongfang Li, Nicholas J. Loman, Moritz Hentschke, Wentong Chen, et al. “Open-Source Genomic Analysis of Shiga-Toxin–Producing E. Coli O104:H4.” New England Journal of Medicine 365, no. 8 (August 25, 2011): 718–24. https://doi.org/10.1056/NEJMoa1107643.

Hendriksen, Rene S., Lance B. Price, James M. Schupp, John D. Gillece, Rolf S. Kaas, David M. Engelthaler, Valeria Bortolaia, et al. “Population Genetics of Vibrio Cholerae from Nepal in 2010: Evidence on the Origin of the Haitian Outbreak.” MBio 2, no. 4 (September 1, 2011): e00157-11. https://doi.org/10.1128/mBio.00157-11.

**Falush, Daniel. “Bacterial Genomics: Microbial GWAS Coming of Age.” Nature Microbiology 1, no. 5 (May 2016): 16059. https://doi.org/10.1038/nmicrobiol.2016.59.**

Earle, Sarah G., Chieh-Hsi Wu, Jane Charlesworth, Nicole Stoesser, N. Claire Gordon, Timothy M. Walker, Chris C. A. Spencer, et al. “Identifying Lineage Effects When Controlling for Population Structure Improves Power in Bacterial Association Studies.” Nature Microbiology 1 (April 4, 2016): 16041. https://doi.org/10.1038/nmicrobiol.2016.41. (with supplement)

Valenzuela-Muñoz, Valentina, Jacqueline Chavez-Mardones, and Cristian Gallardo-Escárate. “RNA-Seq Analysis Evidences Multiple Gene Responses in Caligus Rogercresseyi Exposed to the Anti-Salmon Lice Drug Azamethiphos.” Aquaculture 446 (September 1, 2015): 156–66. https://doi.org/10.1016/j.aquaculture.2015.05.011.

#### Lecture 2
See Lecture 1:

Stratton, Michael R., Peter J. Campbell, and P. Andrew Futreal. “The Cancer Genome.” Nature 458, no. 7239 (April 2009): 719–24. https://doi.org/10.1038/nature07943.

Guidelines:
“Norwegian Sequencing Centre. Guidelines for Completion of Illumina Sample Submission Form.,” n.d. https://www.sequencing.uio.no/forms/samplesubmissionforms.html. https://www.sequencing.uio.no.

Articles:
Metzker, Michael L. “Sequencing in Real Time.” Nature Biotechnology 27, no. 2 (February 2009): 150–51. https://doi.org/10.1038/nbt0209-150.

Rhoads, Anthony, and Kin Fai Au. “PacBio Sequencing and Its Applications.” Genomics, Proteomics & Bioinformatics 13, no. 5 (October 2015): 278–89. https://doi.org/10.1016/j.gpb.2015.08.002.

Muinck, Eric J. de, Pål Trosvik, Gregor D. Gilfillan, Johannes R. Hov, and Arvind Y. M. Sundaram. “A Novel Ultra High-Throughput 16S RRNA Gene Amplicon Sequencing Library Preparation Method for the Illumina HiSeq Platform.” Microbiome 5, no. 1 (July 6, 2017): 68. https://doi.org/10.1186/s40168-017-0279-1.

Macosko, Evan Z., Anindita Basu, Rahul Satija, James Nemesh, Karthik Shekhar, Melissa Goldman, Itay Tirosh, et al. “Highly Parallel Genome-Wide Expression Profiling of Individual Cells Using Nanoliter Droplets.” Cell 161, no. 5 (May 21, 2015): 1202–14. https://doi.org/10.1016/j.cell.2015.05.002.

#### Lecture 3
Schurch, Nicholas J., Pietá Schofield, Marek Gierliński, Christian Cole, Alexander Sherstnev, Vijender Singh, Nicola Wrobel, et al. “How Many Biological Replicates Are Needed in an RNA-Seq Experiment and Which Differential Expression Tool Should You Use?” RNA (New York, N.Y.) 22, no. 6 (2016): 839–51. https://doi.org/10.1261/rna.053959.115.

Liu, Yuwen, Jie Zhou, and Kevin P. White. “RNA-Seq Differential Expression Studies: More Sequence or More Replication?” Bioinformatics (Oxford, England) 30, no. 3 (February 1, 2014): 301–4. https://doi.org/10.1093/bioinformatics/btt688.

Busby, Michele A., Chip Stewart, Chase A. Miller, Krzysztof R. Grzeda, and Gabor T. Marth. “Scotty: A Web Tool for Designing RNA-Seq Experiments to Measure Differential Gene Expression.” Bioinformatics 29, no. 5 (March 1, 2013): 656–57. https://doi.org/10.1093/bioinformatics/btt015.

Sims, David, Ian Sudbery, Nicholas E. Ilott, Andreas Heger, and Chris P. Ponting. “Sequencing Depth and Coverage: Key Considerations in Genomic Analyses.” Nature Reviews. Genetics 15, no. 2 (February 2014): 121–32. https://doi.org/10.1038/nrg3642.

Conesa, Ana, Pedro Madrigal, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, et al. “A Survey of Best Practices for RNA-Seq Data Analysis.” Genome Biology 17, no. 1 (January 26, 2016): 13. https://doi.org/10.1186/s13059-016-0881-8.



#### Lecture 4

La Trobe Institute for Molecular Science, Melbourne, Australia, Thomas Shafee, Rohan Lowe, and La Trobe Institute for Molecular Science, Melbourne, Australia. “Eukaryotic and Prokaryotic Gene Structure.” WikiJournal of Medicine 4, no. 1 (2017). https://doi.org/10.15347/wjm/2017.002.

Langmead, Ben. “Introduction to the Burrows-Wheeler Transform and FM Index,” n.d., 12. http://www.cs.jhu.edu/~langmea/resources/bwt_fm.pdf 

Miller, Jason R., Sergey Koren, and Granger Sutton. “Assembly Algorithms for Next-Generation Sequencing Data.” Genomics 95, no. 6 (June 2010): 315–27. https://doi.org/10.1016/j.ygeno.2010.03.001.

P43, 45, 51 Link: http://www.cbcb.umd.edu/research/assembly_primer

Schatz, Michael C., Arthur L. Delcher, and Steven L. Salzberg. “Assembly of Large Genomes Using Second-Generation Sequencing.” Genome Research 20, no. 9 (September 2010): 1165–73. https://doi.org/10.1101/gr.101360.109.

Roberts, Richard J., Mauricio O. Carneiro, and Michael C. Schatz. “The Advantages of SMRT Sequencing.” Genome Biology 14, no. 6 (July 3, 2013): 405. https://doi.org/10.1186/gb-2013-14-6-405.
