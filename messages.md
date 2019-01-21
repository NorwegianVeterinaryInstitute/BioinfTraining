# Messages

This webpage contains information from outgoing messages, such as homework,
practical information, etc.

The information here is organized in reverse chronological order.

Information that will not be included:

 * sensitive information
 * names of people
 * email addresses
 * locations of files on the internal network
 
## 2017-12-01 Homework

  * You can find almost all information about the course 
  [on this webpage](https://norwegianveterinaryinstitute.github.io/BioinfTraining/).
  Please check the Messages part to find this week's information/homework.
  * For the unix part, we got to the "Redirecting Input" part of the lesson.
  * You will need to create an account on the Abel computing cluster. To
  do that, you need to go [to this link and fill out the form that is
  mentioned there](https://github.com/NorwegianVeterinaryInstitute/organizational/wiki/Abel-User-Guide). Note: use
  the same username as you chose for your account on your virtualbox.<br/>
  **DEADLINE: Tuesday 5th.**
  * Regarding Thursday 7th: please [fill out this doodle to choose which
  session you will be attending](https://doodle.com/poll/sxt5hmub94nhegyv). <br/>
  **DEADLINE: Tuesday 5th.** 
  Note, if you choose to attend either the 11th or the 12, you commit to having
  completed the entire unix exercise on your own
  * Please start finding papers and figures and examples of what you will be 
  wanting to be able to do after the course is done. We will use that to figure
  out what to teach, but also for teaming people up.
  
## 2017-12-14 Homework

  * Refer to the [session's page](data_pre_processing.md) for more information regarding the homework (and example SLURM scripts)
  * Use trimmomatic to trim/remove adapters and low quality reads in _Br_R1.fq.gz_ and _Br_R2.fq.gz_ (or/and _Ed_R1.fq.gz_ and _Ed_R2.fq.gz_)<br/>
  **DEADLINE: Tuesday Jan 9th.**

1. Remember that you are working with paired end data (Change _SE_ to _PE_). 
2. There are two input files and four output files.
3. Use _TruSeq3-PE-2.fa_ instead of _TruSeq3-SE.fa_ since we are dealing with paired end reads.
4. Change _MINLEN_ parameter to _36_.
5. Use appropriate value for _CROP_ (check the fastqc output for raw reads and use the correct value).
6. Also, remember to change _#SBATCH --mem-per-cpu=4Gb_ to _#SBATCH --mem-per-cpu=12Gb_. This is a bigger job and needs more memory (12Gb instead of 4Gb).
