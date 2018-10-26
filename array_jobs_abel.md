#Processing multiple datasets in Parallel.

There are many situations where you have many samples and you can make a script that processes them one after the other. When you have a stand-alone computer that is the only option you have. But doing for instance bacterial genome assemblies for several hundred isolates will take a lot of time, even with a fast computer.

However, when you have access to a Computing cluster such as [Abel](https://www.uio.no/english/services/it/research/hpc/abel/), then you can harness the power of the cluster by spreading out your assembly jobs, by running the assemblies in parallel. That substantially reduces the waiting time for yourself.

#### NOTE
Each user on Abel can only use 400 cpu's simultaneously, that is 200 jobs using 2 CPUs per job. If you want to run more jobs, any job after 400 cpu's are started has to wait until one of the earlier jobs has finished.

### The solutions.

There are three possible methods where you can use parellel processing to analyze many samples.

1. Start slurm jobs with a ``for loop``, where each iteration uses a different dataset.
2. Use the SBATCH command: ``#SBATCH --array=0-99``
3. Use the ``Arrayrun`` command.

Below you will find an explanation of the three methods and a link to a small example script to run fastqc on your illumina sequence data. There are also pros

### 1. Start slurm jobs with a "for loop"
This method is the easiest to setup. What is needed is that you create a slurm job that can take an input file.  Inside the script that input file is called using the following variable ``${1}``.

So then the command to start the slurm jobs becomes:

```
  for file in *.fastq; do
    sbatch fastqc_script.slurm $file
  done
```

When you have X files than this loop will request X jobs in one go. This is good when you have few files to analyze. There is however a little danger here, when using a many, many files it could be that you load in one go many jobs onto Abel, that will keep other peoples jobs from starting. That is not social behaviour. So think about it.

Example script: [for_loop_example_script.slurm](array_example_scripts/for_loop_example_script.slurm)

### 2. Use #SBATCH --array

This Sbatch command is added to the slurm script, and it will start as many jobs as requested. For instance:

```
 #SBATCH --array=0-89
```
This will start ``90 jobs`` in the following manner: request 10 jobs and then sleep several minutes. Then check if all 10 have started. If so request another 10. If not, request as many so that only 10 jobs are waiting. It will continue with this until all 90 jobs are started. This a more social behaviour to request resources on Abel.

```
Q: What happens if the number of requested jobs is not equal to the number of input datasets?
```

Inside the Slurm script you define which datasets are analyzed for each job that is requested. A possible method on how to do that is shown here:

```
#Generate array of unique datasets
FILES=($(ls -1))

# counting number of datasets
echo "Number of samples:" ${#FILES[*]}

#Generate name of current sample
#[$SLURM_ARRAY_TTASK_ID] is the element index. Increases with one per job, i.e. one new sample for every job)
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Current input-file name is" $CURRENT_FILE
```
With the variable $CURRENT_FILE you can then start for instance a job to run Fastqc

Here is an example script: [sbatch_array_example_script.slurm](array_example_scripts/sbatch_array_example_script.slurm).

### 3. Use the ``Arrayrun`` command

This latest command is used in a ``master``script and it will start an X amount jobs using a ``worker``script to analyze a dataset. If functions in a similar way as the previous example, but the set-up is slightly different.

For instance:
```
# running the arrayrun command with tasks id's ranging from 0-MAX_ID
arrayrun 0-$MAX_ID workscript_fastqc_example.slurm
```

In place of ``$MAX_ID``, you could also define a number like in the ``SBATCH --array`` script, however, when using the variable ``$MAX_ID`` the master script automatically knows how many datasets to analyze.

Here is how you set it up:

**MASTER script:**

```
# collecting all the PAIRED dataset files into an array called FOLDERS
FILES=($(ls -1))

# Than we count the number of files to determine the number of tasks
NUM_TASKS=${#FILES[*]}

# determine the maximum number of tasks
MAX_ID=$((NUM_TASKS - 1))

# running the arrayrun command with tasks id's ranging from 0-MAX_ID
arrayrun 0-$MAX_ID ../slurm_scripts/workscript_fastqc_example.slurm
```

The worker script also needs to "know" which datasets to analyze. So there you also need to identify which dataset is used in the current job that is started.

**WORKER script:**

```
# The variable FILES contains all the files for analysis
FILES=($(ls -1))

#create a variable using the filename and combine it with the TASK_ID from arrayrun
CURRENT_FILE=${FILES[$TASK_ID]}

echo current dataset is $CURRENT_FILE
```

This set-up allows one to analyze an unknown amount of datasets. At least you do not have to determine it yourself.

Here are the example files:

* [array_submission_fastqc_example.slurm](array_example_scripts/array_submission_fastqc_example.slurm)
* [workscript_fastqc_example.slurm](array_example_scripts/workscript_fastqc_example.slurm)
