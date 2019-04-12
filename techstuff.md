# Technical tricks and hits

### Abel cluster command line tricks
* [Using screen to avoid dropped connections](#Using-screen-to-avoid-dropped-connections)
* [Conda virtual environments](#Conda-virtual-environments)
* [Obtaining multiple genomes from NCBI database](#Obtaining-multiple-genomes-from-NCBI-database)
* [Processing multiple datasets in parallel](#Processing-multiple-datasets-in-Parallel)
* [Sharing / downloading data with filesender2](#Downloading-data-from-filesender2)

### Virtualbox biolinux tricks
* [Installing Jupyter notebook and Rstudio](#Installing-Jupyter-notebook-and-Rstudio)
* [Changing your login shell](#Changing-your-login-shell)
* [Sharing data via a shared folder](#Sharing-data-via-a-shared-folder)

### Working with different R in several conda environment with one Rstudio
* [One Rstudio for several conda R environments on your PC](#One-Rstudio-for-all-condas)


## [Using screen](#using-screen) to avoid dropped connections

We are using abel. Once we run a program on abel (or any remote computer), that
program will only keep running a long as the terminal is connected to it, that is,
as long as you have the computer and the ssh connection up. This is a bit
inconvenient since we're working on laptops and might want to pack that one down
without interrupting the connection. There is a program on unix computers that
is called `screen` which will let us do that.

### How to start

1. Log onto abel
2. Type in `screen`. You will get a new prompt, but that's it
3. Start whatever program that yo want to have running.
4. Depress `Control-a`, followed by `d`. Your screen will now say something
about a screen being `detached`.

You can now end your ssh session, and your program will keep running.

### How to see and resume your sessions

You can see a list of screens (you can have multiple), type in `screen -ls`.
You will see something along the lines of `23921.pts-123.login-0-1`. The number
before `pts` is what I refer to below.
To connect to a specific screen, type in `screen -DR <number>`.

Note, since there are multiple login nodes on abel. Thus, if you use screens,
you might want to check which one you're on when you start your screen (use
the command `hostname`). Then, when you log in again, you can check again
which one you're on. If you're on the wrong one, use ssh to log into the
other one.

For more, [check out this link](https://www.tecmint.com/screen-command-examples-to-manage-linux-terminals/)


## Installing Jupyter notebook and Rstudio

We are using the Jupyter notebook and Rstudio within the course. Installation
instructions were sent out in email, but we are repeating them here, to have
them gathered in one place. Note, these are installed inside the virtual machine,
but are on the virtual machine, not abel.

  * [How to install Rstudio](install_Rstudio.md)
  * [How to install Jupyter notebook](install_anaconda.md)


## Changing your login shell

When we use a unix terminal, we have a window with a program inside it. This program
interprets the commands we type in and runs the commands. This program is called a
shell. There are several different shell programs. Ubuntu the way we have it set up
right now, will give you a program called `zsh`. In some cases, it is more beneficial
to have a program that is called `bash` inside our terminals instead. One of these
cases is for running jupyter notebooks. This is because when we installed anaconda,
which is how we installed the notebook program, it put some config options into the
config files for `bash`. Thus if we use `bash`, these options will then be activated.
This means that we should change which shell we have. To do this, run the recipe below.
Note: this is something that you do only once!

1. Open a termininal inside your virtual machine. This should be on the vm, not on abel.
2. Run the following command: `chsh -s /bin/bash`
3. It will ask you for a password, this is your vm password, not your abel password nor
   the manager password.
4. You will then get the prompt again, and it will look like nothing happened.
5. Log in and out again
6. You should now be able to type in `jupyter notebook` in a terminal, and you should
   get the notebook opening in a web browser window.
   
   
## Conda virtual environments
On a computer we can install a lot of different software packages. When you do that on the Abel cluster, or on your own Windows or Macintosh machines, it often happens that software requires additional software in order to function properly. These we call "dependencies". For example, some software requires the python version 2.7 while other software requires python version 3.5 or higher. 

On many computing clusters this is solved by a process called "Sandboxing", where a system is set-up that allows one to load a specific set of software by loading a "module".  The module file contains a list on which software to load to set-up the environment is such a way that you can run for instance the SPAdes assembler. 

Note however, that when you load multiple different modules, it can happen that one version loads python 2.7, while another loads 3.5. At such moments, your software becomes "confused" and tries to run a python script with the wrong python version, and it will not work. In such situations, it can be convenient when you do not have to worry about dependencies having conflicts without having to think about the settings of the system you are running. A way of solving this is to use Virtual environments. The purpose of a virtual environment is to create a space where only software is allowed that does not create internal conflicts due to differences in the dependencies needed. For instance, only python version 2.7 is allowed and not python version 3.5, or vice versa. And if you for some reason need to switch python version, it is only a matter of changing the active environment.

For more on python virtual environments check here: * [Python virtual environments](https://realpython.com/python-virtual-environments-a-primer/)

In recent years using virtual environments has improved and now multiple system exists that helps users to manage virtual environments. On Abel we use the conda system, and see for more here: [Conda virtual environments](https://conda.io/docs/user-guide/overview.html). In order to use this it is needed to set-up the conda system so you can use it to run special software, such as the bifrost pipeline, or ncbi-genome-download.

#### How to set up conda for project nn9305K
In order to access the installed conda environments within the conda project you need to modify a file located in your "home" area on Abel, e.g. the directory where you are when you log into Abel.

1. Check if you have one of the following files in your login directory on abel  (`/usit/abel/u1/YOUR_USERNAME`) :
	
	`.bash_login` or `.bash_profile`.  

	The command to use is: ``ls -a``

2. When you have one of these files (and only _**one**_ is needed !!!) , then open the file with nano to add the lines indicated in point 4, else go to point 3.
3. When you do not have one of these files in your login directory, then copy the example file to your login directory. Make sure you are in your login directory on abel before running the following command:

		rsync /work/projects/nn9305k/samplefiles/bash_login .bash_login
		
4. 	Then add the following lines when they are not present in the file with the command: `nano .bash_login

		export PATH="/work/projects/nn9305k/src/anaconda3/bin:$PATH"
		
		. /work/projects/nn9305k/src/anaconda3/etc/profile.d/conda.sh

5. save and close the file.
6. log out of abel and login again.
7.  To check if conda is present in memory type: `conda info`, that should give output on screen that starts with the line:    

			active environment : None
	The rest of that overview is a summary of the current settings for the conda environment.
	
7. To see which conda environments are present type: `conda env list`.
8. Loading the `ncbidown` environment to download genomes from the NCBI databases:
		
		conda activate ncbidown 
This should give you a prompt that starts with the following:  ```(ncbidown) ```.

9. Deactivate the environment: `conda deactivate`.

10. Now one last check, in some cases it is needed to activate your conda environment in a slightly different manner, so type now on the command line: `source activate ncbidown`. 

	If this again activates the `ncbidown`environment, than your conda set-up was correct. You can deactivate this environment with either: `source deactivate` or `conda deactivate`.
	
Now you are all set to use the different conda environments available on the nn9305k project.

## One Rstudio for several conda R environments on your PC
To avoid multiple installations of Rstudio in several environments, it is possible to define which R (in which conda environment) will be used for interpretation in R studio:

[Rstudio Article](https://support.rstudio.com/hc/en-us/articles/200486138-Changing-R-versions-for-RStudio-desktop)

For ubuntu

```
cd ~
nano .bashrc
# add the following line to .bashrc
export RSTUDIO_WHICH_R="which R"
# save and close .bashrc
# close and restart console or:
source .bashrc
```
You can now activate the conda environment (with the R version and packages you want to work with) `conda activate <env_name>`or `source activate <env_name` depending on your conda version 

From within the activated conda environment start rstudio with `rstudio &`
You can check that the appropriate R version is used in R studio with  `R.version.string`


## Sharing data via a shared folder
We have been working for some time with the biolinux virtualmachine on our Windows laptops. One thing that makes life a lot easier is when we can share files between the Windows host and the linux guest system.

  * [Setting up a shared folder](folder_sharing.md)

## Processing multiple datasets in parallel

### The problem
You have  several hundred datasets that you want to analyze using a method of choice. Examples of these are:
* Assembly of shotgun datasets from several hundred bacterial isolates
* Hundreds of metagenomic samples that need to be quality controlled and classified.
* Iteration of parameters for testing a model you are developing.
* etc...

[Processing multiple jobs in parallel](array_jobs_abel.md)

## Obtaining multiple genomes from NCBI database

Downloading one genome from the NCBI database is relatively easy and can be done with help of a webbrowser and the following webpages: 

* [https://www.ncbi.nlm.nih.gov/genome/microbes/](https://www.ncbi.nlm.nih.gov/genome/microbes/)
* [https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/](https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/)

However, in order to download multiple genome from the NCBI webpages takes more effort. See also the section above on [Conda virtual environments](#Conda-virtual-environments). Conda needs to be set-up correctly for this to work.

* [How to download multiple genomes](genome_downloads.md)

## Downloading data from filesender2

In a lot of projects it is important to be able to share sequence data. The University of Oslo has a system in place that can help you share data between yourself and collaborators in the rest of the world. You find it at: [https://filesender2.uio.no/](https://filesender2.uio.no/). It is very easy to use and the filesize seems unlimited. What is not so easy to do is to get the data that was uploaded to filesender, downloaded straight to Abel.

Here is my recipe for when people send me large sequence data and I want to use it on Abel.

1. Login to filesender2 and send out a request for data, unclick the box to allow for more uploads per link (in case your collaborators fail to put everything in one tar file.
2. Collaborator uploads one or more files to the fileserver2 and you get a notification in your mailbox.
3. In the notification mail that data is uploaded, Follow the link to the data for download in your webbrowser. It opens a page on the filesender webpage. 
4. Go to the file you want to download and "right" click the download box and copy the link address.
5. In your Terminal, go to abel and the location where you want to store the data. Somewhere on `/projects/nn9305k/`:
6. Then type in the terminal: `wget -O myfile.data` and paste the link address from 4 at the end but within quotations marks. like this: ``` wget -O large_genome.gz "YOURLINK" ``` 
7. Than the data is downloaded to a file with the filename given in 6.

Note, that you can download multiple file in one go as a zip file by scrolling to the bottom of the file list and copy the link produced by the box: `Download as single (.zip) file`




