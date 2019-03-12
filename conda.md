# Using software installed in conda

## What is conda?

Conda is a package (=software) and environment management system. An environment
is a container of package(s).

Think of an environment as a bubble. When you are in your bubble you do not to
be disturbed by anything outside your bubble.

Conda allows the creation of several environments, that do not interact with each
other. Please read the explanations in [techstuff](techstuff.md#conda-virtual-environments)

![alt text](/figures/conda.svg)

**In other words**:
Creating an environment name=X1 to run a particular software is extremely useful when
the software you want to work with depends (=a dependency) on one or several other
softwares=Y of a specific version to function.

If you want to use a software that depends on other versions of software=Y you can install those
in different environments (eg. name=X2). Then softwares installed in each  won't interfere
with each other.

There are several several packages installed within several conda enviroments for you to work with on Abel. See below.

> PS: be sure that you did follow at first time login the instructions: "[On first time login](https://github.com/NorwegianVeterinaryInstitute/organizational/wiki/Abel-User-Guide)" as it automatically configure of where conda is to be found. If you did it - you do not need to do it again.

## Working with softwares installed in conda environments

- activating conda: `conda activate`
- viewing the environment list: `conda env list`
- listing packages installed in a specific environment: `conda list -n <environment_name>`
  > allows you to see which version of software are installed in each environment)
- activating an environment: `source activate <environment_name>` : enables working with softwares within the activated environment
  > if you want to look at the packages _while the environment is activated_ you need to use: `conda list` the list can be long: if you want look only at part of the list: `conda list <you_want_to_look_at>` ex. `conda list trim*` give me the packages installed that begins with trim
- deactivating an environment (or deactivating conda) `conda deactivate`
  > if an environment is activated: you need to **deactivate twice**: first to deactivate the environment and then to deactivate conda
- getting help: `conda --help`

1. request for resources (if slurm.script > steps below in the script)
2. activate conda and environment
3. your commands (programm you want to use with parameters)
4. deactivate conda and environment

> Oops: bifrost pipeline: ask for resources when you submit the running.scripts.

## Installing packages on conda: **contact Karin**

> if packages you want to work with are neither available on Abel nor conda:
> you need to contact Karin

## Going further

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/index.html)

Links provided in [techstuff](techstuff.md#conda-virtual-environments)
