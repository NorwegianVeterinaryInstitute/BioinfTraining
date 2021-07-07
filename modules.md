
> put in clear notes from Thomas presentation 16/09/2019 On making our own modules

You can create modules in your own work folder on Abel: $HOME equivalent to `/usit/abel/u1/username`

NB: to see the path represented by the $HOME variable do `echo $HOME`

example:
```
cd $HOME
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
./seqk # to launch from the directory
```
To be able to load those modules you have to options

1) Add the path where the software is to your login script `.bash_login`

$PATH:...

2) create a module for this software and load the module when the software is needed

# Create and use your own modules on Abel 

1) create a directory `modules` on your own area: $HOME
```
cd modules
ls
```
All the files listed are text files describing of where the computing ressource will find the the softwares that are installed in modules

A "module file" contains details about applications, description of what the module is, commands that will be launched when you load this module file: in other words: everything that will be loaded and added to the $PATH
> What is the $PATH ? its where the computer looks for softwares to exectute ...documents? 

Eve: NB: see in the HTP advancesd course -> now remember was things about modules

2) in your .bash_login:

Export path of softwares and modules we can use. 


Module environment 

3) Launching module - identical as Abel own modules

`module launch module_name/version`

> version if you want to use a version different than the default version

# Links
[environment modules documentation](https://modules.readthedocs.io/en/latest/index.html)

> for those with experimental ideas ... can do on ubuntu PC -> will use that to test first
