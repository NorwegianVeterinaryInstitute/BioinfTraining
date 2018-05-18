# Technical tricks and hits


## Using screen to avoid dropped connections

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

## sharing data via a shared folder
We have been working for some time with the biolinux virtualmachine on our Windows laptops. One thing that makes life a lot easier is when we can share files between the Windows host and the linux guest system.

  * [Setting up a shared folder](folder_sharing.md)
