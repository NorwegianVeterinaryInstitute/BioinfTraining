### How to setup and interactive job on abel
In order to run an interactive job on abel it is good to first start the little linux program `screen`. After starting screen we start-up a qlogin to request computing time.

### Steps to run an interactive job on abel

* start a `screen`(see below on how to do that)
* Next ask for a qlogin, for instance:
  ```
  qlogin --account=nn9305k --mem-per-cpu=3800M --cpus-per-task=4
   --time=10:0:0
  ```
* once the job is accepted `source`it to active it
  ```
  source /cluster/bin/jobsetup
  ```
* now you can load module, run programs etc...
* if you need to detach (you might want to go home, but the program is still busy) then do:
  ```
  ctrl-a d
  ```
* reattach by using the command: `screen -r `
* if the program is done that use the command `exit`to stop the qlogin and then again `exit` to deactivate the screen.


### A few basic commands for using screen:
To start a screen type on the commandline:
```
screen
```
detaching from a running screen:
```
ctrl-a d
```

re-attaching to a running screen:
```
screen -r
```
check if there are screens running
```
screen -ls
```
this might give you an output that looks a bit like this:
```
There are screens on:
        2477.pts-0.server1      (Detached)
        2522.pts-0.server1      (Detached)
2 Sockets in /var/run/screen/S-root.
```
then re-attach using the command:
```
screen -r 2477.pts-0.server1
```

and finally stop a running screen by typing:
```
exit
```

Here is a nice [tutorial](https://www.howtoforge.com/linux_screen) on using screen.
