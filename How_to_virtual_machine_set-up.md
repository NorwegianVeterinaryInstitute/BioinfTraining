## Setting up and running the Biolinux Virtual machine on a windows 10 laptop.

The bio-linux version we are going to use is: `Bio-linux version 8.0.7.` and below you find the steps to set-up the bio-linux virtualmachine correctly.


#### *This first section (steps 1-5) is only needed as a control to check your hardware set-up*

1. start-up your computer and press `F1` as soon as possible. Now we enter the BIOS system of our laptop.
2. Go to section: `Virtualization`
3. Check if the option `Virtualization` is enabled.  
4. Check if the option `VT for Direct I/O`  is enabled.
5.  When both options are enabled, you can press `F10` to save and exit the bios. 
6. Then normal start-up is proceeded.

#### *Installation with a normal start-up of the laptop*
1.  Login to your computer and create on the Data drive `D:` a folder called: `virtual_machines`
2.  When you receive  one of the memory sticks copy the two files to the directory: `virtual_machines`
		
		Files: `bio-linux-8.0.7.ova`  &  `IntroductionToBio-Linux8_Dec2015.pdf`
3. Go to the desktop and click the `Oracle VM virtual box` icon

	This opens up a window that is the main controller window for virtual machines.
4. select "File" from the top menu, and then "Import Appliance"
	
	A window opens, and click the folder box, which opens a window to select files. 
5. Go to disk `Data D:` and into the folder: `virtual_machines`. 
6. Now select the file `bio-linux-8.0.7.ova`, and click `open`.
7.  In the import window click: `Next`

	Now we are at the appliance settings, which we are going to increase to make better use of your laptops hardware, but leave room for the windows systems to function normally.
8. Now click `import`, this will use the image to create the virtualmachine disk on your harddisk.
	
	After this we are back at the VM main controller window. 
9. Select the Bio-linux-8.0.7, by clicking once. Than click the yellow icon "settings".

	This opens up a window, to change the configuration of your virtual machine. The section we see is called "General"
10. Click the tab: `Advanced`
11.  Set `Shared Clipboard` to: `Bidirectional`

    Now we can copy something from the windows machine and drag it the Linux machine, and vice versa
    
12. Now click the box called `System`. We now see three tabs `Motherboard, Processor` and `Acceleration`.
	
	`Motherboard` includes how much of the memory from the laptop we can reserve for the Linux Virtualmachine. The amount of RAM memory on your laptop is `32 Gb`, the maximum is indicated at `Base-memory`.
13. We increase the `Base-memory` to the maximum of the green bar. (I estimate it to be between 24000 MB and 28000 MB)
14. Click the tab `Processor`
15. Your computer is equipped with 4 cpu's or cores. We set the cpu to `2`. 
16. Than we click `Ok` and we are done with the configuration.

You are now inished with the configuration of the virtual machine

_______
_______


## Start-up of the installed Bio-linux virtualmachine
This section explains how to start-up the bio-linux machine and adjust settings so you can work with a norwegian keyboard outlay. In addition does it point you to the documentation that comes with the bio-linux virtual machine.

In the VM main controller window you see the virtual machines installed on your laptop.

1. Double click the Bio-linux virtualmachine, or click the green arrow button called `start`.
	
	`Now the virtual machine is activated and running.`

	The first window that pops up, is an update window. 

2. For now we do not upgrade the ubuntu version, and click `Don't upgrade` or `Ask me later`

	The Virtualbox program also shows two pop-ups. 	
	
3. The top one tells you about the `Auto Capture Keyboard` option. Dismiss that box by clicking the cross on the right.
4. The bottom one explains about the `mouse pointer integration`. Dismiss also.
5. Look at the top-right of the image. There is an icon called `En` That is the current keyboard layout. Click the button.
6.  move to `Text entry setting`. Click on the `+` sign, search for `Norwegian`, and add that to the list.
7. Close the window. Now click again on `En` and change that to `Norwegian`. 
	
	Now we can use a Norwegian keyboard.

	Next we want to explore the Virtualmachine before we start working with this "Linux" computer.

	On the the desktop you see two icons.

	```
	– Bio-Linux documentation  
	– Sample data (This contains sample data for a lot of the software pre-installed on this virtual linux computer.)

	```

8. Double click the Bio-Linux documentation. 
9. Next double click the introductory tutorial icon.

	Now you have a readme file, a PDF with an introduction to bio-linux, and some archive files (end at tar.xz) with data for different tutorials.

10. open the pdf: `introductionToBioLinux8_jan2015.pdf.`

	We will go through the first few pages to familiarize ourselfs with the bio-linux operating systems. (Pages 1 to 7, until Root-directory)


	Now you are ready to start exploring the Linux computer and now we can explain about the file-system on linux computers.


## Creating a normal account on your virtual machine

The acount you are now using on your virtual machine is the system administrator
account. This user has permission to pretty much anything. Thus, using that account
for your work could become risky if you're not complettely sure of what you're
doing. Thus we'll have you make another one.

1. Click on the 'Settings' button in the Dash. You will see a 'Users' icon, click
on that.
2. On the window that now appears, click 'Unlock' in the upper right corner. You
will be asked for a password to do this, this password is 'manager'.
3. Click on the small plus sign under the square on the left side. You should now
get up a new window.
4. Select account type 'Standard', and write in your full name. You will next 
select a username. If you have had an account on the UiO system at any point, 
please use that. If not, create a username that combines some letters from your
first and your last name.
5. Click 'Add'
6. While your new acount is selected in the menu, click on the password field on
the right side, where it says 'Account disabled'. You will now get up a new window
where you will fill in a password. Please: do not use your vetinst password or 
any other existing password. Also, do not use Norwegian characters or special 
characters. Click 'Change'.
7. Log out from the system manager account by clicking on the 'Gears' (settings) 
icon in the upper right corner, and select 'Log out'.
8. You will now get a login screen. Select your new account, and type in your
password, and you're ready to go.


