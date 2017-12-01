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

### Finished with the configuration of the Virtualmachine

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




