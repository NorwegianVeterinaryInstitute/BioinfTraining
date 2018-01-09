Setting up the Linux Virtual machine for the lecture: `Programming with R`

Start the Linux Virtual Machine (VM) and make sure your are logged in as a "user" and not as the "administrator". The Linux VM will automatically log you in as a administrator, when you start it up, so when that is the case, logout, and login again with your username and password.

Open-up Firefox in the Linux VM and copy the link below into the url address box and follow the instructions below.

* Download the correct version of R-studio for the Linux ubuntu version.
Copy this link to firefox: [https://download1.rstudio.org/rstudio-1.1.383-amd64.deb] (https://download1.rstudio.org/rstudio-1.1.383-amd64.deb)

* Open with Ubuntu Software centre (wait until you get a window saying: rstudio )
 
* Click install, this will set-up R-studio.

Next we need to download the dataset needed for this lecture:

* Create on your Linux Virtual Machine a folder called: `r-novice-inflammation`
* Next copy this link to the Firefox url box in your Linux VM: [r-novice-inflammation-data.zip] (http://swcarpentry.github.io/r-novice-inflammation/files/r-novice-inflammation-data.zip)
* the file needs to be saved in the folder you just created.
* We then unzip the downloaded file. Open a terminal and move to the desktop:
	`cd Desktop`
* move to the folder `r-novice-inflammation`
* unzip the file:  `r-novice-inflammation-data.zip`  with the command: `unzip r-novice-inflammation-data.zip`
* check with `ls` that a folder called: `data` is created.

Now we are set for the lecture.



##Acknowledgement

This lecture was prepared following the lecture from software carpentry:

* [R for reproducible Scientific Analysis](http://swcarpentry.github.io/r-novice-gapminder/)

