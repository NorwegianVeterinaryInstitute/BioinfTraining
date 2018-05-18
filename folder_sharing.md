# How to set-up a shared folder on a biolinux virtualmachine

We have been working for some time with the biolinux virtualmachine on our Windows laptops.
One thing that makes life a lot easier is when we can share files between the Windows host and the linux guest system.

Below you find instructions on how to set-up a shared folder that can be used to transfer files between the host and the guest system and vice versa. The first part will be done on your windows host and the second part is inside the virtualmachine.

### Setting up the host
* Create a folder on your Windows host that you want to use for sharing files. I call it `VM_shared_folder` and I have it in the directory: `temp` at the root of the `C:`` drive.
* Now allow the folder to be shared. I use a right mouse click and set this folder to be shared.
* Start virtualbox and select the virtualmachine for which you want to set-up a shared folder.
* Open `Settings` and go to the box `Shared Folders`.
* On the right you see a little folder sign with a green `plus`sign on it, click that.
* Add your folder: `VM_shared_folder`, to the `Folder Path`box.
* Click the `Auto-mount` and `Make Permanent` boxes.
* Close the window by clicking `OK`. Your shared folder is now set-up on the host.

Start your Virtualmachine

### Setting up the Virtualmachine

* Make sure you are logged in as the system manager
* open a terminal and type the following
    ```
    cd /media
    ```
    In this folder you will find a folder starting with `sf_` and than the folder name you used. In my case it is:

        sf_VM_shared_folder

    This folder is accessible by the system manager only, so we need to add our own username on the biolinux to the group of people that can access this folder.
* Type on the commandline:
    ```
    sudo adduser YOUR_USER_NAME vboxsf
    ```
* now login with your user account and go the folder: `/media/sf_VM_shared_folder`
and check that you can access the folder. The easiest way to test that sharing is working is to create a text document in this folder.
* the final step is to make a symbolic link from your desktop to the shared folder.
Go to your desktop folder while using the terminal, and then type:
    ```
    ln -s /media/sf_VM_shared_folder ./
    ```

Now you are all set. Happy sharing.
