# How to create a textfile with mothur commands

The easiest way to create a mothur commands text file is to copy / paste your commands into a text file when you are working your way through the MiSeq SOP pipeline, when processing a new sample. However, that takes time. Since mothur is creating a logfile each time one uses mothur, you can use the logfile to extract the commands in one go.

### When you have only one log file present, then run:
```
grep -F "mothur > " mothur.NUMBER.logfile > commands.txt
```
A typical logfile filename looks like this: `mothur.1518098970.logfile`

### When there is multiple log files present, then run:
```
grep -F "mothur > " *.logfile >> commands.txt
```
OR better

```
for file in *.logfile; do
	grep -F "mothur > " $file >> commands.txt;
done
```
Why is the later one nicer?

### cleaning up the commands.txt file
The commands.txt file shows the following 

```
mothur > make.contigs(file=stability.files, processors=2)
```
But is should look like this

```
make.contigs(file=stability.files, processors=2)
```

You can modify the text file using a text editor and then use `search & replace`
or you can use the linux command `sed` like this. 

```
sed -e 's/mothur > //g' commands.txt > commands_edit.txt
```

note that you can write this command also like this:

```
sed -e 's@mothur > @@g' commands.txt > commands_edit.txt
```

In this command you can use any symbol as a seperator between the text to fine and the replacement.

#### Now we can used the commands_edit.txt file as input for mothur and it will then run the command.

You run the file like this:

```
mothur commands_edit.txt
```





