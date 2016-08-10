# Topic 2: Introduction to unix command line and R

### Navigating the directory structure
There are a few basic commands needed to move around the unix directory structure
* **pwd** : _Print working directory_. This tells you where you currently are in the directory structure.
* **mkdir** : _Make directory_. This makes a new directory where you specify.
* **cd** : _Change directory_. This moves you to a new directory.
* **ls** : _List directory contents_. This tells you what is in a directory. It has options to order or include more information. A good starting command is "ls -thor"
* **rmdir** : _Remove directory_. Removes an empty directory.
* **man** : _Show manual_. Lets you view the manual for the command.

Now lets try it out.
```bash
> pwd
/home/owens
> mkdir cats
> ls
cats
> cd cats
> pwd
/home/owens/cats
```
**Exercise 1**: 
* Make a new directory in your home directory named "exercise_1". 
* Inside that directory make another directory named "inside_1". 
* Move from your home directory directly to the inside_1 directory.
* Move back to your home directory
* From your home directory look at the directory contents of inside_1.
* Remove both new directories.


<details> 
  <summary>**Answer 1**  </summary>
   ```bash
    > mkdir exercise_1
    > mkdir exercise_1/inside_1
    > cd exercise_1/inside_1
    > cd /home
    > ls exercise_1/inside_1
    > rmdir exercise_1/inside_1
    > rmdir exercise_1
```
</details>

### Writing and viewing text files
* **echo** : Prints text you give it to STDOUT. 
* **cat** : Prints all text of a file to STDOUT.
* **less** : Lets you view text line by line, without loading the whole file (quit using q).
* **vi** : Keyboard only minimal text editor (quit using :wq). [Guide here](https://www.cs.colostate.edu/helpdocs/vi.html)
* **emacs** : Another Keyboard only minimal text editor, but more expandable) (quit using ctrl-x, ctrl-c). [Guide here](http://mally.stanford.edu/~sr/computing/emacs.html)
* **seq** : Prints sequential numbers up to your specified value.
* **mv** : _Move_. Moves a file, also used for renaming.
* **rm** : _Remove_. Removes a file.
Lets see how they work.

```bash
> echo "Bioinformatics is fun!" 
> echo "Bioinformatics fun file!" > tmp.txt
> cat tmp.txt
> seq 100
> seq 100 > numbers.txt
> less numbers.txt
> vi numbers.txt
> emacs numbers.txt
> mv numbers.txt 100numbers.txt
> rm 100numbers.txt
```

**Exercise 2**: 
* Make a new file named exercise_2.txt with the text "This course is fun!".
* Print exercise_2.txt to the screen.
* Open exercise_2.txt in a text editor and add a second line that says "Text editing is awesome."
* Rename exercise_2.txt to exercise_2_done.txt.


<details> 
  <summary>**Answer 2**  </summary>
   ```bash
    > echo "This course is fun!" > exercise_2.txt
    > cat exercise_2.txt
    > vi exercise_2.txt
    > mv exercise_2.txt exercise_2_done.txt
```
</details>

### Pipes and redirection
A key feature of command line use is piping the output of one command to the input of another command. This means that large files can be analyzed in multiple scripts without having to write to disk repeatedly. 
##### *Key terms*
* **STDOUT** : Standard out. The regular output of a script. Can be directed to a file with ">" or "1>", and directed to another command with "|" .
* **STDIN** : Standard in. The input of STDOUT, piped in from another command.
* **STDERR** : Standard error. The error output of a script. Generally prints to screen but can be saved to a file using "2>".

##### *sed*
*Stream editor*. It parses and transforms text using regular expressions. Very powerful, but most easily used to reformat files based on patterns.\
**Examples**: 
* Replace all instances of "1" with "one". 
  * seq 10 | sed s/1/one/g
* Replace lines that only have "1" with "one".
    * seq 10 | sed s/^1$/one/g
 ##### *grep*
*Search using regular expression*. This command searches for patterns and prints lines that match your pattern.\
**Examples**:
* Print all lines with "1".
    * seq 10 | grep 1
* Print all lines without "1".
    * seq 10 | grep -v 1
* Print all lines with "1" and "0"
    * seq 10 | grep 1 | grep 0
* Print all lines with "1" or "2"
    * seq 10 | grep "1\|2" 

**Exercise 3**:
* Print the even numbers up to 100.
* Remove all numbers divisible by 10.
* Add "!" after every number ending in 2.
* Print only numbers with "!" or "3".
* Save the resulting file to exercise_3.txt

<details> 
  <summary>**Answer 2**  </summary>
   ```bash
    > seq 2 2 100 | grep -v 0 | sed "s/2$/2\!/g" | grep '\!\|3' > exercise_3.txt
```
</details>

### Running commands in background
Often you will run commands that take hours or days to finish. If you run them normally your connection needs to be maintained for the whole time which can be impossible. Using _screen_ allows you to keep a screen open while you're logged out and then reconnect to it without loss of data at a later time. 

**Cancel command** = ctrl-c. This will cancel a script that is currently running.
Example: 
```bash
> seq 1000000
ctrl-c to cancel
```
##### *screen* commands:
* **screen** : Opens a new screen session. Only use when you don't have any running screens.
* **screen -r** : Reattaches a detached screen session. Use this when re-entering a previously used screen.
* **screen -d** : Detaches a currently attached screen you aren't using now. Use this when you forgot to formally detach from your screen last time you used it.
* **Commands within screen**
    * ctrl-a, c : Opens a new screen window.
    * ctrl-a, n : Moves to the next screen window.
    * ctrl-a, p : Moves to the previous screen window.
    * ctrl-a, d : Detaches from screen session.
    * exit : Closes current screen window.

**Exercise 4**:
* Open a new screen session.
* Make a new window.
* In the new window print numbers 1 to 10000000 
* Move back to your old screen.
* Periodically check on the number screen to check when it is done.
* When counting is done, close the original empty window.
* Detach from the screen session.
* Reattach to the same screen session
* Close the screen session entirely
 
<details> 
  <summary>**Answer 4**  </summary>
   ```bash
   > screen 
   ctrl-a, c
   > seq 10000000
   ctrl-a, n
   ctrl-a, p
   > exit
   ctrl-a, d
   > screen -r
   > exit
 ```
</details>







