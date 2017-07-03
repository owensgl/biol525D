# Topic 2: Commandline and R
In this tutorial we will be going through parts of several software carpentry workshops.  Right now, just code along with me, but after the lesson you can go through the tutorial on your own to clarify any understanding problems. 
## [TUTORIAL](http://swcarpentry.github.io/shell-novice/)

### Other Tutorials
Here are some good tutorials if you're interested in learning a programming language
* [Python](https://www.codecademy.com/learn/python)
* [Perl](http://www.perl.com/pub/2000/10/begperl1.html)
* [R](http://swirlstats.com/)

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

**Exercise 1**:
* Print the even numbers up to 100.
* Remove all numbers divisible by 10.
* Add "!" after every number ending in 2.
* Print only numbers with "!" or "3".
* Save the resulting file to exercise_3.txt

<details> 
<summary>
<b>Answer 1</b>
</summary>

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
#### *Byobu*:
[Guide to Byobu](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-byobu-for-terminal-management-on-ubuntu-16-04)

First you have to install byobu on your system. In ubuntu you can use the following command:
```bash
>sudo apt-get install byobu
```
Byobu can create multiple levels.
* **Session**: A running instance of byobu. You can have multiple of these and when you start byobu you select which session you want to run. You can also switch between sessions. Sessions will continue existing and running on your computer until you shut them down. You may want multiple sessions if you connect to with different screen sizes. 
* **Window** : A session can have multiple windows. You can easily toggle between windows using F3 and F4. If you start a command in a window and then detach the session or switch windows, the command will continue running. Generally when you are working, you will have multiple windows open for different tasks (e.g. testing a script, editing that script, looking for files).
* **Panes** : A window can have multiple panes. Panes split your window into multiple panes. These are functionally windows, but exist together on your screen. Useful if you want to observe multiple things at one (e.g. watch cpu usage while running a script).

#### Commands in Byobu
* **byobu** : Opens byobu and attaches a session. If you have multiple sessions you will have to select which session to attach.
* **F2** : Creates a new window.
* **F3** : Toggles through your windows.
* **F8** : Renames the current open window in the list.
* **F7** : Lets you view scrollback history in the current window.
* **SHIFT+F2** :  Creates a horizontal pane.
* **CTRL+F2** : Creates a vertical pane.

**Exercise 2**:
* Open a new byobu session.
* Make a new window.
* In the new window print numbers 1 to 10000000 
* Move back to your old window.
* Periodically check on the number screen to check when it is done.
* When counting is done, close the original empty window.
* Detach from the session.
* Reattach to the same byobu session
* Close the byobu session entirely
 
<details> 
<summary><b>Answer 2</b></summary>
  
   ```bash
   > byobu 
   F2
   > seq 10000000
   F3
   F3
   > exit
   F6
   > byobu
   > exit
 ```
 
</details>

 

### Daily Assignments
1. What is one task you'd rather use an R script instead of a shell script? Why? What is one task you'd rather use a shell script, instead of an R script? Why?
2. Why is piping directly between programs faster than writing each consecutive output to the disk? Explain using information about computer hardware.




