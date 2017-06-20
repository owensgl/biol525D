# Topic 2: Commandline and R
In this tutorial we will be going through parts of several software carpentry workshops. The tutorial can be found [here](http://swcarpentry.github.io/shell-novice/). Right now, just code along with me, but after the lesson you can go through the tutorial on your own to clarify any understanding problems. 

### Tutorials
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
  <summary>**Answer 1**  </summary>
   ```bash
    > seq 2 2 100 | grep -v 0 | sed "s/2$/2\!/g" | grep '\!\|3' > exercise_3.txt
   ```
</details>
<details>
<summary>CocoaPods</summary>
Add the following line to your `Podfile`:

```ruby 
pod 'YourAwesomeLibrary'
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

**Exercise 2**:
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
  <summary>**Answer 2**  </summary>
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

#### *Byobu*
This program is similar to screen, but has some advantages over screen especially in keeping track of windows.
You can find it [here](http://byobu.co/)


Lastly, we will go through a [short tutorial](http://swcarpentry.github.io/r-novice-gapminder/01-rstudio-intro/) on R. 

###Daily Assignments
1. What is one task you'd rather use an R script instead of a shell script? Why? What is one task you'd rather use a shell script, instead of an R script? Why?
2. Why is piping directly between programs faster than writing each consecutive output to the disk? Explain using information about computer hardware.




