# Topic 2: Basic commands tutorial

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
<summary> <b>Answer 1</b>  </summary>
  
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
* **nano** : A minimal text editor with the easiest learning curve. (quite using ctrl-x). [Guide here](https://wiki.gentoo.org/wiki/Nano/Basics_Guide)
* **vi** : Keyboard only minimal text editor with a steep learning curve. (quit using :wq). [Guide here](https://www.cs.colostate.edu/helpdocs/vi.html)
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
<summary><b>Answer 2</b>  </summary>

   ```bash
    > echo "This course is fun!" > exercise_2.txt
    > cat exercise_2.txt
    > vi exercise_2.txt
    > mv exercise_2.txt exercise_2_done.txt
```
</details>


