# Topic 2: Introduction to unix command line and R

### Navigating the directory structure
There are a few basic commands needed to move around the unix directory structure
* *pwd* : _Print working directory_. This tells you where you currently are in the directory structure.
* *mkdir* : _Make directory_. This makes a new directory where you specify.
* *cd* : _Change directory_. This moves you to a new directory.
* *ls* : _List directory contents_. This tells you what is in a directory.

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


<details> 
  <summary>**Answer 1**  </summary>
   ```bash
    > mkdir exercise_1
    > mkdir exercise_1/inside_1
    > cd exercise_1/inside_1
    ```
</details>




