#How to troubleshoot a new program:
* Are all the parameters and files spelled correctly? Look for typos.
* Does it work with example data? If not the program wasn't installed sucessfully.
* Does my input data look like the example data? 
* Can I make a minimum dataset that works? Start with just one site and two individuals. Then when that works, keep adding data until it stops working or you fix it.

#How to test a shell script as you build it:
* Run with the minimum number of samples. Consider taking a subset of reads for your samples so they run extra fast.
* Use "exit" to partition off parts of your script to run one at a time.	
* Comment out things you don't want to rerun after you have made them work.
* Use if commands to only run commands you need. Make the shell check to see if the output file already exists and then don't run the command if that is so.
