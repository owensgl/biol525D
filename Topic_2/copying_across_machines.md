---
title: "Topic 2: Copying files across machines"
---

These instructions are part of [Topic 2](./#copying-files)

### Copying files between machines

You can copy files from your computer to your server account and (vice versa) in a multitude of ways.

Follow instructions that pertain to your software and operating system.

With the terminal (linux/mac/mobaxterm):
-------------------

**To copy a small set of files you can use the `scp` command.**

   > *Note:* If your `~/.ssh/config` has an alias for serveruser@serverhost, called b525, substitute "serveruser@serverhost" with just "b525" below

   Examples:

       # creates a copy of three files in ~/path/on/server/
	   # the destination directory must exist
       scp ~/mydir/file1 ./file2 /path/to/file3 serveruser@serverhost:path/on/server/

       # creates a copy of the file in the home directory on the server (nothing after :)
	   scp myfile.txt serveruser@serverhost:

       # the server path is relative to the home directory, but an absolute path can be
	   # provided
	   scp myfile.txt serveruser@serverhost:/scratch/myuser/

       # for Moba, your C:\ drive is in /mnt/c/.
	   # moba also creates a folder called MyDocuments inside your moba home
	   # which points to your windows documents folder.
	   scp ~/MyDocuments/my_file.txt serveruser@serverhost:path/on/server/

You can reverse the arguments to copy files from the server to your laptop. The last argument is always the destination.

**To copy a large set of files (or a small set of large files)**

For large files, or synchronizing content between computers, best use
`rsync`, which can skip files that are already on the server side
(based on timestamps). There are also options to transfer based on
content differences (see `man rsync`).

To copy a whole local directory `localdir` recursively on the server, into `~/dest/`,
you would use:

    # Note: the destination has a trailing '/'. This means "into dest/".
	# X/1.txt ends up in ~/dest/X/1.txt
    rsync -v --progress -rlpt localdir serveruser@serverhost:dest/

    # Note: the destination has no trailing '/'. This means "X" is renamed to "dest"
	# X/1.txt ends up in ~/dest/1.txt
	rsync -v --progress -rlpt localdir serveruser@serverhost:dest

Options used in the example:

   - `-r`: recursive
   - `-l`: copy symlinks as symlinks
   - `-p`: keep permissions (e.g. executables remain executables, read-only remain read-only, etc.)
   - `-t`: preserve timestamps
   - `-v`: verbose

Using cyberduck (mac/windows)
--------------


   If you connect to your assigned server using SFTP, you can browse to the destination folder and use drag and drop.
   Take care to save your file with Unix line endings if you can configure your editor to do so.

1. Using the sftp browser in MobaXterm

   Once you're connected in a user session, on the left there is a pane showing your files. You can use drag and drop
   and create new files there.

Using Psftp (from the Putty suite for windows)
-----------------

Start psftp.exe (<kbd>windowskey</kbd>+<kbd>r</kbd> `psftp` <kbd>enter</kbd>).

The commands below are typed at the psftp prompt:

Connect to the server with the `open` command:

    open serverhost

> *Note*: If you have a saved putty profile for the course, e.g. `b525`, you can simply do:
>
>             open b525
>
> and your username, host, and your preset key/password settings will be loaded from the named profile.

With the SFTP protocol, you adjust the source folder and destination folder independently, then you issue
either a PUT to upload a local file (to the remote), or a GET to download a remote file (to the local computer).

1. Change the local folder with the `lcd` command. eg. for the user `myself`'s My Documents:

        > lcd c:\users\myself\documents

1. Change the folder on the remote server with `cd`. e.g. for folder `~/my_files/`:

        > cd my_files

1. You can list the files in the _remote_ directory with `ls`. You can list the files on the local
   directory with `!dir`.

1. You can print the current _remote_ directory with `pwd` (print working directory), and the _local_ remote directory
   with `lpwd`.

1. To copy a file from the _local_ folder to the _remote_ folder, you use `put`. e.g.

        > put local_file.txt

   See `help put` to show more options. It supports recursive copies.

1. To copy a file from the _remote_ folder to the _local_ folder, you use `get` e.g.:

        > get remote_file.txt

   See `help get` to show more options. It supports recursive copies.

1. If you need to create new directories you can use the `mkdir` (on the remote, and
   `!mkdir` on the local side.

1. If you intend to copy more files later on, you can leave the connection open.

1. If you're done transferring things, you can issue the `quit` command.
