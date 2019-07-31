Biol525D
=======

We recommend that you browse the course content from the following webpage:

>> https://owensgl.github.io/biol525D/


This repository contains all the files used to generate the website.

## OBTAINING ALL THE FILES ON THIS SITE

To access the course content offline, you may download an up-to-date
snapshot archive of the site content from this location:
https://github.com/owensgl/biol525D/archive/master.zip

But the recommended approach is to use `git` which tracks changes and supports incremental updates:

    git clone https://github.com/owensgl/biol525D.git

If the course content changes, you can update your local copy by going into to the **biol525D** directory that was created by the previous command and invoking the command:

    git pull

## MAKING CHANGES TO THE SITE

Site content can be edited directly from github's builtin
editors. Every time a commit is made to the master branch, github
rebuilds the live version of the site content. Changes can take around
30seconds to take effect.

If you wish to perform a large amount of tweaks, or apply larger
structural changes to the site, you may speed up the iterative test
development cycle, by following instructions in [building the site
locally](./build-site-locally.md) to create a test environment local
to your computer. You would tweak the site to your satisfaction and
push to github to publish your work.

## LICENSE AND COPYRIGHT

Copyright (C) 2015 S. Evan Staton, Sariel Hubner, Sam Yeaman

Modified work (c) 2016, 2017, 2018 Gregory Owens, Kathryn Hodgins

Modified work (c) 2019 Gregory Owens, Kathryn Hodgins, JS Legare

This program is distributed under the MIT (X11) License, which should
be distributed with the package. If not, it can be found here:
http://www.opensource.org/licenses/mit-license.php

The license file is [here](./LICENSE)
