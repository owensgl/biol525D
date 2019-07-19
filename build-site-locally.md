
## instructions to run the jekyll stuff from your local machine

Pushing changes to github will force github pages to regenerate pages of the website.
If you want a quicker turnaround time for visually inspecting changes, you have to build
the site locally.

 1. Once: install RVM

      # based on https://rvm.io/rvm/install .
            gpg --keyserver hkp://pool.sks-keyservers.net --recv-keys 409B6B1796C275462A1703113804BB82D39DC0E3 7D2BAF1CF37B13E2069D6956105BD0E739499BDB
	          curl -sSL https://get.rvm.io | bash -s stable --ruby

 1. Once: install Ruby (> 2.1.0)

      source ~/.rvm/scripts/rvm

      # this ruby version matches the one in .ruby-version
            rvm install ruby-2.4.5

      # set as default (optional)
            rvm --default use 2.4.5

 1. Once: install additional ruby gems in your custom ruby install

      cd <repodir> && rvm use
      gem install bundler
      bundle install # (installs stuff from ./Gemfile)

 1. Each time: When working on the site, use the same version of ruby with the project's Gemset

      source ~/.rvm/scripts/rvm
      cd <repodir> && rvm use

 1. The rest, building, deploying, is jekyll specific
     See https://help.github.com/en/articles/setting-up-your-github-pages-site-locally-with-jekyll

    View local site: `bundle exec jekyll serve`

    (this makes the site available for you to browse at localhost:4000/, and will regenerate the site
     when files on your computer change)

