
## instructions to run the jekyll stuff from your local machine

 1. Once: install RVM
      gpg --keyserver hkp://pool.sks-keyservers.net --recv-keys 409B6B1796C275462A1703113804BB82D39DC0E3 7D2BAF1CF37B13E2069D6956105BD0E739499BDB
      curl -sSL https://get.rvm.io | bash -s stable --ruby

 1. Once: install Ruby (> 2.1.0)

      source ~/.rvm/scripts/rvm

      rvm install ruby-2.4.5
      # set as default
      rvm --default use 2.4.5

 1. Once: install additional ruby gems in your custom ruby install

      gem install bundler
      bundle install # (installs stuff from ./Gemfile)

 1. Each time: When working on the site, use the same version of ruby

      source ~/.rvm/scripts/rvm
      rvm use 2.4.5

 1. The rest, building, deploying, is jekyll specific
    See https://help.github.com/en/articles/setting-up-your-github-pages-site-locally-with-jekyll

    View local site: `bundle exec jekyll serve`
