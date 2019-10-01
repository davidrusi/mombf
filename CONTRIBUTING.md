# `mombf` contributor guide

You can contribute to `mombf` by raising issues about bugs, suggest
improvements to the documentation and vignettes, etc. You can also contribute
by sending code directly to `mombf`, however, before sending a Pull Request
to `mombf` check out the Issue tracker to see
possible contributions. If you would like to add a new feature, please
submit an issue first with your feature request to see if it fits with
what we have in mind.

To submit a Pull Request, fork `mombf`'s GitHub repository, fix some bug, improve the docs...
and submit a Pull Request on GitHub.

# `mombf` developer guide
This section of the Contributing guide is intended for _developers_
of `mombf` in the sense that they are expected to release new versions
of the package, which is done via R-forge. All other contributions can
be done with GitHub forks and Pull Requests.

## Get Started

`mombf` uses both GitHub and R-forge. GitHub is used to track
changes, issues and handle pull requests, and R-forge is used to
compile system specific binaries and interact with CRAN. In order
to create a folder tracking both platforms, start by creating a
repository from R-forge:

    $ git svn clone svn+ssh://user@scm.r-forge.r-project.org/svnroot/mombf/pkg/mombf/

This will create a folder called `mombf` containing the R package. If `git
svn` does not work, try installing `git-svn`, in Debian/Ubuntu this can be
done with `sudo apt install git-svn`. Check that the output of `ls mombf` is
similar to:

    ChangeLog  CONTRIBUTING.md  data  DESCRIPTION  INDEX
    LICENSE  man  NAMESPACE  R  src  vignettes

Enter into the `mombf` folder and check that `.git/config` contains something
similar to:

    [svn-remote "<rforge/svn>"]
        url = path/or/url/to/svn/repository
        fetch = :refs/remotes/<r-forge/git>-svn

Rename the master branch:

    $ git branch -m r-forge_local

Set the local repository to track the GitHub repository too, and create a
branch following GitHub's master:

    $ git remote add origin https://github.com/davidrusi/mombf.git
    $ git checkout master

Finally, checkout and create a feature branch to work on. It is not
recommended to work on master.

## Merge GitHub changes back into R-forge
Once there are enough changes to release a new version of `mombf`, changes in
GitHub must be merged into R-forge:

    $ git checkout r-forge_local
    $ git merge master
    $ git rebase

And finally, push this changes to R-forge:

    $ git svn rebase
    $ git svn dcommit
