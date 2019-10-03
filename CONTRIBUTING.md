`mombf` contributor guide
-------------------------

You can contribute to `mombf` by raising issues about bugs, suggest
improvements to the documentation and vignettes, etc. You can also contribute
by sending code directly to `mombf`, however, before sending a Pull Request
to `mombf` check out the Issue tracker to see
possible contributions. If you would like to add a new feature, please
submit an issue first with your feature request to see if it fits with
what we have in mind.

To submit a Pull Request, fork `mombf`'s GitHub repository, fix some bug, improve the docs...
and submit a Pull Request on GitHub.

`mombf` developer guide
-----------------------

This section of the Contributing guide is intended for _developers_
of `mombf` in the sense that they are expected to release new versions
of the package, which is done via R-forge. All other contributions can
be done with GitHub forks and Pull Requests.

## Background
`mombf` uses both GitHub and R-forge. GitHub is used to track
changes, issues and handle pull requests, and R-forge is used to
compile system specific binaries and interact with CRAN. This means that
there are two _different_ version control systems tracking the _same code_.
Thus, **it is vital to check the branch you are working on is up to date
before working**.

The general idea behind this is the following. R-forge uses svn which is
centralized, therefore, in order to avoid problems with its log, we will
use GitHub to develop in a distributed way _starting from svn HEAD's_. The
steps in the development process can then be summarized as follows. New
features are developed on GitHub, and merged into its master branch (to
simplify everything, it is recommended to squash merge). Right after merging
a change into GitHub master, the change is updated to R-forge (which will give
the new commit a different identifier) and GitHub master branch is eventually
rebased to R-forge master.

If the changes are developed in R-forge, then before merging anything into
GitHub master, it should be checked that GitHub master is up to date with
R-forge.

## Get Started

In order
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
    $ git fetch origin
    $ git checkout master

Finally, checkout and create a feature branch to work on. It is not
recommended to work on master.

## Merge GitHub changes back into R-forge
Once there are enough changes to release a new version of `mombf`, changes in
GitHub must be merged into R-forge. After making sure that local master is up
to date with GitHub master:

    $ git checkout r-forge_local
    $ git merge master

Afterwards, push this changes to R-forge:

    $ git svn rebase
    $ git svn dcommit

Finally, rebase GitHub master so that it shares the same tree as R-forge
repository.

    $ git checkout master
    $ git rebase r-forge_local
    $ git push -f

Note that the push of the rebased master to GitHub must be a forced one
because we are actually updating the commit ids to match the ones in R-forge.
However, it should also be noted that the code is **not** changed, so there
should never be merge problems while rebasing.
