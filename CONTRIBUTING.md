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

## Contents
1. [Background](#Background)
2. [Getting Started](#Getting-Started)
3. [Development on GitHub](#Development-on-GitHub)
   1. [Development on a feature branch](#Development-on-a-feature-branch)
   2. [Development on master](#Development-on-master)
1. [Merge GitHub changes back into R-forge](#Merge-GitHub-changes-back-into-R-forge)
   1. [`push_mombf.sh`](#push_mombfsh)
1. [Syncing the local repository with GitHub and R-forge](#Syncing-the-local-repository-with-GitHub-and-R-forge)
   1. [`pull_mombf.sh`](#pull_mombfsh)
1. [Test suite and code coverage](#Test-suite-and-code-coverage)

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

## Getting Started

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

## Development on GitHub
#### Development on a feature branch
The development on GitHub is based on Pull Requests, which allows to automatically check the package with `R CMD check`, run tests and check code coverage automatically. In order to take full profit of all this features, the procedure should be the following:

1. Make sure `master` is up to date. See [Syncing the local repository with GitHub and R-forge](#Syncing-the-local-repository-with-GitHub-and-R-forge) for details. Generally, it should be running [`pull_mombf.sh`](#pull_mombfsh).
1. Switch to a new branch using `git checkout -b <branch-name> master` and edit and commit there.
1. `git commit` to add the changes to the feature branch.
1. Push the changes to Github with `git push`. You may get an error if this is the first time the branch is pushed to GitHub, but Git itself should provide the correct command, execute Git's recommendation.

Once the changes have been pushed to Github, it is recommended to use the website interface to add any comments on the new feature, check that continuous integration builds are run without errors, that code coverage is acceptable and merge the new feature into the master branch using the button "Squash and merge".

#### Development on master
For simple changes, it may be better to directly commit to `master`. In this case, the workflow would be simpler and faster because no checks are performed before merging and GitHub's web interface is not used.

1. Make sure `master` is up to date. See [Syncing the local repository with GitHub and R-forge](#Syncing-the-local-repository-with-GitHub-and-R-forge) for details. Generally, it should be running [`pull_mombf.sh`](#pull_mombfsh).
2. Work on feature
3. `git commit` to add changes to `master`
4. Follow the steps in [Merge GitHub changes back into R-forge](#Merge-GitHub-changes-back-into-R-forge). Generally, it should be running [`push_mombf.sh`](#push_mombfsh)


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

#### `push_mombf.sh`
For convenience, all these commands can be run at once using `push_mombf.sh`. As an extra precaution, it runs `pull_mombf.sh` automatically before executing any commands.

    $ bash scripts/push_mombf.sh

## Syncing the local repository with GitHub and R-forge
Before starting to work it is crucial to make sure our local branches are up to date with their remote repositories. That is, `f-forge_local` has all the changes in R-forge and `master` has all the changes in GitHub. Update `r-forge_local` with R-forge changes using `git svn rebase` (equivalent to `svn update`):

    $ git checkout r-forge_local
    $ git svn rebase

Update `master` with GitHub changes:

    $ git checkout master
    $ git fetch origin
    $ git rebase origin/master

When updating a local branch we may want to pull the changes from its remote branch or merge the updates in master to the branch. To sync the local branch with its remote, follow the same steps described right above changing `master` by the branch name. To merge the changes in master into the master branch, run the following commands **after** having updated master:

    $ git checkout local_branch
    $ git merge master

#### `pull_mombf.sh`

For convenience, a bash script, `pull_mombf.sh` has been added to convert all this instructions to a single line of code. It will always update `r-forge_local` and `master` when called:

    $ bash scripts/pull_mombf.sh

Moreover, it can also update a local branch with its remote:

    $ bash scripts/pull_mombfs.sh -p local_branch

Or merge master into a local branch:

    $ bash scripts/pull_mombf.sh -m local_branch

The argument `local_branch` can be ommited, which will result in `local_branch` being set to the current branch. For example, if we are on the branch `extend_contributing`, the following command:

    $ bash scripts/pull_mombf.sh -p

Will pull from their original remotes `r-forge_local`, `master` and `extend_contributing`, printing an output similar to:

```
Updating r-forge_local with changes from R-forge repo
...
Current branch r-forge_local is up-to-date.

Updating master with changes from GitHub
...
Current branch master is up-to-date.

Updating extend_contributing with changes from GitHub
...
Already up-to-date.

Returning to original branch extend_contributing
Already on 'extend_contributing'
```

## Test suite and code coverage
Test are run automatically whenever someone sends a Pull Request to mombf to check whether the changes would break any of mombf's functionalities. They are also automatically run on any commit pushed to GitHub's master branch to let users know the code in GitHub is installable and usable just by looking at the badge in the [README.md](https://github.com/davidrusi/mombf#mombf). On every successful test suite execution, a coverage report is directly generated by `codecov`.

Moreover, tests can be run locally to ensure local changes do not break anything. The easiest way is using `devtools::test()`. Running tests locally requires having `testthat` and `patrick` installed. Code coverage can also be checked locally using `covr::report()`, which requires having R package `covr` and `gcov` installed.
