#!/bin/bash

# store current branch
CURRENT_BRANCH="$(git branch | grep \* | cut -d ' ' -f2)"

# make sure repository is up to date
bash scripts/pull_mombf.sh

# merge GitHub changes into r-forge_local
printf "Merging master into r-forge_local\n"
git checkout r-forge_local
git merge master

# Update R-forge remote
printf "Pushing r-forge_local to R-forge repo\n"
git svn rebase
git svn dcommit

# Update Github and rebase so that both trees match
printf "Making master match r-forge_local\n"
git checkout master
git rebase r-forge_local

printf "Pushing master to GitHub repo\n"
git push -f

# return to original branch
printf "\nReturning to original branch ${CURRENT_BRANCH}\n"
git checkout "$CURRENT_BRANCH"
