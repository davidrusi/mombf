#!/bin/bash

FORCE=false
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -f|--force)
        FORCE=true
        shift # past argument
        shift # past value
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# store current branch
CURRENT_BRANCH="$(git branch | grep \* | cut -d ' ' -f2)"

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
if [ "$FORCE" = true ] ; then
    git push -f
else
    git push
fi

# return to original branch
printf "\nReturning to original branch ${CURRENT_BRANCH}\n"
git checkout "$CURRENT_BRANCH"
