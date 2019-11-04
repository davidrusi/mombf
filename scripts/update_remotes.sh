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
git checkout r-forge_local
git merge master

# Update R-forge remote
git svn rebase
git svn dcommit

# Update Github and rebase so that both trees match
git checkout master
git rebase r-forge_local
if [ "$FORCE" = true ] ; then
    git push -f
else
    git push
fi

# return to original branch
git checkout "$CURRENT_BRANCH"
