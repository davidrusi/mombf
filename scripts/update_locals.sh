#!/bin/bash

set -e # fail on first error

PULL=false
MERGE=false
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -p|--pull)
        PULL_BRANCH="$2"
        PULL=true
        shift # past argument
        shift # past value
        ;;
        -m|--merge)
        MERGE_BRANCH="$2"
        MERGE=true
        shift # past argument
        shift # past value
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# store current branch
CURRENT_BRANCH="$(git branch | grep \* | cut -d ' ' -f2)"

# Update r-forge_local
printf "Updating r-forge_local with changes from R-forge repo\n"
git checkout r-forge_local
git svn rebase

# Update Github and rebase so that both trees match
printf "\nUpdating master with changes from GitHub\n"
git checkout master
git fetch origin
git rebase origin/master


if [ "$PULL" = true ] ; then
    printf "\nUpdating ${PULL_BRANCH} with changes from GitHub\n"
    git checkout "$PULL_BRANCH"
    git pull
fi

if [ "$MERGE" = true ] ; then
    printf "\nMerging master changes into ${PULL_BRANCH}\n"
    git checkout "$MERGE_BRANCH"
    git merge master
fi

# return to original branch
git checkout "$CURRENT_BRANCH"
