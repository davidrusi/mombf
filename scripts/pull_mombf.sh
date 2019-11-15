#!/bin/bash

set -e # fail on first error

# store current branch
CURRENT_BRANCH="$(git branch | grep \* | cut -d ' ' -f2)"

PULL=false
MERGE=false
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -p|--pull)
        if [[ -n "$2" ]]; then
          PULL_BRANCH="$2"
          shift # past value
        else
          PULL_BRANCH="$CURRENT_BRANCH"
        fi
        shift # past argument
        PULL=true
        ;;
        -m|--merge)
        if [[ -n "$2" ]]; then
          MERGE_BRANCH="$2"
          shift # past value
        else
          MERGE_BRANCH="$CURRENT_BRANCH"
        fi
        shift # past argument
        MERGE=true
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

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
    printf "\nUpdating %s with changes from GitHub\n" "${PULL_BRANCH}"
    git checkout "$PULL_BRANCH"
    git pull
fi

if [ "$MERGE" = true ] ; then
    printf "\nMerging master changes into %s\n" "${PULL_BRANCH}"
    git checkout "$MERGE_BRANCH"
    git merge master
fi

# return to original branch
printf "\nReturning to original branch %s\n" "${CURRENT_BRANCH}"
git checkout "$CURRENT_BRANCH"
