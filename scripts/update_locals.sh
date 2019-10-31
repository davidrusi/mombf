#!/bin/bash

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

# Update r-forge_local
git checkout r-forge_local
git svn rebase

# Update Github and rebase so that both trees match
git checkout master
git fetch origin
git rebase origin/master


if [ "$PULL" = true ] ; then
    git checkout "$PULL_BRANCH"
    git pull
fi

if [ "$MERGE" = true ] ; then
    git checkout "$MERGE_BRANCH"
    git merge master
fi
