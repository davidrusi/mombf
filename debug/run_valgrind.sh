#!/bin/sh
R -d "valgrind --log-file=valgrind.log --leak-check=full" --vanilla < valgrind.R
