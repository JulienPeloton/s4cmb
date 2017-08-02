#!/bin/bash

## Script to measure the coverage of the test suite (via doctest).
## Launch it using ./coverage
## and open the html files under the folder htmlcov/
## Skip xpure.py as it is not really part of the pipeline
for i in s4cmb/*.py
do
    if [ $i == "s4cmb/xpure.py" ]
    then
        echo 'skip' $i
    else
        coverage run -a --source=s4cmb $i
        #coverage report $i
        #coverage html $i
    fi
done
