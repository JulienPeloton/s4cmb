#!/bin/bash

## Script to measure the coverage of the test suite (via doctest).
## Launch it using ./coverage
## and open the html files under the folder htmlcov/
## Skip xpure.py as it is not really part of the pipeline
for i in s4cmb/*.py
do
    coverage run -a --source=s4cmb --omit=s4cmb/xpure.py $i
    #coverage report $i
    #coverage html $i
done
