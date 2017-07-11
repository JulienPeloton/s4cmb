#!/bin/bash

## Script to measure the coverage of the test suite (via doctest).
## Launch it using ./coverage
## and open the html files under the folder htmlcov/
for i in instrument/*.py
do  
 echo 'coverage for '$i 
 coverage run $i
 coverage report $i 
 coverage html $i
done

for i in time_ordered_data/*.py
do
 echo 'coverage for '$i
 coverage run $i
 coverage report $i
 coverage html $i
done
