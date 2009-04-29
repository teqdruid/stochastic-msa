#!/bin/sh
##
## inflmsa.sh
## 
##  Runs the msa program with "good" options for InflB
##


./msa $1 -k 10 -a 250 -b .05 -m 20 -hg 50 -prog 5 .5 -out $2 -pout $3