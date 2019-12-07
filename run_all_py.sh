#!/bin/bash

for i in ./*/*/*.py ; do j=`dirname $i` ; k=`basename $i` ; cd $j ; ls $k ; ./$k ; cd - ; done
