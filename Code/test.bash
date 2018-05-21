#!/bin/bash
if [ -z $1 ]
then
	echo "Usage : $0 nomDesImages"
fi
make clean
make mrt
mkdir img/$1
./mrt img/$1/$1\0 0
./mrt img/$1/$1\1 1
./mrt img/$1/$1\2 2
./mrt img/$1/$1\3 3
exit 0
