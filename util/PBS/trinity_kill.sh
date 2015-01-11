#!/bin/bash

# much slower than perl version
# read jobnumbers in current directory (and directory passed as variable) and kill them (start at the last job and move to oldest)

if [ $1 ]; then
	if [ -e "$1/jobnumbers.out" ];then
		tac "$1/jobnumbers.out" | 
		while read line
		do
			echo "$1/jobnumbers.out": Stopping $line
			qdel -W force $line
			qdel -W force $line >/dev/null 2>/dev/null
		done
		rm -f "$1/jobnumbers.out"
	fi
	exit 0
fi

if [ -e jobnumbers.out ];then
		tac jobnumbers.out | 
		while read line
		do
			echo jobnumbers.out: Stopping $line
			qdel -W force $line
			qdel -W force $line >/dev/null 2>/dev/null
		done
		rm -f jobnumbers.out
else
	echo "No previous jobs found"
fi
