#!/bin/bash

cd $HOME/downloads/
PATH=/bin:/usr/bin

line=1
max_proc=1 #numero de downloads ?
list_file="$HOME/downloads/todo.txt"
#prog="/usr/local/bin/wget"
prog="wget"

while true
do
	while true
	do
		proc=`ps -f -u $USER | grep -c $prog`
		# grep is in the list too
		let proc--

		lines=`grep -c "" $list_file`
		echo "Proc: $proc / $max_proc	Line: $line / $lines"

		[[ $proc -ge $max_proc || $line -gt $lines ]] && break

		params=`grep -n "" $list_file | grep "^$line:" |
		                               sed -e "s/^$line://"`
		echo $params | tee -a done.txt archive.txt

		# ignore empty lines
		if [ "$params" ]; then
			$prog -b $params
			sleep 3
		fi
		let line++
	done

	echo "Waiting..."
	sleep 10
done