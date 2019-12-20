#!/usr/bin/env bash


if [ "$1" == "grid" ]; then
   shift
	 echo $@
   echo "GRID SUBMISSION"
	 $@ 
fi

if [ "$1" == "shell" ]; then
   shift
	 echo $@
	 last_id=$#
   last_element=${@:last_id}
   echo "SHELL SUBMISSION " $last_element
   sh $last_element &
fi
