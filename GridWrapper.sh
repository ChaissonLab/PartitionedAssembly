#!/usr/bin/env bash


if [ "$1" == "grid" ]; then
   shift
	 echo $@
	 $@
fi

if [ "$1" == "shell" ]; then
   shift
	 echo $@
	 last_id=$#
   last_element=${@:last_id}
   bash $last_element
fi
