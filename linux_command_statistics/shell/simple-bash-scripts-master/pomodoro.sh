#!/bin/bash

# A simple shell script to use as a pomodoro app.
# The first argument is the focus time length.
# The second argument is the break length.
# Made by Kiailandi (https://github.com/kiailandi)

wseconds=${1:-25}*60
pseconds=${2:-wseconds/300}*60

while true; do
	date1=$(($(date +%s) + $wseconds))
	while [ "$date1" -ge $(date +%s) ]; do
		echo -ne "$(date -u --date @$(($date1 - $(date +%s))) +%H:%M:%S)\r"
	done
	notify-send "Break" "Time to walk and rest"
	read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
	date2=$(($(date +%s) + $pseconds))
	while [ "$date2" -ge $(date +%s) ]; do
		echo -ne "$(date -u --date @$(($date2 - $(date +%s))) +%H:%M:%S)\r"
	done
	notify-send "Work" "Time to get back to work"
	read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
done
