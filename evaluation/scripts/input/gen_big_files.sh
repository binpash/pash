#!/bin/bash

rm -f 10G.txt
touch 10G.txt
for (( i = 0; i < 10; i++ )); do                                               
  cat 1G.txt >> 10G.txt                                                           
done

rm -f 100G.txt
touch 100G.txt
for (( i = 0; i < 10; i++ )); do                                               
  cat 10G.txt >> 100G.txt                                                           
done

rm -f all_cmds_x1000.txt
touch all_cmds_x1000.txt
for (( i = 0; i < 100; i++ )); do                                               
  cat all_cmds_x10.txt >> all_cmds_x1000.txt                                                           
done

