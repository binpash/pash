#!/bin/bash
# This script is part of a study on OASA's Telematics
# Diomidis Spinellis and Eleftheria Tsaliki
# https://insidestory.gr/article/noymera-leoforeia-athinas

# Hours each bus is active each day

# Records are day, hour, line, bus
<input.csv sed 's/T\(..\):..:../,\1/' | awk -F, '
!seen[$1 $2 $4] { seen[$1 $2 $4] = 1; hours[$1 $4]++; bus[$4] = 1; day[$1] = 1; }
END {
   PROCINFO["sorted_in"] = "@ind_str_asc"
   for (d in day)
     printf("\t%s", d);
   printf("\n");
   for (b in bus) {
     printf("%s", b);
     for (d in day)
       printf("\t%s", hours[d b]);
     printf("\n");
   }
}' >5a.txt
