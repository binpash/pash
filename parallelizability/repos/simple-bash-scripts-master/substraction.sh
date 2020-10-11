#!/bin/bash
echo .Enter the First Number: .
read a
echo .Enter the Second Number: .
read b
x=$(($a - $b))
echo $a - $b = $x
