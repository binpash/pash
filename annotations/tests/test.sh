#!/bin/bash

for test_case in cmd-*.sh; do
  ./$test_case || exit
done


