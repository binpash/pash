#!/bin/bash

ls -lt ../../reports/ | head -n 2 | tail -n 1 | awk '{print $NF}'
