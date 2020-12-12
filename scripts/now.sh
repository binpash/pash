#!/bin/bash

ls -lt ../../reports/ | head -n 2 | tail -n 1 | cut -d ' ' -f7-10
