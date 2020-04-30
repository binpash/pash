#!/bin/bash

func() { read -r distro setup; echo $distro $setup; }

export -f func

cat ../evaluation/usecases/shellcheck/temp_input.txt |
    func
