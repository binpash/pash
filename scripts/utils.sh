#!/usr/bin/env bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

# another solution for capturing HTTP status code
# https://superuser.com/a/590170

eexit(){
  echo $1 'please email pash-devs@googlegroups.com'
  exit 1
}

nargs(){
  echo $# $1 $2
}

rm-files(){
  echo "${@}"
  rm -r "${@}"
  exit 0
}

append_nl_if_not(){
  ## Adds a newline at the end of a file if it doesn't already end in a newline.
  ## Used to prepare inputs for PaSh.
  if [ -z "$1" ]; then
    echo "No file argument given!"
    exit 1
  else
    if [ ! -f "$1" ]; then
      echo "File $1 doesn't exist!"
      exit 1
    else
      tail -c 1 "$1" | od -ta | grep -q nl
      if [ $? -eq 1 ]; then
        echo >> "$1"
      fi
    fi
  fi
}


# TODO: This logic appears verbatim in a few places.
# Deduplicate by having other scripts call this one.
infer_unix_like_distro() {
    if type lsb_release >/dev/null 2>&1 ; then
	distro=$(lsb_release -i -s)
    elif [ -e /etc/os-release ] ; then
	distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
    fi

    printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]'
}


use_alternative_debian_gpp10() {
    apt-get install software-properties-common -y
    add-apt-repository ppa:ubuntu-toolchain-r/test -y
    apt-get install g++-10
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
    update-alternatives --set g++ /usr/bin/g++-10
}
