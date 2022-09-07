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

install_deps_source_setup() {
    # move to the input directory
    cd input/
    # check if there are dependencies
    if [ -e install-deps.sh ]; then
        echo "Installing dependencies"
        bash install-deps.sh
    fi
    # source the setup file
    # it contains the fetch dataset function
    # and the export variable function for IN, IN_PRE
    source setup.sh
    # fetch the dataset
    setup_dataset $1 > /dev/null
    cd ..
}
#########################
# The command line help #
#########################
usage() {
    echo "Usage: `basename $0` [option...] -- shell script to build PaSh"
    echo
    echo "   -h, --help                 Show this help message"
    echo "   -o, --opt-agg              Install g++-10 and switch to it as main compiler. Build the optimized c++ aggregators (run with sudo)"
    echo "   -s, --show-deps            Show all the required dependencies (does not setup/deploy PaSh nor its dependencies)"
    echo "   -e, --install-eval         Install all the dependencies needed for reproducing the evaluation figures (uses sudo, only for Ubuntu/Debian currently)"
    echo
    exit 1
}

##########################################
# Install all the required libraries and #
# dependencies for PaSh evaluation       #
##########################################
install_eval_deps() {
    echo "Installing evaluation dependencies (needs sudo)"
    # needed for majority of the benchmarks (not available in docker instances)
    sudo apt-get install unzip
    paths="$(find $PASH_TOP/evaluation/benchmarks -name install-deps.sh)"
    for f in $(echo $paths); do
        path=$(dirname $(readlink -f $f))
        cd $path
        bash install-deps.sh
        cd - > /dev/null
    done
    echo "Generating PDF plots of the evaluation results is optional and requires R-packages"
    echo "Follow Installation Guide from: $PASH_TOP/evaluation/eval_script/README.md"
}

##########################################
# parse and read the command line args   #
##########################################
read_cmd_args() {
    # Transform long options to short ones
    for arg in "$@"; do
        shift
        case "$arg" in
            "--opt-agg")        set -- "$@" "-o" ;;
            "--show-deps")      set -- "$@" "-s" ;;
            "--install-eval")   set -- "$@" "-e" ;;
            "--help")           set -- "$@" "-h" ;;
            *)                  set -- "$@" "$arg"
        esac
    done

    while getopts 'opsreh' opt; do
        case $opt in
            # passthrough the variable to the Makefile for libdash
            o) export optimized_agg_flag=1 ;;
            s) export show_deps=1 ;;
            r) export show_eval_deps=1 ;;
            e) export install_eval=1 ;;
            h) usage >&2 ;;
            *) echo 'Error in command line parsing' >&2
                exit 1
        esac
    done
}


confirm_installation_works() {
# If installation  hasn't worked, it is difficult to provide a meaningful
# action---apart from checking the logs, which are currently not comprehensive
  echo "Confirming installation works.."
  set +e
  $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh
  if [ $? -ne 0 ]; then
    echo "Something failed, please check logs"
  fi
  set -e
}


append_pash_to_rc() {
  # FIXME: The right files are not entirely clear yet
  rc_configs=(~/.shrc ~/.bashrc  ~/.zshrc ~/.cshrc ~/.kshrc) # add more shell configs here
  for config in "${rc_configs[@]}"
  do
      ## if the config exists
      ## check if it contains an old entry of Pash
      if [ -e "$config" ]; then
          # get the shell name
          shell_name=$(echo $(basename $config) | sed 's/rc//g' | sed 's/\.//g')
          echo "Do you want to append \$PASH_TOP to $shell_name ($config) (y/n)?"
          read answer
          if [ "$answer" != "${answer#[Yy]}" ] ;then 
              tmpfile=$(mktemp -u /tmp/tmp.XXXXXX)
              # create a backup of the shell config
              cp $config ${config}.backup
              # remove all the entries pointing to PASH_TOP and PATH
              grep -ve "export PASH_TOP" $config > $tmpfile
              mv $tmpfile $config
              path_ans=0
              # check if PATH contains PASH_TOP reference
              # we need to store it in a variable otherwise is messes up with the
              # existing environment
              var=$(grep -e "export PATH" $config | grep -e '$PASH_TOP') || path_ans=$?
              # if the return code is 0 -> there is a reference of $PASH_TOP in 
              # PATH, remove it
              if [ "$path_ans" == 0 ]; then
                  # remove previous references to PASH_TOP from PATH
                  grep -v 'export PATH=$PATH:$PASH_TOP' $config > $tmpfile
                  mv $tmpfile $config
              fi
              ## there isn't a previous Pash installation, append the configuration
              echo "export PASH_TOP="$PASH_TOP >> $config
              echo 'export PATH=$PATH:$PASH_TOP' >> $config
          fi
      fi
  done
}

isDocker(){
    local cgroup=/proc/1/cgroup
    test -f $cgroup && [[ "$(<$cgroup)" = *:cpuset:/docker/* ]]
}

isDockerBuildkit(){
    local cgroup=/proc/1/cgroup
    test -f $cgroup && [[ "$(<$cgroup)" = *:cpuset:/docker/buildkit/* ]]
}

isDockerContainer(){
    [ -e /.dockerenv ]
}
