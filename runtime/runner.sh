#!/usr/bin/env bash
#
# Runner - Executes preprocessed shell scripts.
#
# This is a standalone executable that runs a shell script with
# appropriate bash flags and arguments.
#

set -e

# Default values
script_path=""
shell_name=""
declare -a script_args=()
allexport_flag="+a"
verbose_flag=""
xtrace_flag=""
debug_level=0

# Usage function
usage() {
    cat << EOF
Usage: runner.sh SCRIPT_PATH SHELL_NAME [SCRIPT_ARGS...] [OPTIONS]

Execute a preprocessed shell script with bash.

Arguments:
  SCRIPT_PATH         Path to the preprocessed script to execute
  SHELL_NAME          Name to use for \$0 in the executed script
  SCRIPT_ARGS         Arguments to pass to the script

Options:
  -a                  Enable the allexport shell option
  +a                  Disable the allexport shell option (default)
  -v                  Print shell input lines as they are read
  -x                  Print commands and their arguments as they execute
  -d, --debug LEVEL   Set debug level (0=none, 1=info, 2=debug)
  -h, --help          Show this help message
EOF
    exit 0
}

# Logging function
log() {
    if [ "$debug_level" -ge 1 ]; then
        echo "PaSh: $*" >&2
    fi
}

# Parse arguments
if [ $# -lt 2 ]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Get positional arguments first
script_path="$1"
shell_name="$2"
shift 2

# Parse remaining arguments (flags and script args)
while [ $# -gt 0 ]; do
    case "$1" in
        -a)
            allexport_flag="-a"
            shift
            ;;
        +a)
            allexport_flag="+a"
            shift
            ;;
        -v)
            verbose_flag="-v"
            shift
            ;;
        -x)
            xtrace_flag="-x"
            shift
            ;;
        -d|--debug)
            if [ -z "$2" ] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
                echo "Error: --debug requires a numeric argument" >&2
                exit 1
            fi
            debug_level="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            # Everything else is a script argument
            script_args+=("$1")
            shift
            ;;
    esac
done

# Verify script exists
if [ ! -f "$script_path" ]; then
    echo "Error: Script not found: $script_path" >&2
    exit 1
fi

# Log execution details
log "Runner starting..."
log "Script path: $script_path"
log "Shell name: $shell_name"
log "Arguments: ${script_args[*]}"
log "Flags: a=$allexport_flag, v=$verbose_flag, x=$xtrace_flag"
log "----------------------------------------"

# Build bash command
bash_cmd="/usr/bin/env bash"

# Add flags
bash_flags="$allexport_flag $verbose_flag $xtrace_flag"

# Execute the script
# Uses: bash [flags] -c "source script" shell_name [args...]
log "Executing: $bash_cmd $bash_flags -c \"source $script_path\" $shell_name ${script_args[*]}"

exec $bash_cmd $bash_flags -c "source $script_path" "$shell_name" "${script_args[@]}"
