##
## A library of shell functions that can be used to 
## easily create a building/dependency installing/input 
## downloading scripts.
##



##
## This function checks if all the files in the arguments exist
## It returns 0 if all files exist, or 1 otherwise
##
files_exist_done_check()
{
    for file in "$@"; do
        if [ ! -f "$file" ]; then
            return 1
        fi
    done
    return 0
}

##
## This function checks if number of files in a sequence of directories
## is correct.
## Returns 0 if number is correct, or 1 otherwise
##
number_of_files_in_dir()
{
    local expected_number=$1
    local actual_number=$(ls "${@:2}" | wc -l)
    if [ $expected_number -eq $actual_number ]; then
        return 0
    else
        return 1
    fi
}

##
## This function executes a single idempotent step only if its check fails
##
## Requirements:
## - The step needs to be idempotent
## - The check needs to also check file sizes if there is concern of non-idempotence or failed download
##
execute_step()
{
    local step_fun=$1
    local step_done_check_fun=$2
    local step_desc=${3:-"Execution step"}

    # shellcheck disable=SC2086
    if ! eval $step_done_check_fun; then
        echo "$step_desc is not done, executing..."
        # shellcheck disable=SC2086
        eval $step_fun
        # shellcheck disable=SC2086
        eval $step_done_check_fun || { echo "ERROR: $step_desc failed!"; exit 1; }
    fi
    echo "$step_desc completed."
}

## Issues:
##
## - An overarching problem is that these take time in general, 
##   and therefore testing them out is not really feasible.
## - Another problem is that by doing that manually, 
##   we cannot get completely fine-grained. For example, we could
##   only copy the missing file _a la_ Rattle, instead of running
##   the whole step.
## - Another problem is that idempotence checking is hard to do manually.
## - Another issue is that generating the checks is cumbersome and error-prone.
##   Users need to think whether they need file_exists/number_of_files/size checks,
##   and if they are downloading, they need to first download and then determine the check.
## 