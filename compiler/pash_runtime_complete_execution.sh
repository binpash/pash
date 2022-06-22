#!/bin/bash

##
## Completes execution by measuring and logging execution times and restoring state
##

##
## (6)
##

pash_exec_time_end=$(date +"%s%N")

## TODO: Maybe remove the temp file after execution

## We want the execution time in milliseconds
if [ "$pash_output_time_flag" -eq 1 ]; then
    pash_exec_time_ms=$(echo "scale = 3; ($pash_exec_time_end-$pash_exec_time_start)/1000000" | bc)
    pash_redir_output echo "Execution time: $pash_exec_time_ms  ms"
fi

## Source back the output variables of the compiled script. 
## In all cases we should have executed a script
pash_redir_output echo "$$: (7) Recovering BaSh variables from: $pash_output_var_file"
source "$RUNTIME_DIR/pash_source_declare_vars.sh" "$pash_output_var_file"

## Save the previous `set` state to a variable
pash_redir_output echo "$$: (7) Reading current BaSh set state from: ${pash_output_set_file}"

pash_redir_output echo "$$: (7) Current BaSh set state: $(cat "$pash_output_set_file")"
## WARNING: This has to happen after sourcing the variables so that it overwrites it
pash_previous_set_status=$(cat "$pash_output_set_file")

export pash_input_args
pash_redir_output echo "$$: (7) Arguments (might) have been updated to be: $pash_input_args"

## Propagate the `set` state after running the script to the outer script
## TODO: Maybe move this to the end to avoid spurious failures
pash_redir_output echo "$$: (7) Current PaSh set state: $-"
source "$RUNTIME_DIR/pash_set_from_to.sh" "$-" "$(cat "$pash_output_set_file")"
pash_redir_output echo "$$: (7) Reverted to BaSh set state before exiting: $-"

pash_redir_output echo "$$: (7) Reverting last BaSh exit code: $pash_runtime_final_status"
(exit "$pash_runtime_final_status")
