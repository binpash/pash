#!/bin/bash

##
## Completes execution by measuring and logging execution times and restoring state
##

if [ "$PASH_DEBUG_LEVEL" -ne 0 ] && [ "$pash_avoid_pash_runtime_completion_flag" -ne 1 ]; then
    ##
    ## (5)
    ##
    pash_redir_output echo "$$: (5) BaSh script exited with ec: $pash_runtime_final_status"
    pash_redir_output echo "$$: (5) Current BaSh shell: $pash_previous_set_status"
    pash_redir_output echo "$$: (5) Reverted to PaSh set state to: $-"

    ## Prepare a file for the output shell variables to be saved in
    # pash_output_var_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")
    # # pash_redir_output echo "$$: Output vars: $pash_output_var_file"
    # source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_output_var_file"

    ## Prepare a file for the `set` state of the inner shell to be output
    pash_output_set_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")
    pash_redir_output echo "$$: (5) Writing current BaSh set state to: $pash_output_set_file"
    echo "$pash_previous_set_status" > "$pash_output_set_file"

    ##
    ## (6)
    ##

    ## Source back the output variables of the compiled script. 
    ## In all cases we should have executed a script
    ## TODO: Is this actually necessary? The variables are already set due to source
    ## TODO: If it is not necessary then we can get away without even saving the variables to a file or saving the set state to a file
    # pash_redir_output echo "$$: (7) Recovering BaSh variables from: $pash_output_var_file"
    # source "$RUNTIME_DIR/pash_source_declare_vars.sh" "$pash_output_var_file"

    # ## Save the previous `set` state to a variable
    # pash_redir_output echo "$$: (7) Reading current BaSh set state from: ${pash_output_set_file}"

    # pash_redir_output echo "$$: (7) Current BaSh set state: $(cat "$pash_output_set_file")"
    ## WARNING: This has to happen after sourcing the variables so that it overwrites it
    # pash_previous_set_status=$(cat "$pash_output_set_file")

    export pash_input_args
    pash_redir_output echo "$$: (7) Arguments (might) have been updated to be: ${pash_input_args[@]}"

    ## Restore the set state from a file because it has been rewritten by sourcing variables
    export pash_previous_set_status="$(cat "$pash_output_set_file")"
fi