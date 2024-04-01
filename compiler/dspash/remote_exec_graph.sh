ir_file=$1

worker_manager_client="$DISH_TOP/pash/compiler/dspash/worker_manager_client.py"
python3 "$worker_manager_client" "$ir_file" "$declared_functions" "$DSPASH_SOCKET" "$DSPASH_REEXEC_SOCKET"
