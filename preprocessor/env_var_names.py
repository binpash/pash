##
## Variable names used in the pash runtime
##


def loop_iters_var() -> str:
    return "pash_loop_iters"


def loop_iter_var(loop_id: int) -> str:
    return f"pash_loop_{loop_id}_iter"
