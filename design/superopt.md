
The Unix shell has a widely diverse set of uses and thus PaSh's parallelization attempt has to be aware of several constraints.
For example, short-running one-liners, possibly typed by the user during an interactive session of the shell, should not see any slowdown.
The same script may be run on vastly different environments ranging from low-capability devices to multiprocessor servers used for compute-intensive workloads.

0. The basic model we have discussed takes into account the
parallelism factor: run PaSh with different widths on test inputs to
identify the best configurations.

1. On top of that, if (1) the input is not large enough, or if (2) the
pipeline/program does not execute for long enough, the superoptimizer
should say that it does not even make sense to run the compiler —
possibly not even the full superoptimizer. This is solvable by first
running the sequential version on, say, [the input, half the input,
quarter the input] and trying to extrapolate how long the full
sequential run would take.

2. If speedup from several stages is not good enough, it might make
sense to keep them sequential. For example, if a sort is followed by
`tr`, `tr`, it might make sense to not put a `split` after collecting
the results of `sort`.

3. For specific hardware configurations (say, 2x CPUs), it might not
even make sense to parallelize a program. i.e., there are programs (or
individual commands) that will give speedup above a certain
parallelism factor (lower bound) and only below a different one (upper
bound).

4. We could assume the developer gives the superoptimizer a global
"optimization budget" in terms of hours: for example what can we do
(optimize) in, say, a minute, an hour or overnight? Some optimization
checks are easier than others, and some are more expensive than
others—so we should prioritize the optimization axes and start with
ones that are cost-efficient (this could lead to a nice algorithm).

5. We could also assume we can spend some time (5.1) checking commands
for side-effects (via some form of system tracing), and (5.2)
synthesizing some of the parallel aggregators for commands in the
program that fall outside the GNU Coreutils  and POSIX subsets. Both
are highly parallelizable tasks—so if we assume a fixed budget (see
4), how many resources can we spend on (5.1) and (5.2)?

6. Later we might want to model the fixed runtime costs of a
distributed infrastructure — i.e., edges in the dataflow graph may
also model network channels (rather than pipes), when nodes in the
graph run on different physical computers. This would inform the
viability of different parallelization configurations (like --width),
and even allow some form of fusion for the distributed case (where we
co-locate multiple stages).

Generally, I have focused here on making the "fast-path" faster, where
if we're not improving things by a lot we should immediately abort and
run the sequential code.
This will look great in the evaluation, because we can show the system
never slows down anything, but also that when it improves things, it
improves them by a lot.
