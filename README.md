# sdsh
Scalable Distributed Shell or something!

The key insight here is that the key shell abstractions (streams and pipes) are trivially distributed in a scalable and fault-tolerant manner---for example, similar to Spark's RDD's.

# Possible Homework Breakdown

* Create a few categories---trivially distributable (e.g., pure functions), might require some coordination, 

* Break down most primitives 

  * maybe some shells (such as Plan 9's `rc` shell) have primitives that are naturally more amenable to distribution than others.

* What about variable sharing (mostly write once / read many times)

* Other constructs such as `if` and `for`?

* Need to Redirect input / output streams

* We need to look into an extensibility. What if the user, in their environment, have access to primitives not "known" to the shell?

# A Few Shell Examples

Examples moved to [scripts](./scripts).

# (Some) Related Work

In shell scripting:
* Attempts on distributed shells: `rc`, 
* Parallel commands: GNU's `parallel`, 
* Tons of parallel scripting: Lsf, Swift, ..

Outside shell scripting:
* [pipelines vs. CSP etc.](https://swtch.com/~rsc/thread/)
* [Pipeline minutiae](https://en.wikipedia.org/wiki/Pipeline_(Unix))
* [RDDs](https://www.usenix.org/system/files/conference/nsdi12/nsdi12-final138.pdf)
* [Smoosh papers](http://shell.cs.pomona.edu/)
* [dgsh](https://github.com/dspinellis/dgsh)

# Evaluation

There are several directions for evaluating these ideas:

* Find  "scripts" out in the  wild and "auto-convert" them  to parallel. Compare
performance (and, possibly, correctness) over sequential runs.

*  Take distributed  / MR  / Spark  / Stream  programs that  map to  our subset;
convert  them  to  sequential  shell   scripts  (manually),  and  then  back  to
distributed scripts (w/ our tool), and compare.

# Notes

## K: On distributing user written commands

It seems that a very important feature of a distributed shell, would
be to correctly distribute any command that might have been written by
the user. By default we should never distribute a command that hasn't
been analysed to make sure that it can be distributed safely.

It might be feasible to study some of the standard shell commands
(grep, ls, ...) and implement a distributed version of them (either
custom or on top of the standard implementation). However what happens
with user-written commands? In order to deal with this, we should come
up with some abstraction, that represents how to 'lightly' distribute
each user written command. This could be given to the system by the
user using some configuration file, so that our implementations knows
how and when to distribute a user-written command. 

Is this a necessary feature for a distributed shell? If so, what is a
nice abstraction that represents how "much" can one command be lightly
distributed. Has something like that been implemented in any other
distributed shell?

By lightly distributed, I mean that the command is run several times
as is, the only extension being that the input is split and the output
is merged appropriately.

# Future Ideas

> Look into literature  on streaming sorts (log steps); the  first few things do
> not have  to wait  for the whole  thing---FPGA merge-sorts and  how do  you go
> through the various stages
> (Andre')

Since this is a relatively simplified and well-understood model, we could easily
take into account scheduling aspects---i.e.,  where to schedule pipeline stages,
based on available hardware. For example, I might have GPUs / FPGAs, or a set of
tightly-coupled processors.

