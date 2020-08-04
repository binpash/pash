OSDI '20 Paper #361 Reviews and Comments
===========================================================================
Paper #361 PaSh: Light-touch Data-Parallel Shell Processing


Review #361A
===========================================================================

Paper summary
-------------
I really enjoyed reading this paper on PaSh and strongly support its acceptance.

The intuition in this paper is that the sequential operations in a script can be converted into a data-flow graph that can be analyzed to extract parallelism, and then transformed into an equivalent script that runs faster without requiring any edits or rewrites from a script developer. PaSh supports command developers via a set of lightweight type annotations that can be used to indicate what operations are parallelizable.

The authors prioritize correctness in their work, and have added comprehensive support for PaSh by supporting POSIX and GNU Coreutils commands via a library. Further, the authors have developed a runtime component that supports common operations in scripts such as data aggregation, splitting, merging, etc. 

The evaluation methodology and comprehensiveness is neat. I was surprised to learn about the 37 challenges shared by Bell Labs as the Unix50 legacy, and that PaSh could execute each pipeline as-is. Further, I found it interesting that PaSh could parallelize execution and improve performance without any developer effort.

Strengths
---------
This paper is a fun read as it introduces an interesting problem, describes an innovative solution, and sprinkles useful nuggets of wisdom throughout the text. For instance, Section 2.1 introduces readers to the UNIX philosophy and how composition allows us to leverage multiple modular tools to solve a complex problem.

* The authors introduce four parallelizability classes that capture all POSIX and GNU Coreutils commands. Further, the authors explicitly note which commands serve as synchronization operations that bound parallelization.

* Section 5 delves into multiple challenges encountered when realizing the design in earlier sections. I especially like how exit calls are handled without causing deadlocks.

* The authors have publicly released their code, benchmarks and evaluation at https://github.com/3487fc/pash

Significant Weaknesses
----------------------
Writing a review for this paper was somewhat challenging as there's no clear comparisons for this work. As such, I note some ways that the paper can be further improved but none of the points below are major critiques of the work:

* I found Section 3.2 and 4.2 difficult to read. After some reflection, I believe that Section 3.2 was difficult to read as the text doesn't quite explain why the authors made decisions around extensibility as they did. For instance, I expect that PaSh maintains a state machine where it tracks which command is being executed and what flags are passed so it appropriately track the "parallelizability" of the operation -- this isn't mentioned in the text.

* I believe that I found Section 4.2 difficult to read as I had to brush up on graph transformations. I'm unsure how the authors can further simplify or explain the points here, but I wished to bring it up. That said, please don't drop this narrative as I think it is useful to include in the paper.

Comments for author
-------------------
Nits:

* Figure 2 does not seem to be referenced from the text

* The conclusion's last sentence seems truncated

Novelty
-------
3. New contribution

Experimental methodology
------------------------
5. Excellent

Writing quality
---------------
4. Well-written

Overall merit
-------------
5. Strong accept (Excellent paper, I will fight strongly for it)

Reviewer expertise
------------------
2. Some familiarity



Review #361B
===========================================================================

Paper summary
-------------
PaSh is an auto-parallelizing POSIX shell. To do so with arbitrary
POSIX commands, PaSh requires a few annotations on these commands that
describe each command's parallelization class out of 4 possible
classes, depending on their command line arguments. PaSh creates a
data-flow graph from parallelizable regions of shell scripts, which it
auto-parallelizes according to command parallelization classes. It
then generates a new, parallel shell script from the parallel
data-flow graph. The authors implement PaSh, provide annotated
versions of GNU Coreutils and POSIX commands, and evaluate PaSh using
a suite of 12 common UNIX one-liners, 37 UNIX pipelines from UNIX's
50th anniversary, NOAA weather analysis, Wikipedia web indexing, and a
bioinformatics shell script. PaSh can parallelize almost all shell
scripts in the suite with speed-up of up to 61.1x at a COST of 2.

Strengths
---------
PaSh is extremely useful. What serious OS researcher wouldn't want it?
A large number of data analysts also use shell scripts that can
benefit from parallelization.  Hence, PaSh is useful beyond the OS
community.

This is a well-written, thorough study of the many gnarly
parallelization problems that lurk in half a century of evolution that
is the POSIX shell environment. The paper takes the study lessons and
applies them to a practical solution.

Significant Weaknesses
----------------------
The complexity of annotating commands is under-evaluated. Since the
usefulness of PaSh hinges on how straight-forward it is to incorporate
commands for parallelization, a more thorough analysis, ideally
quantitatively (complexity of annotations, person-time it took to
develop annotations for Coreutils and POSIX), would be insightful.

Comments for author
-------------------
Kudos to the authors for delivering such a useful tool! The analysis
of parallelization problems with UNIX commands is very informative,
too. I am happy to accept this paper and only have minor comments.

One might debate the novelty of the presented techniques. Using
data-flow graphs and commutativity to determine and optimize parallel
patterns in programs and scripts is not new. Neither is the idea that
data-parallel pipelines may be constructed from shell programs. For
example, Dryad [1] (in particular, the "Nebula" scripting language
described in Section 7.1 of the Dryad paper) was using the same
mechanism and it should get a treatment in the related work section of
this paper. However, I think PaSh's thorough application of these
techniques to the domain of POSIX shell scripts and UNIX commands is
interesting and new.

You did a very thorough analysis of all the factors that affect the
parallelizability of UNIX commands. You found that command line
arguments can change the parallelization behavior of commands and
potentially put them in another class. However, I wonder if the
analysis if complete. For example, could the behavior of commands not
be influenced by the setting of environment variables or in
configuration files? I did not find these mentioned in the
parallelizability study.

References:

[1] Michael Isard, Mihai Budiu, Yuan Yu, Andrew Birrell, and Dennis
Fetterly. 2007. Dryad: distributed data-parallel programs from
sequential building blocks. In Proceedings of the 2nd ACM
SIGOPS/EuroSys European Conference on Computer Systems 2007 (EuroSys
’07). Association for Computing Machinery, New York, NY, USA, 59–72.

Novelty
-------
3. New contribution

Experimental methodology
------------------------
4. Good

Writing quality
---------------
4. Well-written

Overall merit
-------------
4. Accept (Good paper, I will advocate for it)

Reviewer expertise
------------------
3. Knowledgeable



Review #361C
===========================================================================

Paper summary
-------------
The paper presents PaSh, a system for automatically converting a
pipeline of shell scripts into a data-parallel execution plan for
faster completion.  This work demonstrates that the map-reduce data
processing field is sufficiently mature that it can be incorporated
into commodity tools.

Strengths
---------
+ It's an OS paper! (really! it's a paper on improving one of the
  primary ways we interact with the OS!)

+ The paper marks a signficant milestone in the maturity of our
collective understanding of data parallel processing.  

+ The paper demonstrates, somewhat surprisingly, the potential to
treat shell scripts as a processing framework on par with MapReduce,
Spark, etc.

Significant Weaknesses
----------------------
- the paper is "only" an application of data flow processing
  techniques to the shell script environment.

- PaSh cannot operate purely off the shelf, individual commands must
  be properly annotated in order for the system to safely parallelize
  them.

Comments for author
-------------------
This is an interesting paper which simultaneously strikes me for the
audacity of being "simply" an improved shell script and for
encapsulating the effort to bring parallel data processing to the
(shell-script knowing) masses.

My first response to this paper was "it's a shell script, no way".  My
second response was "they've applied DFG to shell scripting, so what?".

The longer I sat with the paper, the more I realized that (a) I never
would have thought to do that, (b) I never would have expected it to
work, and (c) shell scripts probably have the highest surface
area/lowest bar for adoption of	all of the data-parallel frameworks.

The paper is well written and clear.  Technically, it thoroughly
describes the conceptual steps needed to map posix shell scripts into
a DFG.  The clarity of presentation is such that none of the technical
steps stand out as surprising, and at each point it is very easy to
see how (and why) the specific choice was made. And you've managed to
do this for shell scripts!  The one major thing that is missing for me
is the secret sauce.  Where is the core technical novelty?  Or does it
lie simply in the realization that you understand DFG well enough to
be able to apply it to the shell environment?

I would like to see a bit more treatment of how the DFG relates to
that of large distributed parallel processing frameworks (e.g., the
MapReduce lineage).

One piece notably missing from the presentation is the overall
complexity of the annotations required to make PaSh work.  The
language (as shown in the appendix) is small, but it is not clear how
difficult it is to annotate commands correclty, or the performance
impact of being overly conservative in writing the annotations.

nit:  the final sentence of the paper cuts off in the middle.

Novelty
-------
2. Incremental improvement

Experimental methodology
------------------------
4. Good

Writing quality
---------------
4. Well-written

Overall merit
-------------
5. Strong accept (Excellent paper, I will fight strongly for it)

Reviewer expertise
------------------
2. Some familiarity



Review #361D
===========================================================================

Paper summary
-------------
This paper describes a system for annotating shell commands with their parallelism properties and a compiler that then uses these properties generates DAG based parallelism, exploiting those properties, to improve performance through parallelism of shell scripts.  It demonstrates improvements on parallel hardware on large inputs.

Strengths
---------
* Practical system that would be a good addition to standard shell distributions

Significant Weaknesses
----------------------
* No new research results or insights on scripts
* Missing critical details on what the compiler actually does to the code

Comments for author
-------------------
> All (most?) of these scripts were first developed in 1970-1980s or thereabouts, and apparently have not been updated to reflect the reality of parallel hardware platforms on which they run.  Thus, even though hardware parallelism has been ubiquitously available since SMT was introduced in IBM and Intel platforms over 20 years ago, apparently shell commands remain sequential. One choice is for developers to introduce parallelism into each shell command individually, and then configure/tune them based on the platform when they are installed.  This more traditional approach would produce highly tuned commands, but would lack task parallelism between scripts.

kk: The point is not that the commands would be highly tuned but lack
task parallelism. If one could highly tune all commands to be parallel
then they can achieve great performance. Better than Pash. The space
that we are targeting is when this is not feasible to rewrite an old
command to be parallel.


> Thus the potential advantage of this paper is to show that scripts are complex and there is lots of task parallelism that can be exploited between shell commands, unfortunately this paper does not deliver insight or results on this key advantage of their approach.  Unfortunately Figure 8, does not show any results at all on what the compiler actually does to the code. In compiler papers, it is standard to show which optimizations are actually achieving the results.

kk: This might a presentation fault, but the compiler actually
transforms scripts as described in section 4, 5. It performs the
transformations and then adds the primitives to compose everything
together efficiently.  We report different results for different
optimizations enabled (split, no split, no eager) to show the benefit
from each.

> The main results in Figure 7 show that each application requires a slight different approach to achieve the best speed ups, e.g., sort uses "parallel" and spell uses "parallel + B spilt", which shows that to achieve the parallelism potential requires customizing the compiler to each command.  Furthermore, even though the input sizes are very large, e.g., 1 GB, most commands don't scale, and thus the compiler will need to learn that too and stop at the right point.

kk: This might be a presentation fault again. parallel + automatic
split is the best for all. The reason why we don't show it for
e.g. grep, is that it has the exact same performance as the parallel
(as it doesn't add any splits).

> The paper promises insights by studying the Unix50 scripts, but there are no tables or results which explore the nature of these scripts, nor any insight about them. 

kk: Not sure how to respond to that

Many script inputs are small, in fact, the Unix50 scripts were in many cases made artificially bigger (10 GB) to be able to show speedups. The only thing we learn from section 6.2 is that parallelism in some cases improves the performance of these scripts.

The general approach of annotating methods with their parallelism properties and using compiler transformations to create DAG parallelism is well known, although the application of this approach to shell commands is a new application.

Compiler Optimizations for Scalable Parallel systems, Pande and Agrawal (Ed), 2001,   See Section V, Ch. 18. 
Automatic generation of DAG parallelism, Cytron et al. 1989.
A general approach to mapping of parallel computation upon multiprocessor architectures, Kim and Browne, 1988.

kk: We have to read these papers.

Novelty
-------
1. Published before or openly commercialized

Experimental methodology
------------------------
2. Poor

Writing quality
---------------
3. Adequate

Overall merit
-------------
1. Reject (Serious problems, I'll argue against this paper)

Reviewer expertise
------------------
4. Expert



Review #361E
===========================================================================

Paper summary
-------------
This paper introduces PaSh, which is a compiler for shell scripts that analyzes their dependency structure and transforms them into semantically equivalent scripts that execute with a higher degree of parallelism. To achieve this, PaSh introduces an annotation scheme for builtin shell commands and utilities, and assigns each command (or command a particular set of command-line arguments) a "parallelizability class" that indicates how it may be rewritten into a data parallel form (e.g. distributing a stateless function across splits of an input file, or inserting map-reduce style parallelism for more complicated operations like `wc`). The compiler transforms the AST of an original shell script into the parallel form based on identifying groups of parallelizable operators, and additionally inserts partitioning and buffering operators to improve the effective parallelism of execution. 

PaSh is evaluated on several sets of shell scripts downloaded from the web, achieving a range of speedups up to ~60$\times$ for a pure data-parallel job on a 64-CPU server, and up to 15$\times$ for a set of solutions to Unix puzzles published by Nokia Bell Labs.

Strengths
---------
The main stength of this paper is the ability of PaSh to handle already-existing shell scripts (instead of introducing yet another DSL for parallel programming), and the demonstrated performance improvements on those scripts.

Significant Weaknesses
----------------------
The main contribution of the paper is an application of existing ideas in dataflow/stream processing to the domain of POSIX utilities, and it is not clear that this meets the novelty bar for OSDI. Many of the potentially interesting aspects (such as the annotation language, the code transformation, the parallelization policy) are not discussed in sufficient detail.

Comments for author
-------------------
Thank you for submitting this paper to OSDI. While I found the high-level motivation for the paper compelling, and the evaluation results encouraging, I think that the paper needs substantial revisions and additional work before it would meet the level for OSDI. Please take the following suggestions into account when revising your paper for a future submission.

The most novel contribution in the paper appears to be the annotation record language for describing the parallellizability and argument-dependent semantics for a shell command. Unfortunately, the paper does not provide any details about what it takes to write one of these: the example for `comm` is inscrutable (what do the tuples on the right-hand sides mean?), and the grammar in Appedix A is incomplete. I looked at the additional material for more details and found [this comment](https://github.com/3487fc/pash/blob/f8729eb15fc3b2f0e26cdf8b18d1ecd9fe509b65/compiler/command_categories.py#L5-L11):

```
###
### This file contains the logic that determines the category and the
### input outputs of a command. In the future, this should be replaced
### by a parser of the command classification DSL.
###

## TODO: Make the DSL
```

Since this DSL is one of the main stated contributions of the paper, the apparent fact that it didn't exist at the time of writing is a red flag. It would be genuinely interesting to understand the design of the DSL after it has been fully designed, since the space of shell command inputs is highly unstructured with a long tail of different styles, and designing a language to represent all of these styles is a real challenge. It would also be interesting to see how you define the aggregation functions for "parallelizable pure" commands 

The strategy for extracting parallelism in PaSh is very conservative. For example, I was a little surprised to see `cp` consigned to the "side-effectful" category without any further discussion. Did you encounter scripts where `cp` (or `mv`) was applied to perform a side effect that was "local" or "private" to a parallel activity, e.g. by operating on temporary files or in a temporary directory? Is there some way you could capture this pattern in PaSh, e.g. by modeling the state of the file system in your annotation language? It seems like such a feature would be useful to deal with external programs that do not use stdout (or write to a named argument file) for their output. 

The other aspect I would like to understand is how you translate potentially data-parallel `for` loops in the input script. The example in Figure 1 has a loop over years that seems like a natural candidate for parallelization, but it is not clear how or whether PaSh will attempt to execute the individual iterations in parallel (before combining the output in order). Are there any further restrictions if the loop bounds are given as arguments at runtime (or computed values), instead of in the compilation pass? On a related note, in Table 2, it was not obvious that increasing the degree of parallelism would increase the number of nodes in the graph and have such a dramatic impact on compile time: did you consider using an abstraction in the dataflow graph to represent variable amounts of fan-in and fan-out?

There are a few policy questions that are unevaluated in the paper. For example, using the `relay` node unwisely could exhaust memory, and it wasn't clear how much buffering it performs. (Would a similar outcome be achievable by changing the size of the system pipe buffer? Could you design a better version of `mkfifo` that supports finer-grained configuration here?) Also, did you investigate allowing different static parallelism factors for different parts of the graph, or creating a scheduler that dynamically adjusts the amount of parallelism in different parts of the graph (which might require a different dataflow representation for parallelism)? Finally, how do you propagate size information between nodes so that you can use the optimized implementation of `split` (e.g. is this part of the annotation record DSL)?

Usually I would applaud precision in a systems paper, but description of the dataflow graph model in Section 4 mostly provides details that are mostly irrelevant. I expect that most readers will be familiar with the basic concepts of I/O streams in a modern operating system, and have probably written a few shell scripts in their time, so it's unclear what additional value is obtained from the automata-theory-style definitions for a stream of characters. Similarly, the formal definition of statelessness in terms of a semigroup homomorphism is fine, but not terribly informative. Even with the precision, there seemed to be a disconnect between $D^*$ as a set of *finite* words, and the possibility that an input stream can be unbounded. The description of the aggregation function $agg$ as taking the *concatenation* of two mapped inputs seemed incorrect in the case of `sort`, since the merging version of `sort` definitely does not concatenate its separate inputs. If you developed the formal definitions further - say, into a proof, or a formal semantics - there might be more value here, but it seemed like a poor use of space. 

There are a few more style nits as well:
* Compounding the impression that the paper was rushed, the conclusions end mid-sentence at the end of page 12.
* In the related work section you claim without evident that "most [streaming dataflow] systems perform optimizations that do not preserve semantics". It's not at all clear what you're claiming here.
* Words like "astonishingly", "deceiving" and "unprecedented" are unscientific and rather hyperbolic.
* In a few places you use the term "task-parallel" in a sense that wasn't familiar to me (usually "tasks" are *independent* chunks of work, as in the independent map or reduce tasks in a MapReduce job, whereas the `|` operator induces what is typically called "pipeline-parallel").
* In the "overcoming lazieness" section, the two `grep` commands appear to be missing an input file (or are both mysteriously grepping stdin?).

Novelty
-------
3. New contribution

Experimental methodology
------------------------
4. Good

Writing quality
---------------
2. Needs improvement

Overall merit
-------------
1. Reject (Serious problems, I'll argue against this paper)

Reviewer expertise
------------------
3. Knowledgeable



Review #361F
===========================================================================

Paper summary
-------------
PaSh is a tool that statically analyses and rewrites a Unix shell script to a
parallel version of the original script. To do this, it relies on a database of
annotations/metadata for common POSIX commands, identifying their degree of
parallelisability given different command options. It then translates a script
to a dataflow graph, applying transformations, and rewrites it to a parallel
form by inserting named pipes, background tasks, and splitter/aggregator
processes where required. An evaluation with a range of scripts shows that PaSh
is effective at extracting parallelism from existing pipelines.

Strengths
---------
The topic is perfect for the OSDI audience, the approach (indeed, even the
very problem setting) are novel, and the results obtained are impressive.

Significant Weaknesses
----------------------
Niche applicability: do users actually write shell scripts any more? If they do,
do they consist of complex data-processing pipelines?

Comments for author
-------------------
Thank you for the submission. I enjoyed reading this paper, and suspect the OSDI
audience will too. That said, I suspect the actual audience for the tool will be
rather niche. The premise espoused in the introduction seems to be that there's
a large class of users who write shell scripts but are not experienced
developers. Let's just say I'm skeptical about that claim. Maybe it was true in
the 80s, but today's casual programmer or data analyst is more likely writing
spreadsheet macros than shell scripts. Section 2.1 reads like a requiem for the
glory of the Unix shell; where you (and I!) may see as the "deceiving
succinctness" of a pipeline of tr, cut, grep and sed commands, the average user
is more likely to see an unintuitive, un-debuggable mess.

S3.1: I'm not sure `finger` is a good example these days. How about `ps` or `uptime`?

S3.2: as a purely practical matter, where are the annotations for a command
stored? There must be an extra metadata file, or something like it -- where do
you put it?

S5.1: the definition of parallelisable regions means that many of the throwaway
shell scripts I write would not be parallelisable, because (for ease of
development and debugging) I tend to write them as a series of steps where each
step writes output to a temporary file or environment variable. The example you
give (`cat | grep "foo" > f3 && sort f3`) is a great illustration of this
choice. Ignoring error handling (and, really, who can claim they write shell
scripts with complete and correct error handling?), it's straightforward to
remove the temporary file and move the sort command into the pipeline. Preventing
this is a conservative choice, but might be nevertheless useful.

S6.2: how do you know that the Unix50 solutions were "written by non-experts"?

S6.2: the awk anecdote here is telling -- simpler commands are parallelisable,
but anything with internal state (awk etc.) is not. This is unfortunate, as the
simple commands often have a highly restricted functionality. awk is much nicer
for dealing with input columns than cut or tr, because it does sane field
splitting by default. Maybe this motivates for some knowledge or transformation
of simple awk scripts in the tool?

It's probably worth mentioning PowerShell in the related work. There's no
auto-parallelisation, but the latest version has a parallel foreach that can be
used in a pipeline. Since PowerShell pipelines pass records with named fields
(rather than byte streams) a lot of the fiddly bits of Unix pipelines (e.g.
line-by-line versus all-at-a-time input and record splitting) are avoided. See:
https://devblogs.microsoft.com/powershell/powershell-foreach-object-parallel-feature/
https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_pipelines

The conclusion was cut off mid-sentence. I'm guessing this is why the submission
appears to have been flattened to bitmaps which broke my usual reviewing
workflow, as well as all the internal links. Please avoid that in the future.

Novelty
-------
3. New contribution

Experimental methodology
------------------------
3. Average

Writing quality
---------------
4. Well-written

Overall merit
-------------
4. Accept (Good paper, I will advocate for it)

Reviewer expertise
------------------
3. Knowledgeable
