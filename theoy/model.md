01/01
-----------

Dataflow Graph: processes (vertices), channels (edges).
* Processes:
    * Memory:  can have arbitrary memory (up to any _b_) because they can swap; they can have no input (resp. output) channels, called data sources (resp. sinks). 
    * Monotonicity: processes are monotonic (monotonicity means you can't retract output---more input leads to more output.)
  (**continuity?**)
* Channels:
    * maintain FIFO ordering.
    * come into two types: (i) bounded channels that only allow a single element to travel; (ii) unbounded (or bounded by _b_?) that allow arbitrary---up to an arbitrarily large _b_ (similar to _ε_ in privacy?).
    * mixed-synchrony: reading from _any_ channel (both bounded and unbounded) blocks; writing is non-blocking.
* There is a correspondance between processes and channels---a process can serve as a channel and a channel can be a process. Given unbounded-memory processes, one could simulate unbounded channels by adding a "dummy" process in the middle (what's blocking then?); so, just from the channel ID you can't identify whether a process is writing to a blocking or a non-blocking channel---and this is true because `/x/y` might be a socket (finite) or a file (infinite) and you don't know!

Data dependence and Well-behavior:
* Well-behavior is this one-in-one-out property (Arvind) (see `grep -c`; not sure how this works with self-initializing)
* Static dataflow (Edward Lee) says streams bounded in space.
* They are data-independent, but we're data-dependent.

Feedback Loops:
* Feedback loops are possible, if not unavoidable because of aliasing---and many times welcome to express interesting patterns.
* Also, all major dataflow models have feedback (including Kahn Networks and Jack Dennis' dataflow that partly inspired Unix).
* Good news is, we have _β_-bounded channels:-)
* And our feedback is bounded---we have loops but can think of them as epochs.

Input-out Determinism (or Determinacy):
* Both channels and processes are deterministic. Channels are deterministic by definition, as they are FIFO; processes are deterministic because they are sequential programs viewable as functions from sequences to sequences---they always produce the same output for the same input.
* Indeed, behavior is *time-independent*!
* There are three sources of non-determinism which we can exclude for now: (i) a non-deterministic shuffling processes that can be plugged into the graph---i.e., `shuf`, (ii) a few channels such as `/dev/random`, `/dev/urandom` etc. (ii) a `$RANDOM` environment variable
* None of these classes is truly random (they are not used by cryptographers), but it also doesn't seem to affect pipeline scheduling.
* **what about multiple readers on a pipe?**

Labels and tagging: Wrapping channels and processes with labels allows processes to operate on different parts of the stream in parallel, aware of the (size) missing data. Tag-aware merging can be left for the end of the pipeline: _i.e._, you could think of {v_1, v_2, ∅, v_4, ∅, v_6, ...} as {<1, v_1>, <2, v_2>, <3, ∅>, <4, v_4>, <5, ∅>, <6, v_6>, ...}, or equivalently, with holes {<1, v_1>, <2, v_2>, <4, v_4>, <6, v_6>, ...}. However, wrapping primitives with tags a la `--->[untag [grep] tag]--->` (i.e., inserting tags in the data-plane) is difficult, because (i) will mess up the semantics (e.g., `grep -v` can filter more lines than it should), and (ii) performance overhead of wrapping / unwrapping individual elements. To solve this, we do not augment _everything_ with a sequence number; instead, we place metadata only on the _boundaries_ of each "substream", on every split operation, so that we know how to merge.

01/08
-----

* denotational semantics for stream PL?
* can we have a `paste` `join` command?
* We said we will all look into the commands / examples and identify distributed zip.

* Could we have the equivalent of π-calculus or CCS for Unix streams, i.e. ,the Unix-stream calculus?

