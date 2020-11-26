
The key insight here is that the key shell abstractions (streams and pipes) are trivially distributed in a scalable and fault-tolerant manner---for example, similar to Spark's RDD's.

# Possible Homework Breakdown

* Create a few categories---trivially distributable (e.g., pure functions), might require some coordination, 

* Break down most primitives 

  * maybe some shells (such as Plan 9's `rc` shell) have primitives that are naturally more amenable to distribution than others.

* What about variable sharing (mostly write once / read many times)

* Other constructs such as `if` and `for`?

* Need to Redirect input / output streams

* We need to look into an extensibility. What if the user, in their environment, have access to primitives not "known" to the shell?


[RaftLib](https://github.com/RaftLib/RaftLib)
