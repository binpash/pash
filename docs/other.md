## Issues

Commands:

* BSD commands are (often, slightly) different from GNU
* Commands generally have spaces of variable size that affect `agg` and `diff`
* Our runtime, `split` and `eager`, is affected by POSIX vs. optimized Linux choices
* aggregators could be written in (stream-based) Rust
* Commands are context dependent --- for example, `grep` might come from `alias grep = "awk"`. Generally, since PaSh's annotations (contrary to types) have no notion of scope, 
  we have to think about this.

## Broader Impact

Ideas about broader impact can be found in this [document](https://docs.google.com/document/d/1XUAXr-Wt44Z2LLIN4OtK6FAlk-KOCHAs-_tWbKoJQGI/edit?usp=sharing).

