## Notes

If we solve some of these questions below we can include them as
pieces of wisdom in the language description.

### Notes/Questions regarding the language design

- Properly handle `-` as a file argument in commands. - means
  stdin for many commands. Again this can be an option in the
  annotation.
  + I added an option "stdin-hyphen" to indicate this.

- What is the proper and correct way to separate non flag arguments
  from flag arguments. We should make sure that our separation is
  adequate.

- Can multi-input commands (such as general comm, diff) be classified
  as pure?

- How exactly should we handle combined flags (e.g. -13 is the same as
  -1 and -3). Should the annotation include both? Or just separate
  flags and then our annotation component takes care of that? Is this
  convention necessary? If it is a known convention maybe we can
  include it in the annotation of the command (meaning that when a
  command annotation contains an X field, then we can combine flags
  etc. This would make the annotation language simpler.

- How should we handle `grep`'s inputs-outputs? At the moment I am
  doing it with a predicate that checks the args size but this might
  not be the best way.
