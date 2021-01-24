## Annotations

This is a reference guide for annotations and how to define them.

### Supported Keys

This part of the document describes the supported keys of annotation objects.

#### `command` -- command name

This denotes the command name.

#### `cases` -- sequence of annotation cases

This denotes the sequence of annoation cases, namely a sequence of predicates on command arguments, each of which indicates the inputs and outputs of the command as well as its parallelizability class.

#### `comment` -- comment

Used for text comments in the annotation.

#### `options` -- annotation options

This contains a sequence of options that affect how an annotation is interpreted.
Example values include:
 - "stdin-hyphen", which means that `-` denotes reading from stdin
 - "empty-args-stdin", which means that if the command was not given any input argument files, then it reads its input from stdin

#### `short-long` -- short-long option correspondence

This contains a correspondence between relevant short and long options. It is used to simplify annotation writing (by allowing the user to not duplicate all predicates for both long and short option of the argument).

TODO: Often long options are given their value with `=` and not as a separate argument. Add another option to handle that.

TODO: Extend this to be option-defs, defining options completely.

It should have a list of objects with the following fields:
- short
- long
- value (none, optional, mandatory)
- short value separator ' ' | '=' | ...
- long value separator ' ' | '=' | ...

Defining options properly is necessary if a command reads from non-option arguments. Otherwise it fails (e.g. `cut -d ' ' -f 1` might be considered to read from ' ' and 1).

#### `options-with-arguments` -- options that have a mandatory value

This contains all the options that have a mandatory argument.

TODO: What should we do with optional values?