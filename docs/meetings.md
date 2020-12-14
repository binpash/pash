
### Meeting Agenda 12/14

* Started the OSX port; but hit on two types of problems — PaSh's runtime commands are not POSIX. We can fix this, but it will require resources outside . Until then, it might make sense to focus on our measurements on Linux machines. More updates next week.
* Started a simple CI server. The focus is to monitor correctness, not performance: It runs all our benchmarks with small inputs, and the smoosh suite.
* The experience of setting up PaSh on various machines is that it's insanely difficult due to its Python and OCaml dependencies. We need to solve this, because it's really affecting our ability to recruit students to work with the three of us on the project.
* We have two new benchmark sets: (i) modern aliases from GitHub (4,820,352 aliases -> 211,982 with pipes), which involve no expansion and are expected to be interactive (ii) NLP pipelines from the solutions from "unix for poets"



**Issues with dependencies**

Two issues: [73](https://github.com/andromeda/pash/issues/73), [74](https://github.com/andromeda/pash/issues/74)

Problems with Python 3.8
* The fact that we can only run with a single python version is a huge impediment
* https://github.com/pypa/pip/issues/5367
* https://stackoverflow.com/questions/58758447/how-to-fix-module-platform-has-no-attribute-linux-distribution-when-instal

Problems with OCaml
* https://github.com/janestreet/install-ocaml — no repos for debian, requires downloading manually pre-built opam binary
* when trying to install OCaml with the opam binary, you hit bwrap — an external sandboxing dependency not part of OCaml.
* The related opam issue from 2018 (https://github.com/ocaml/opam/issues/3424) even has a couple of posts from Xavier Leory !
