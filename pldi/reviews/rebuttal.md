Notes on rebuttal

About networked-distribution experiments: without a location-aware distributed filesystem, the speedup is lowered significantly due to data movement across nodes.
(How to say "there is still speedup, but not as high.")
This is why did not report networked execution; a dish filesystem would be a significant contribution in its own right.
Dish currently comparable to classic distributed operating systems such as Sprite and Amoeba that allow multiple processing servers while storing data on a single networked server; 
  the big difference, of course, is that Dish solves pipeline distribution _without_ the need for a new operating system, simply by (carefully) re-writing pipelines.
