# PaSh Documentation
Quick Jump: [using pash](#using-pash) | [videos](#videos--video-presentations) | [papers](#academic-papers) 

## Using PaSh

The following resources offer overviews of important PaSh components.

* Short tutorial: [introduction](./tutorial.md#introduction), [installation](./tutorial.md#installation), [execution](./tutorial.md#running-scripts), and [next steps](./tutorial.md#what-next)
* Annotations: [parallelizability](../annotations#main-parallelizability-classes), [study](../annotations#parallelizability-study-of-commands-in-gnu--posix), [example 1](../annotations#a-simple-example-chmod), [example 1](../annotations#another-example-cut), [howto](../annotations#how-to-annotate-a-command)
* Compiler: [intro](../compiler#introduction), [overview](../compiler#compiler-overview), [details](../compiler#zooming-into-fragments), [earlier versions](../compiler#earlier-versions)
* Runtime: [split](../runtime#stream-splitting), [eager](../runtime#eager-stream-polling),  [cleanup](../runtime#cleanup-logic),  [aggregate](../runtime#aggregators)
* Scripts: [one-liners](#common-unix-one-liners), [unix50](#unix-50-from-bell-labs), [weather analysis](#noaa-weather-analysis), [web indexing](#wikipedia-web-indexing)

## Videos & Video Presentations

The following presentations offer short PaSh introductions:

* PaSh: Light-touch Data-Parallel Shell Processing (EuroSys'21: [10-minute](https://www.youtube.com/watch?v=6eSZ1yirnU8) | [20-minute](https://www.youtube.com/watch?v=Je-jt2pNA-s))
* PaSh: A parallelizing shell (Athens PL Seminar 2020: [~40 minutes](https://www.youtube.com/watch?v=UAkfruEvLTk&list=PLdrM8z9GiOahvmZsPn1CXf4EVjy8OA9aq&index=11&t=76s))
* PaSh: A parallelizing shell (POPL 2021 SRC Teaser: [3 minutes](https://www.youtube.com/watch?v=3uqYJo1v1E0))

## Academic Papers & Events

Academic papers, presentations, and other events related to PaSh.

**An Order-aware Dataflow Model for Parallel Unix Pipelines**  
Shivam Handa, Konstantinos Kallas, Nikos Vasilakis, Martin Rinard  
[pdf](https://arxiv.org/pdf/2012.15422.pdf) | acm | [event](https://icfp21.sigplan.org/)
```bibtex
  ..coming soon..
```

**The Future of the Shell: UNIX and beyond**  
Michael Greenberg, Konstantinos Kallas, Nikos Vasilakis  
[pdf](https://fut-shell.github.io/panel-summary.pdf) | acm | [event](https://sigops.org/s/conferences/hotos/2021/#program)
```bibtex
@inproceedings{panel:hotos:21,
  author = {Greenberg, Michael, and Kallas, Konstantinos, and Vasilakis, Nikos},
  title = {The Future of the Shell: UNIX and beyond},
  year = {2021},
  isbn = {9781450384384},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3458336.3465296},
  doi = {10.1145/3458336.3465296},
  booktitle = {Proceedings of the Workshop on Hot Topics in Operating Systems},
  pages = {1–6},
  numpages = {6},
  location = {Ann Arbor, MI, USA},
  series = {HotOS '21}
}
```

**UNIX Shell Programming: The Next 50 Years**  
Michael Greenberg, Konstantinos Kallas, Nikos Vasilakis  
pdf | acm | [event](https://sigops.org/s/conferences/hotos/2021/#program)
```bibtex
@inproceedings{pash:hotos:21,
  author = {Greenberg, Michael, and Kallas, Konstantinos, and Vasilakis, Nikos},
  title = {UNIX Shell Programming: The Next 50 Years},
  year = {2021},
  isbn = {9781450384384},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3458336.3465294},
  doi = {10.1145/3458336.3465294},
  booktitle = {Proceedings of the Workshop on Hot Topics in Operating Systems},
  pages = {1–6},
  numpages = {6},
  location = {Ann Arbor, MI, USA},
  series = {HotOS '21}
}
```

**PaSh: Light-touch Data-Parallel Shell Processing**  
Nikos Vasilakis*, Konstantinos Kallas*, Konstantinos Mamouras, Achilles Benetopoulos, Lazar Cvetković  
[arxiv](https://arxiv.org/pdf/2007.09436.pdf) | acm | video
```bibtex
@inproceedings{pash:eurosys:21,
  author = {Vasilakis, Nikos and Kallas, Konstantinos and Mamouras, Konstantinos and Benetopoulos, Achilles and Cvetkovi\'{c}, Lazar},
  title = {PaSh: Light-Touch Data-Parallel Shell Processing},
  year = {2021},
  isbn = {9781450383349},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3447786.3456228},
  doi = {10.1145/3447786.3456228},
  booktitle = {Proceedings of the Sixteenth European Conference on Computer Systems},
  pages = {49–66},
  numpages = {18},
  keywords = {automatic parallelization, pipelines, Unix, source-to-source compiler, POSIX, shell},
  location = {Online Event, United Kingdom},
  series = {EuroSys '21}
}
```
