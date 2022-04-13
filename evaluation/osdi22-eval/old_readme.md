##  "Artifact Reproduced" badge: 

To reproduce the results presented in the paper:
1. Follow the [installation instructions](https://github.com/binpash/pash/blob/main/docs/tutorial/tutorial.md#installation) to set up PaSh on your infrastructure
2. Make sure you have the `PASH_TOP` environment variable set up pointing to the top-level directory in the repo: `echo $PASH_TOP`.
3. Some of the benchmarks have additional library/package dependencies that may require `sudo` installation. You may install those dependencies manually prior to the script execution (so they don't interrupt normal execution):
`
$PASH_TOP/evaluation/benchmarks/web-index/input/install-deps.sh
$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/install-deps.sh
`
4. For the plot generation, additional packages are required (optional):
```sh
# install the R environment and libraries
# follow the guide here https://cran.r-project.org/bin/linux/ubuntu/

# install ggplot2 for R
# needed for generating the pdf
sudo R
> install.packages('ggplot2', dep = TRUE)
> q()
> n
```
4. Run the following commands:

```sh
cd $PASH_TOP/evaluation/eval_script
# run all the benchmarks 3 times
bash run_all.sh
```

This script will run all the benchmarks included in the evaluation and gather all the execution results for:
```
Classics
Unix50
Covid-mts
NLP
AvgTemp
WebIndex
Dependency Untrangling
```

```sh
## OPTIONAL
# it will format the generated results for the plot generation
# the final source data will be called final.csv
bash gen_data.sh
# it will generate the figures for all the results
./generate_charts.R final.csv
```
The results include the sequential baseline (running `bash`), `PaSh JIT` and full `PaSh AOT` performing all transformations and flag configurations. The scripts are configured to apply with a parallelism `--width` of `16`.

### Notes 

- All the input data has been scaled down with respect to the ones shown in the paper so that they execute in a reasonable amount of time.

- Note on hardware and software requirements: Running these scripts will require significant disk space (>100GB), it will take considerable time (several hours), and will need many CPUs (ideally more than 32). 


### Expected Results
- All the execution logs may be found at `$PASH_TOP/evaluation/eval_results/run{1,2,3}`

- A collection of all the experiments with all the configurations are available for future reference  [here](https://docs.google.com/spreadsheets/d/1flDa2H7FplJBq7JiKAN_7aGfEV3SaPVzbTLrGLjNIpA/edit#gid=0)
