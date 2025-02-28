# PaSh Setup and Execution Guide

This guide walks you through **setting up and running PaSh**, including:
- Installing dependencies
- Running `setup-pash.sh` to install required tools
- Exporting `PASH_TOP`
- Running the `setup-annotations.sh` script to configure package installation with an optional `--local-annotations` flag

---

## **🔹 Step 1: Clone PaSh Repository (If Not Already Cloned)**
If you haven't already cloned PaSh, do so with:

```sh
git clone https://github.com/binpash/pash.git
cd pash
```

## **🔹 Step 2: Run setup-pash.sh to Install Dependencies and Build Runtime Tools**
PaSh includes a setup script that installs necessary dependencies and compiles missing runtime binaries.

Run:

```sh
./scripts/setup-pash.sh
```

Before using PaSh, set the `PASH_TOP` environment variable to the directory where PaSh is stored:

```sh
export PASH_TOP=/path/to/pash
```

To make this persistent across terminal sessions, add it to your shell configuration file:

```sh
echo 'export PASH_TOP=/path/to/pash' >> ~/.bashrc
source ~/.bashrc
```

## **🔹 Step 3: Run setup-annotations.sh to Configure Package Installation**
PaSh provides a script to configure dependencies, including an option to install local annotations.

To install dependencies normally:

To install local annotations:
```sh
./setup-local.sh
```

This will:
- Clone the `annotations` repository as a **sibling** to `pash` (i.e., in the same parent directory), or in a specified path if provided.
- Update `requirements.txt` to add the local `annotations` repo.

## **🔹 Step 4: Run PaSh**
Once everything is set up, test PaSh with:

```sh
cd pash/evaluation/intro
time $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh
```

To run it with local annotation:
```sh
cd pash/evaluation/intro
time $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh --local-annotations
```

You can also run the full PaSh test suite:
```sh
$PASH_TOP/evaluation/tests/test_evaluation_scripts.sh
```
For more details, check the test script here:  
[PaSh Evaluation Tests](https://github.com/binpash/pash/blob/main/evaluation/tests/test_evaluation_scripts.sh)


