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

## **🔹 Step 2: Run setup-pash.sh and Run setup-annotations.sh**
PaSh includes a setup script that installs necessary dependencies and compiles missing runtime binaries.

Run:

```sh
./pash/scripts/distro-deps.sh
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

Setup the local annotations assuming you have `PASH_TOP` set:
```sh
./setup-annotations.sh
```

You can also provide the path to the PaSh repo directly if `PASH_TOP` is not set. Use -f to force overwrite of the annotations directory if it already exists.
```sh
./setup-annotations.sh ~/pash -f
```

This will:
- Clone the `annotations` repository as a **sibling** to `pash` (i.e., in the same parent directory), or in a specified path if provided.

## **🔹 Step 4: Run PaSh**
Once everything is set up, test PaSh with:

To run it with local annotation:
```sh
cd pash/evaluation/intro
time $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh --local-annotations-dr /path/to/annotations
```

For more details, check the test script here:  
[PaSh Evaluation Tests](https://github.com/binpash/pash/blob/main/evaluation/tests/test_evaluation_scripts.sh)

## **🔹 Step 5: Add your own annotations**
[Link to annotations README](https://github.com/binpash/annotations/blob/main/README.md) on how to add your own custom annotations

