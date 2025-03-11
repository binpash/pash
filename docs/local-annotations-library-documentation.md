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

Setup the local annotations:
```sh
./setup-annotations.sh
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

# How to Add a New Command to PaSh's Annotation System

This guide explains how to manually create InputOutputInfo and ParallelizabilityInfo generators for a new command in PaSh.


### **Step 1: Add the Command to the Name Mapping**

PaSh uses a dictionary to map shell command names to Python class names. To register a new command:

1. Open **`AnnotationGeneration.py`** (or the relevant file where `DICT_CMD_NAME_TO_REPRESENTATION_IN_MODULE_NAMES` is defined).
2. Add an entry for the new command:

   ```python
   DICT_CMD_NAME_TO_REPRESENTATION_IN_MODULE_NAMES = {
       ...
       "<command-name>": "<ClassRepresentation>",  # Add your new command here
       ...
   }
### **Step 2: Create the InputOutputInfo and ParallelizabilityInfo Generator files**

Each command requires an InputOutputInfo and ParallelizabilityInfo generator, which determines how it handles input and output files.

Navigate to:

```sh
pash_annotations/annotation_generation/annotation_generators/
```

Create a two files named:

 ```sh
InputOutputInfoGenerator<ClassRepresentation>.py
ParallelizabilityInfo<ClassRepresentation>.py
```

If your command is "cat-wrapper", the file should be:

```sh
InputOutputInfoGeneratorCatWrapper.py
ParallelizabilityInfoCatWrapper.py
```


### **Step 3: Implement the InputOutputInfo and ParallelizabilityInfo Generators**

Inside the newly created files, define a class that **inherits from the appropriate interface**:

- **For Input/Output Behavior**: Inherit from `InputOutputInfoGeneratorInterface`
- **For Parallelization Behavior**: Inherit from `ParallelizabilityInfoGeneratorInterface`

### **InputOutputInfo Generator**
In the **InputOutputInfo generator**, specify how your command processes input and produces output. This includes:
- How the command **reads input**.
- Whether it **writes to stdout** or modifies files in place.
- How **each flag affects input and output behavior**.

For example:
- Commands like `cat` read from **stdin** or files and write to **stdout**.
- Commands like `mv` modify files **in place** without stdout output.
- Commands like `grep` take both **input files and options** that affect behavior.

---

### **ParallelizabilityInfo Generator**
In the **ParallelizabilityInfo generator**, define what parallelization strategies can be applied while maintaining **correct execution**. Consider:
- Whether the command can process input **in independent chunks** (e.g., `sort` can, but `grep` with `-A` or `-B` cannot).
- Whether it can be **executed in parallel** on separate input files.
- Whether it requires **ordering constraints** to maintain correctness.

For example:
- `sort` can process chunks **independently**, then merge results.
- `wc` can process chunks **independently**, and would then sum up the results.
- `cat` with no flags is **stateless**, so the default options work.

By implementing these details, you ensure **efficient parallel execution** while preserving the **functional correctness** of your command.



