# PaSh Setup and Execution Guide

This guide walks you through **setting up and running PaSh** in an isolated environment, including:
- Installing dependencies
- Setting up the virtual environment (`venv`)
- Running `setup-pash.sh` to install required tools
- Exporting `PASH_TOP`
- Running the custom `local-or-pypi.sh` script to configure package installation

---

## **🔹 Step 1: Clone PaSh Repository (If Not Already Cloned)**
If you haven't already cloned PaSh, do so with:

```sh
git clone https://github.com/binpash/pash.git ~/pash
cd ~/pash
```

## **🔹 Step 2: Set Up the Virtual Environment (venv)**
PaSh requires Python dependencies that should be installed inside a virtual environment.

Create and activate a virtual environment:

```sh
python3 -m venv venv
source venv/bin/activate
```


## **🔹 Step 3: Run setup-pash.sh to Install Dependencies and Build Runtime Tools**
PaSh includes a setup script that installs necessary dependencies and compiles missing runtime binaries.

Run:

```sh
./scripts/setup-pash.sh
```



Do not forget to export before using pash: `export PASH_TOP=/home/your-username/pash`
To make this persistent across terminal sessions, add it to your shell configuration file:

```sh
echo 'export PASH_TOP=/home/your-username/pash >> ~/.bashrc
source ~/.bashrc
```

## **🔹 Step 5: Run local-or-pypi.sh to Configure Package Installation**
PaSh provides a script to configure whether dependencies should be installed locally (using a cloned repository) or from PyPI.

To use local dependencies:

```sh
./local-or-pypi.sh local
```

This will:
- Check if `~/annotations` exists.
- If not, clone it from GitHub.
- Update `requirements.txt` to use the local `annotations` repo.

To use PyPI dependencies:

```sh
./local-or-pypi.sh pypi
```


## **🔹 Step 6: Run PaSh**
Once everything is set up, test PaSh with:

```sh
cd ~/pash/evaluation/intro
time $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh
```

### **Expected Output:**

```sh
3712

