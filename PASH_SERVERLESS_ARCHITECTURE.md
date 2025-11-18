# PaSh Serverless (LEASH) Architecture & Integration Guide

**Last Updated:** 2025-11-14
**Purpose:** Complete reference for understanding PaSh serverless compilation, execution, and potential integration points

---

## Table of Contents

1. [PaSh Overview](#pash-overview)
2. [LEASH Serverless Mode](#leash-serverless-mode)
3. [Compilation Flow](#compilation-flow)
4. [Serverless Execution Flow](#serverless-execution-flow)
5. [Holepunch/Pashlib Architecture](#holepunchpashlib-architecture)
6. [Generated Scripts Structure](#generated-scripts-structure)
7. [Integration Points for Reusing Scripts](#integration-points-for-reusing-scripts)
8. [How to Skip Compilation and Reuse Scripts](#how-to-skip-compilation-and-reuse-scripts)

---

## PaSh Overview

**PaSh** = **Pa**rallelizing **Sh**ell - A shell script compiler and runtime that automatically parallelizes Unix shell scripts.

### Core Concept

```
Input: Sequential shell script
       ↓
PaSh Compiler: Analyze dataflow, generate parallel execution plan
       ↓
Output: Distributed execution across workers (EC2 or Lambda)
```

### Two Execution Modes

1. **Local Mode** - Execute on single machine with multiple processes
2. **Serverless Mode (LEASH)** - Execute across EC2 + Lambda workers

---

## LEASH Serverless Mode

**LEASH** = **L**ambda **E**xecution of **A**utomatically-parallelized **Sh**ell scripts

### Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│                      User Input                          │
│                   ./run_leash.sh                         │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│                    pa.sh (entry)                         │
│  Invokes: compiler/pash.py --serverless_exec            │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│              PaSh Compiler (pash.py)                     │
│  - Parse shell script                                    │
│  - Build dataflow graph (IR)                            │
│  - Optimize & parallelize                               │
│  - Generate distributed execution plan                  │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│        ServerlessManager (serverless_executor.py)       │
│  - Generate shell scripts for each node                 │
│  - Upload Lambda scripts to S3                          │
│  - Invoke ec2-handler.py (socket server on EC2)         │
│  - Orchestrate execution with holepunch                 │
└─────────────────────────────────────────────────────────┘
                          ↓
┌──────────────────┐              ┌──────────────────────┐
│  EC2 Scripts     │              │  Lambda Scripts      │
│  (local exec)    │◄────holepunch────►│  (remote exec)   │
│  - Download S3   │              │  - Sort data         │
│  - Split data    │              │  - Process chunks    │
│  - Merge results │              │                      │
└──────────────────┘              └──────────────────────┘
```

### Key Components

1. **`pa.sh`** - Entry point, sets up environment, calls compiler
2. **`compiler/pash.py`** - Main compiler, parses and optimizes scripts
3. **`compiler/serverless/serverless_executor.py`** - ServerlessManager class
4. **`compiler/serverless/ir_helper.py`** - Converts IR to shell scripts
5. **`scripts/serverless/ec2-handler.py`** - Socket server on EC2 for coordination
6. **`aws/invoke-lambda.py`** - Lambda invocation wrapper
7. **`pashlib` (Rust)** - Holepunch NAT traversal library

---

## Compilation Flow

### Step-by-Step Process

#### 1. Entry Point: `pa.sh`

```bash
./pa.sh --serverless_exec \
  --width 2 \
  -c "cat input.txt | sort | wc -l"
```

**What happens:**
- Sets `PASH_TOP=/home/ubuntu/pash`
- Exports environment variables
- Calls: `python3 compiler/pash.py --serverless_exec -c "cat input.txt | sort | wc -l"`

**Key file:** `/home/ubuntu/pash/pa.sh`

---

#### 2. Compiler: `pash.py`

**Location:** `/home/ubuntu/pash/compiler/pash.py`

**Process:**
```python
def main():
    # 1. Parse arguments
    config = parse_args()

    # 2. Parse shell script into AST
    ast = parse_shell(config.command)

    # 3. Build dataflow graph (IR)
    ir_graph = build_ir(ast)

    # 4. Optimize & parallelize
    optimized_ir = optimize(ir_graph, width=config.width)

    # 5. If serverless mode:
    if config.serverless_exec:
        from compiler.serverless.serverless_executor import ServerlessManager
        manager = ServerlessManager(optimized_ir, config)
        manager.execute()
    else:
        # Local execution
        execute_local(optimized_ir)
```

**Output:** Intermediate Representation (IR) - a dataflow graph representing parallel execution plan

---

#### 3. ServerlessManager: Script Generation

**Location:** `/home/ubuntu/pash/compiler/serverless/serverless_executor.py`

**Key method:** `execute()`

```python
class ServerlessManager:
    def __init__(self, ir_graph, config):
        self.ir = ir_graph
        self.width = config.width
        self.bucket = os.getenv('AWS_BUCKET')
        self.session_id = generate_session_id()  # Timestamp

    def execute(self):
        # 1. Generate shell scripts for each IR node
        scripts = self.generate_scripts()

        # 2. Determine which scripts run on EC2 vs Lambda
        ec2_set, lambda_set = self.classify_scripts(scripts)

        # 3. Upload Lambda scripts to S3
        self.upload_lambda_scripts(lambda_set)

        # 4. Start ec2-handler.py (socket server)
        self.start_ec2_handler()

        # 5. Execute main graph script (orchestrator)
        self.run_main_script()
```

**Script classification logic:**

```python
def classify_scripts(self, scripts):
    ec2_set = set()
    lambda_set = set()

    for script_id, script_content in scripts.items():
        # EC2 scripts contain: split, sort -m, s3-put
        if 'split' in script_content or \
           'sort -m' in script_content or \
           's3-put' in script_content:
            ec2_set.add(script_id)
        else:
            lambda_set.add(script_id)

    return ec2_set, lambda_set
```

**Why this classification?**
- **EC2 scripts:** Download input, split data, merge results, upload output
- **Lambda scripts:** Process individual chunks (sort, grep, etc.)

---

#### 4. Script Generation: `ir_helper.py`

**Location:** `/home/ubuntu/pash/compiler/serverless/ir_helper.py`

**Process:**
```python
def ir_to_shell_scripts(ir_graph, session_id):
    """
    Converts IR graph nodes into shell scripts
    """
    scripts = {}

    for node in ir_graph.nodes:
        script_id = node.id  # UUID
        script_path = f"/tmp/scripts/{session_id}/{script_id}.sh"

        # Generate shell commands for this node
        commands = []

        # Example node: Splitter
        if node.type == 'splitter':
            commands.append('#!/bin/bash')
            commands.append(f's3-get s3://{bucket}/{input_key} /tmp/input.txt')
            commands.append(f'split -n l/{width} /tmp/input.txt /tmp/chunk_')
            for i in range(width):
                commands.append(f's3-put /tmp/chunk_{i} s3://{bucket}/chunks/{session_id}/chunk{i}.txt')

        # Example node: Sorter
        elif node.type == 'sort':
            commands.append('#!/bin/bash')
            commands.append(f'CHUNK_ID=$1')
            commands.append(f's3-get s3://{bucket}/chunks/{session_id}/chunk$CHUNK_ID.txt /tmp/input.txt')
            commands.append(f'sort /tmp/input.txt -o /tmp/sorted.txt')
            commands.append(f's3-put /tmp/sorted.txt s3://{bucket}/results/{session_id}/sorted$CHUNK_ID.txt')

        # Example node: Merger
        elif node.type == 'merger':
            commands.append('#!/bin/bash')
            for i in range(width):
                commands.append(f's3-get s3://{bucket}/results/{session_id}/sorted{i}.txt /tmp/sorted{i}.txt')
            commands.append(f'sort -m /tmp/sorted*.txt -o /tmp/final.txt')
            commands.append(f's3-put /tmp/final.txt s3://{bucket}/outputs/{output_key}')

        scripts[script_id] = '\n'.join(commands)

    return scripts
```

**Output:** Dictionary of `{script_id: script_content}`

**Storage location:** Scripts saved to disk and uploaded to S3:
- Local: `/tmp/scripts/{session_id}/{script_id}.sh`
- S3: `s3://{bucket}/sls-scripts/{session_id}/{script_id}.sh` (Lambda scripts only)

---

## Serverless Execution Flow

### Main Orchestration Script

**Generated by:** `ServerlessManager.generate_main_graph_script()`

**Example structure:**
```bash
#!/bin/bash
# Main orchestration script - runs on EC2

PASH_TOP=/home/ubuntu/pash
SESSION_ID=1763144163

# Invoke all nodes in parallel
$PASH_TOP/aws/invoke-lambda.py splitter-id $SESSION_ID ec2 &
$PASH_TOP/aws/invoke-lambda.py sorter1-id $SESSION_ID lambda &
$PASH_TOP/aws/invoke-lambda.py sorter2-id $SESSION_ID lambda &
$PASH_TOP/aws/invoke-lambda.py merger-id $SESSION_ID ec2 &

wait  # Wait for all to complete
```

### Execution Steps

#### 1. Start EC2 Handler

```bash
python3 scripts/serverless/ec2-handler.py &
```

**What it does:**
- Listens on socket (default: `0.0.0.0:9999`)
- Receives script execution requests
- Executes EC2 scripts locally
- Uses **pashlib holepunch** for NAT traversal with Lambda

#### 2. Invoke Lambda Workers

```bash
python3 aws/invoke-lambda.py <script_id> <session_id> lambda
```

**What it does:**
```python
def invoke_lambda(script_id, session_id, exec_type):
    if exec_type == 'lambda':
        # Download script from S3
        script_path = download_from_s3(
            f"sls-scripts/{session_id}/{script_id}.sh"
        )

        # Invoke Lambda function
        lambda_client.invoke(
            FunctionName='pash-lambda-worker',
            InvocationType='Event',  # ASYNC!
            Payload=json.dumps({
                'script_s3_key': f"sls-scripts/{session_id}/{script_id}.sh",
                'session_id': session_id
            })
        )

    elif exec_type == 'ec2':
        # Send to ec2-handler via socket
        send_to_ec2_handler(script_id, session_id)
```

**Key insight:** Lambda invocation is **asynchronous** (`InvocationType='Event'`)

---

## Holepunch/Pashlib Architecture

### The Problem: Lambda NAT Traversal

**Challenge:** Lambda functions are behind NAT. How do EC2 and Lambda communicate directly?

**Solution:** Holepunch NAT traversal using DynamoDB as rendezvous point

### Architecture

```
┌──────────────────────────────────────────────────────┐
│              DynamoDB "rdv" Table                    │
│  (Rendezvous point for peer discovery)              │
│                                                       │
│  Item: {                                             │
│    node_id: "splitter-uuid",                        │
│    ip: "10.0.1.50",                                 │
│    port: 12345,                                     │
│    timestamp: 1763144163                            │
│  }                                                   │
└──────────────────────────────────────────────────────┘
           ▲                              ▲
           │ Register                     │ Register
           │                              │
    ┌──────┴────────┐            ┌───────┴────────┐
    │  EC2 Node     │            │  Lambda Node   │
    │  pashlib      │            │  pashlib       │
    │               │            │                │
    │ 1. Register   │            │ 1. Register    │
    │ 2. Query DDB  │            │ 2. Query DDB   │
    │ 3. Hole punch │◄──TCP────►│ 3. Hole punch  │
    │ 4. Send data  │            │ 4. Recv data   │
    └───────────────┘            └────────────────┘
```

### Pashlib Components

**Location:** `/home/ubuntu/splash-stun-lib/` (separate Rust library)

**Key files:**
- `src/holepunch.rs` - NAT traversal logic
- `src/db_helper.rs` - DynamoDB registration
- `src/network.rs` - TCP connection management

**How it works:**

```rust
// 1. Register self in DynamoDB
fn register_node(node_id: &str, ip: &str, port: u16) {
    dynamodb_client.put_item(
        table_name: "rdv",
        item: {
            "node_id": node_id,
            "ip": ip,
            "port": port,
            "timestamp": now()
        }
    )
}

// 2. Query DynamoDB for peer
fn get_peer_info(peer_id: &str) -> (ip, port) {
    loop {
        item = dynamodb_client.get_item(
            table_name: "rdv",
            key: {"node_id": peer_id}
        )
        if item.exists() {
            return (item.ip, item.port)
        }
        sleep(500ms)
    }
}

// 3. Attempt NAT traversal (simultaneous TCP open)
fn holepunch_connect(peer_ip: &str, peer_port: u16) -> TcpStream {
    // Both sides attempt to connect to each other
    // simultaneously, creating NAT hole
    let stream = TcpStream::connect((peer_ip, peer_port))
        .with_timeout(Duration::from_secs(15))
        .expect("Holepunch failed");

    stream
}
```

### Timing Requirements

**Critical timing window:** ±15 seconds

1. **EC2 script starts** → Registers in DynamoDB → Waits for Lambda
2. **Lambda invoked (async)** → Cold start (200ms-10s) → Registers in DynamoDB
3. **Both nodes query DynamoDB** → Find each other
4. **Simultaneous TCP connect** → Create NAT hole
5. **Data transfer** over direct TCP connection

**Problem:** Lambda cold starts are unpredictable (200ms-10s), often exceeding NAT traversal timeout window.

### Why Holepunch Failed in Our Tests

| Issue | Impact |
|-------|--------|
| **Async Lambda invocation** | Lambda starts seconds after EC2 |
| **Cold start variability** | 200ms-10s unpredictable delay |
| **Timing window** | NAT traversal requires ±15s alignment |
| **No retries** | Single attempt, no exponential backoff |
| **Opaque errors** | Just "Connection timed out" - hard to debug |

**Success rate:** 0% (all 10+ attempts failed)

---

## Generated Scripts Structure

### Typical Script Layout

**Example:** Sort with width=2

```
/tmp/scripts/1762478811/
├── 559c679d-49ed-49c3-9a1f-08c3ce9e5f02.sh  # Splitter (EC2)
├── 7422cb0d-628b-49ae-b29d-387640f949f6.sh  # Sorter 1 (Lambda)
├── 6b9eae70-3757-42a4-a3fb-eee1dc6a73ca.sh  # Sorter 2 (Lambda)
└── 839ba2b4-ed80-4864-a800-04694e745f76.sh  # Merger (EC2)
```

### Splitter Script (EC2)

```bash
#!/bin/bash
# 559c679d-49ed-49c3-9a1f-08c3ce9e5f02.sh

PASH_TOP=/home/ubuntu/pash
FOLDER_ID=$1  # Passed as argument

# Download input from S3
aws s3 cp s3://inout741448956691/oneliners/inputs/1M.txt /tmp/input.txt

# Split into 2 chunks (round-robin lines)
# PaSh uses custom split logic for even distribution
python3 $PASH_TOP/scripts/split_round_robin.py \
  /tmp/input.txt \
  /tmp/chunk_ \
  2

# Upload chunks to S3
aws s3 cp /tmp/chunk_0 s3://inout741448956691/chunks/$FOLDER_ID/chunk0.txt
aws s3 cp /tmp/chunk_1 s3://inout741448956691/chunks/$FOLDER_ID/chunk1.txt

# Signal completion via holepunch to downstream nodes
pashlib_signal_done "splitter-done"
```

### Sorter Script (Lambda)

```bash
#!/bin/bash
# 7422cb0d-628b-49ae-b29d-387640f949f6.sh

PASH_TOP=/home/ubuntu/pash
FOLDER_ID=$1
WORKER_ID=0

# Wait for splitter to signal ready via holepunch
pashlib_wait_for_signal "splitter-done"

# Download chunk from S3
aws s3 cp s3://inout741448956691/chunks/$FOLDER_ID/chunk$WORKER_ID.txt /tmp/input.txt

# Sort
sort /tmp/input.txt -o /tmp/sorted.txt

# Upload result
aws s3 cp /tmp/sorted.txt s3://inout741448956691/results/$FOLDER_ID/sorted$WORKER_ID.txt

# Signal completion
pashlib_signal_done "sorter-$WORKER_ID-done"
```

### Merger Script (EC2)

```bash
#!/bin/bash
# 839ba2b4-ed80-4864-a800-04694e745f76.sh

PASH_TOP=/home/ubuntu/pash
FOLDER_ID=$1

# Wait for both sorters via holepunch
pashlib_wait_for_signal "sorter-0-done"
pashlib_wait_for_signal "sorter-1-done"

# Download sorted chunks
aws s3 cp s3://inout741448956691/results/$FOLDER_ID/sorted0.txt /tmp/sorted0.txt
aws s3 cp s3://inout741448956691/results/$FOLDER_ID/sorted1.txt /tmp/sorted1.txt

# Merge
sort -m /tmp/sorted0.txt /tmp/sorted1.txt -o /tmp/final.txt

# Upload final result
aws s3 cp /tmp/final.txt s3://inout741448956691/oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt
```

---

## Integration Points for Reusing Scripts

### Where to Hook In

```
┌─────────────────────────────────────────────────────────┐
│              Normal PaSh Flow                            │
├─────────────────────────────────────────────────────────┤
│  pa.sh → pash.py → ServerlessManager → execute()       │
│                                                          │
│  ┌──────────────────────────────────────┐              │
│  │ 1. Parse script (parse_shell)        │              │
│  └──────────────────────────────────────┘              │
│           ↓                                             │
│  ┌──────────────────────────────────────┐              │
│  │ 2. Build IR (build_ir)               │              │
│  └──────────────────────────────────────┘              │
│           ↓                                             │
│  ┌──────────────────────────────────────┐              │
│  │ 3. Optimize (optimize)               │              │
│  └──────────────────────────────────────┘              │
│           ↓                                             │
│  ┌──────────────────────────────────────┐              │
│  │ 4. Generate scripts (ir_to_scripts) │ ◄─── CACHE   │
│  └──────────────────────────────────────┘              │
│           ↓                                             │
│  ┌──────────────────────────────────────┐              │
│  │ 5. Execute (ServerlessManager)       │ ◄─── HOOK    │
│  └──────────────────────────────────────┘              │
└─────────────────────────────────────────────────────────┘

         Integration Points:
         ══════════════════════

1. CACHE: Save generated scripts to reuse later
         Location: /tmp/scripts/{session_id}/

2. HOOK:  Skip steps 1-4, directly execute with:
         - Load cached scripts
         - Use S3-based orchestration
         - Skip holepunch
```

### Option A: Modify ServerlessManager

**Goal:** Add `--reuse-scripts` flag to skip compilation

```python
# compiler/serverless/serverless_executor.py

class ServerlessManager:
    def __init__(self, ir_graph, config):
        self.config = config
        self.scripts_dir = config.get('scripts_dir')  # Pre-generated scripts

    def execute(self):
        if self.config.get('reuse_scripts'):
            # Load pre-generated scripts
            scripts = self.load_scripts_from_dir(self.scripts_dir)
        else:
            # Normal flow: generate from IR
            scripts = self.generate_scripts_from_ir()

        # Continue with execution
        self.run_s3_orchestration(scripts)  # NEW: Use S3 instead of holepunch
```

**Usage:**
```bash
pa.sh --serverless_exec \
  --reuse-scripts /tmp/scripts/1762478811 \
  --width 2 \
  --orchestration s3  # NEW flag
```

### Option B: Standalone Orchestrator (What We Built)

**Goal:** Completely bypass PaSh, directly orchestrate pre-generated scripts

**Current implementation:**
- `manual_s3_orchestrator.py` - Full version
- `manual_s3_orchestrator_no_split.py` - No-split version
- `manual_s3_orchestrator_byte_ranges.py` - Byte-range version

**Future enhancement:** Parse PaSh-generated scripts and execute them

```python
# pash_script_executor.py

def execute_pash_scripts(scripts_dir, session_id, width):
    """
    Execute pre-generated PaSh scripts using S3 orchestration
    """
    # 1. Load all scripts
    scripts = load_scripts(scripts_dir)

    # 2. Classify scripts (EC2 vs Lambda)
    ec2_scripts = [s for s in scripts if is_ec2_script(s)]
    lambda_scripts = [s for s in scripts if not is_ec2_script(s)]

    # 3. Convert scripts to our S3-based format
    converted_scripts = convert_scripts_to_s3(scripts)

    # 4. Execute using our orchestrator
    orchestrator = S3Orchestrator(
        scripts=converted_scripts,
        width=width,
        bucket=os.getenv('AWS_BUCKET')
    )
    orchestrator.execute()

def is_ec2_script(script_content):
    """Determine if script should run on EC2"""
    keywords = ['split', 'sort -m', 's3-put', 'merge']
    return any(kw in script_content for kw in keywords)
```

---

## How to Skip Compilation and Reuse Scripts

### Method 1: Cached Scripts Directory

**Scenario:** You've already run PaSh once, scripts are in `/tmp/scripts/{session_id}/`

```bash
# 1. Run PaSh once to generate scripts (save session_id)
./run_leash.sh sort.sh 1M.txt 2

# Scripts saved to: /tmp/scripts/1762478811/

# 2. Copy scripts to persistent location
mkdir -p /home/ubuntu/pash/cached_scripts/sort_width2
cp -r /tmp/scripts/1762478811/* /home/ubuntu/pash/cached_scripts/sort_width2/

# 3. Reuse scripts with our orchestrator
python3 pash_script_executor.py \
  --scripts-dir /home/ubuntu/pash/cached_scripts/sort_width2 \
  --session-id $(date +%s) \
  --width 2 \
  --orchestration s3
```

### Method 2: Modify pa.sh to Skip Compilation

**Add flag to pa.sh:**

```bash
# pa.sh

REUSE_SCRIPTS=""
SCRIPTS_DIR=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --reuse-scripts)
      REUSE_SCRIPTS="--reuse-scripts"
      SCRIPTS_DIR="$2"
      shift 2
      ;;
    # ... other flags
  esac
done

if [ -n "$REUSE_SCRIPTS" ]; then
  # Skip to execution, don't recompile
  python3 compiler/serverless/execute_cached.py \
    --scripts-dir "$SCRIPTS_DIR" \
    --width "$WIDTH" \
    --orchestration s3
else
  # Normal compilation flow
  python3 compiler/pash.py "$@"
fi
```

**Usage:**
```bash
./pa.sh --serverless_exec \
  --reuse-scripts /home/ubuntu/pash/cached_scripts/sort_width2 \
  --width 2 \
  --orchestration s3
```

### Method 3: Direct Script Conversion

**Convert PaSh scripts to our S3 orchestrator format:**

```python
# convert_pash_to_s3.py

def convert_pash_script_to_s3(pash_script, worker_id):
    """
    Convert PaSh-generated script to S3-based Lambda worker
    """
    # Parse PaSh script
    lines = pash_script.split('\n')

    # Extract operation (sort, grep, etc.)
    operation = detect_operation(lines)

    # Generate S3-based Lambda script
    lambda_script = f"""
import boto3
import subprocess
import os

def lambda_handler(event, context):
    s3 = boto3.client('s3')

    # Download chunk
    bucket = event['bucket']
    chunk_key = event['chunk_key']
    s3.download_file(bucket, chunk_key, '/tmp/input.txt')

    # Execute operation: {operation}
    subprocess.run(['{operation}', '/tmp/input.txt', '-o', '/tmp/output.txt'])

    # Upload result
    result_key = event['result_key']
    s3.upload_file('/tmp/output.txt', bucket, result_key)

    return {{'statusCode': 200}}
"""

    return lambda_script

def convert_all_scripts(pash_scripts_dir):
    """Convert entire PaSh script directory"""
    scripts = load_all_scripts(pash_scripts_dir)

    for script_id, script_content in scripts.items():
        if is_lambda_script(script_content):
            converted = convert_pash_script_to_s3(script_content, script_id)
            save_lambda_script(script_id, converted)
```

**Usage:**
```bash
python3 convert_pash_to_s3.py \
  --pash-scripts /tmp/scripts/1762478811 \
  --output /home/ubuntu/pash/s3_lambdas/

# Deploy converted Lambdas
./deploy_converted_lambdas.sh

# Run orchestrator
python3 manual_s3_orchestrator.py \
  --bucket inout741448956691 \
  --input oneliners/inputs/1M.txt \
  --output oneliners/outputs/result.txt \
  --workers 2
```

---

## Implementation Roadmap

### Phase 1: Script Caching (Low Effort)

**Goal:** Cache PaSh-generated scripts for reuse

```bash
# 1. Run PaSh with caching enabled
PASH_CACHE_DIR=/home/ubuntu/pash_cache ./run_leash.sh sort.sh 1M.txt 2

# 2. Reuse cached scripts
./run_leash.sh --reuse-cache /home/ubuntu/pash_cache/sort_1M_width2
```

**Files to modify:**
- `pa.sh` - Add `--cache-dir` flag
- `compiler/pash.py` - Add script caching logic
- `compiler/serverless/serverless_executor.py` - Load from cache

### Phase 2: S3-Based Execution (Medium Effort)

**Goal:** Replace holepunch with S3 orchestration in PaSh

```bash
./run_leash.sh --orchestration s3 sort.sh 1M.txt 2
```

**Files to modify:**
- `compiler/serverless/serverless_executor.py` - Add `S3Orchestrator` class
- `compiler/serverless/ir_helper.py` - Generate S3-compatible scripts
- `scripts/serverless/ec2-handler.py` - Remove holepunch dependency

### Phase 3: Hybrid Approach (High Effort)

**Goal:** Support both holepunch and S3 orchestration

```bash
# Holepunch mode (original)
./run_leash.sh --orchestration holepunch sort.sh 1M.txt 2

# S3 mode (new)
./run_leash.sh --orchestration s3 sort.sh 1M.txt 2

# Auto-select based on reliability
./run_leash.sh --orchestration auto sort.sh 1M.txt 2
```

**Benefits:**
- Fallback to S3 if holepunch fails
- Use holepunch for low-latency workloads
- Use S3 for reliability

---

## Key Files Reference

### PaSh Core
- `/home/ubuntu/pash/pa.sh` - Entry point
- `/home/ubuntu/pash/compiler/pash.py` - Main compiler
- `/home/ubuntu/pash/compiler/parse_to_ast.py` - Shell parser
- `/home/ubuntu/pash/compiler/ir.py` - IR definitions

### Serverless Components
- `/home/ubuntu/pash/compiler/serverless/serverless_executor.py` - ServerlessManager
- `/home/ubuntu/pash/compiler/serverless/ir_helper.py` - IR to shell scripts
- `/home/ubuntu/pash/scripts/serverless/ec2-handler.py` - EC2 socket server
- `/home/ubuntu/pash/aws/invoke-lambda.py` - Lambda invocation wrapper

### Holepunch Library
- `/home/ubuntu/splash-stun-lib/src/holepunch.rs` - NAT traversal
- `/home/ubuntu/splash-stun-lib/src/db_helper.rs` - DynamoDB coordination
- `/home/ubuntu/splash-stun-lib/src/network.rs` - TCP networking

### Configuration
- `/home/ubuntu/pash/runtime/serverless/serverless.yml` - Lambda deployment config
- `/home/ubuntu/pash/runtime/config.yaml` - PaSh runtime config
- `/home/ubuntu/pash/.pash_init` - Environment setup

---

## Environment Setup

```bash
# Required environment variables
export PASH_TOP=/home/ubuntu/pash
export AWS_BUCKET=inout741448956691
export AWS_ACCOUNT_ID=741448956691
export AWS_REGION=us-east-1

# DynamoDB table for holepunch
export DYNAMODB_TABLE=rdv

# Lambda function name
export LAMBDA_FUNCTION=pash-lambda-worker

# S3 paths
export S3_SCRIPTS_PREFIX=sls-scripts
export S3_CHUNKS_PREFIX=chunks
export S3_RESULTS_PREFIX=results
```

---

## Next Steps

1. **Extract and analyze cached scripts** from existing PaSh runs
2. **Implement script parser** to understand PaSh script structure
3. **Create conversion layer** from PaSh scripts to S3 orchestrator
4. **Add caching to PaSh** to avoid recompilation
5. **Contribute S3 orchestration back to PaSh** as alternative execution mode

---

**For new sessions:** Read this file along with `S3_ORCHESTRATOR_README.md` for complete context on both PaSh internals and our S3-based orchestration solution.
