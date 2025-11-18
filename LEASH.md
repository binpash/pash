# Leash - PaSh Serverless Execution System

> **📝 IMPORTANT**: This document describes how Leash works. Whenever you learn something new about Leash, you MUST update this file by appending new information or modifying existing sections to keep it accurate and up-to-date.

---

## What is Leash?

**Leash** is PaSh's serverless execution mode that distributes shell script computation across **EC2 (coordinator)** and **AWS Lambda (workers)** for parallel execution. It transforms dataflow graphs into distributed pipelines where:

- **Heavy computation** (sort, grep, etc.) → Offloaded to Lambda workers
- **Coordination** (split, merge, I/O) → Handled by EC2 coordinator
- **Data transfer** → Streaming via TCP sockets (no intermediate files)

---

## Architecture Components

### 1. Compilation (PaSh Compiler)

**Location**: `/home/ubuntu/pash/compiler/serverless/serverless_executor.py`

**Process**:
1. PaSh analyzes shell script and creates dataflow IR (Intermediate Representation)
2. Compiler identifies parallelizable regions
3. Generates distributed execution plan with width parameter (e.g., width=2 means 2-way parallelism)
4. Creates bash scripts for each node in the computation graph
5. Uploads scripts to S3 at `sls-scripts/{folder_id}/{script_id}.sh`

**Key Function**: `prepare_scripts_for_serverless_exec()` in `ir_helper.py`

### 2. Execution Coordinator (EC2)

**EC2 Handler**: `/home/ubuntu/pash/scripts/serverless/ec2-handler.py`

**Port**: 9999 (TCP socket server)

**Responsibilities**:
- Listens for script execution requests
- Downloads scripts from S3
- Executes bash scripts locally
- Handles split/merge operations
- Manages S3 upload/download of data

**Message Format**:
```json
{
  "ids": ["script-uuid-1", "script-uuid-2"],
  "folder_ids": ["timestamp-folder"]
}
```

**Execution Flow**:
1. Receives JSON request via TCP socket
2. For each script ID:
   - Downloads from S3: `s3://{bucket}/sls-scripts/{folder_id}/{script_id}.sh`
   - Saves to `/tmp/script-{folder_id}-{script_id}.sh`
   - Executes: `/bin/bash /tmp/script-{folder_id}-{script_id}.sh {folder_id}`

### 3. Worker Execution (Lambda)

**Lambda Function**: `/home/ubuntu/pash/runtime/serverless/lambda-function.py`

**Function Name**: `lambda` (deployed to AWS Lambda)

**Invocation Type**: `Event` (asynchronous) or `RequestResponse` (synchronous)

**Payload Format**: Same as EC2 handler

**Responsibilities**:
- Receives parallelizable computation tasks
- Downloads script from S3
- Executes computation (sort, grep, etc.)
- Streams results back via network

**Environment**:
- Python 3.9 runtime
- Custom layer with PaSh runtime tools (`dgsh-tee`, etc.)
- Memory: ~1769 MB (configurable)
- Timeout: Depends on workload

### 4. Invocation Interface

**Invoke Script**: `/home/ubuntu/pash/aws/invoke-lambda.py`

**Usage**:
```bash
python3 invoke-lambda.py {script_id} {folder_id} {type}
# type: "lambda" or "ec2"
```

**Lambda Invocation**:
```python
lambda_client.invoke(
    FunctionName="lambda",
    InvocationType="Event",
    Payload=json.dumps({"ids": [id_], "folder_ids": [folder_id]})
)
```

**EC2 Invocation** (to localhost:9999):
```python
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect(("localhost", 9999))
s.sendall(json.dumps({"ids": [id_], "folder_ids": [folder_id]}).encode("utf-8"))
```

---

## Data Flow Mechanisms

### 1. Network Communication (pashlib)

**Binary**: `/opt/pashlib` (Rust implementation)

**Purpose**: Creates point-to-point TCP connections between nodes for streaming data

**Protocol**: UUID-based send/recv pairs

**Format**:
```bash
/opt/pashlib send*{uuid}*{sender_port}*{num_connections}*{fifo_path} \
              recv*{uuid}*{receiver_port}*{num_connections}*{fifo_path}
```

**Example**:
```bash
# Node A sends to Node B
Node A: /opt/pashlib send*49e60090-...*0*1*/tmp/fifo10
Node B: /opt/pashlib recv*49e60090-...*1*0*/tmp/fifo12
```

**How it Works**:
- Matching UUIDs establish TCP connection
- Data flows from sender's FIFO → network → receiver's FIFO
- No intermediate storage, pure streaming

**Port Discovery**:
- Nodes bind to ephemeral ports
- pashlib prints: `Get local_addr: {ip}:{port} external_addr: {public_ip}:{port}`
- Connections established via AWS internal networking or public IPs

### 2. Named Pipes (FIFOs)

**Purpose**: Inter-process communication within a single node

**Pattern**:
```bash
mkfifo /tmp/pash_xxx/hash/#fifo10
command1 > /tmp/pash_xxx/hash/#fifo10 &
command2 < /tmp/pash_xxx/hash/#fifo10 &
```

**Usage**: Connect processes in pipeline (e.g., download → split → send)

### 3. Buffering (dgsh-tee)

**Binary**: `runtime/dgsh-tee`

**Purpose**: Buffer data between pipeline stages to prevent blocking

**Usage**:
```bash
dgsh-tee -i input_fifo -o output_fifo -I -b 5M
```

**Flags**:
- `-i`: Input FIFO
- `-o`: Output FIFO
- `-I`: Enable input mode
- `-b 5M`: 5MB buffer size

**Why Needed**: Prevents deadlocks when processes have mismatched read/write speeds

### 4. S3 Integration

#### Download: `aws/s3-get-object.py`
```python
s3.download_fileobj(BUCKET, object_key, file_handle)
```
**Usage**: `python3.9 aws/s3-get-object.py "inputs/1M.txt" "/tmp/fifo"`
- Streams S3 object directly to FIFO
- No disk storage
- Uses AWS_BUCKET environment variable

#### Upload: `aws/s3-put-object.py`
```python
s3.upload_fileobj(file_handle, BUCKET, object_key)
```
**Usage**: `python3.9 aws/s3-put-object.py "outputs/result.txt" "/tmp/fifo"`
- Reads from FIFO and streams to S3
- No intermediate files

### 5. Data Splitting (r_split)

**Binary**: `runtime/r_split`

**Purpose**: Round-robin split input into N equal chunks

**Usage**:
```bash
r_split -r input_fifo num_lines output_fifo1 output_fifo2 ... output_fifoN
```

**Example**:
```bash
r_split -r /tmp/fifo2 1000000 /tmp/fifo10 /tmp/fifo14
# Splits 1M lines into 2 outputs (~500K each)
```

**Algorithm**: Reads lines from input, distributes round-robin to N outputs

---

## Computation Patterns

### Pattern 1: Split-Process-Merge (Width=N)

**Use Case**: Parallelizable operations (sort, grep, etc.)

**Structure**:
```
EC2: Split input into N chunks
  ↓
Lambda×N: Process chunks in parallel
  ↓
EC2: Merge results
```

**Example (sort, width=2)**:
```
EC2: Download from S3 → r_split into 2 chunks → Send to Lambdas
Lambda 1: Receive → sort → Send back
Lambda 2: Receive → sort → Send back
EC2: Receive both → sort -m (merge) → Upload to S3
```

### Pattern 2: Binary Merge Tree (Width=4+)

**Use Case**: When N > 2, creates hierarchical merges to reduce coordinator load

**Structure (width=4)**:
```
EC2: Split into 4
  ↓
Lambda×4: Process
  ↓
EC2×2: Merge pairs (2 merge nodes)
  ↓
EC2: Final merge
```

**Advantage**: Distributes merge work across multiple EC2 instances

### Pattern 3: Pipeline Stages

**Use Case**: Multi-stage transformations

**Structure**:
```
EC2: Stage 1 → Lambda: Stage 2 → EC2: Stage 3
```

**Example (word frequency)**:
```
EC2: Download → Split
Lambda×N: tr | tr | sort
EC2: Merge → uniq -c → sort -rn → Upload
```

---

## Script Generation

### Script Structure

Every generated script follows this pattern:

```bash
#!/bin/bash
export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1

# PaSh helper functions (pash_communicate_daemon, etc.)
# ...

# Create temp directories
mkdir -p /tmp/pash_xxx/hash/

# FIFO management
rm_pash_fifos() {
  rm -f /tmp/pash_xxx/hash/#fifo1
  rm -f /tmp/pash_xxx/hash/#fifo2
  # ...
}
mkfifo_pash_fifos() {
  mkfifo /tmp/pash_xxx/hash/#fifo1
  mkfifo /tmp/pash_xxx/hash/#fifo2
  # ...
}

rm_pash_fifos
mkfifo_pash_fifos

pids_to_kill=""

# Actual computation pipeline
{ command1 < fifo1 > fifo2 & }
pids_to_kill="${!} ${pids_to_kill}"

{ command2 < fifo2 > fifo3 & }
pids_to_kill="${!} ${pids_to_kill}"

# Wait for all processes
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}

rm_pash_fifos
( exit "${internal_exec_status}" )
```

### Script Components

1. **Environment Setup**: Set PATH to include PaSh runtime tools
2. **Helper Functions**: Daemon communication (usually no-ops in serverless)
3. **FIFO Creation**: Create all named pipes needed for pipeline
4. **Pipeline Execution**: Launch all commands in background with `&`
5. **Process Tracking**: Collect PIDs for cleanup
6. **Synchronization**: Wait for all processes to complete
7. **Cleanup**: Remove FIFOs and exit with status

---

## Execution Flow Example (sort, width=2)

### Step-by-Step

**Initialization**:
1. User runs: `pa.sh --serverless --width 2 sort.sh`
2. PaSh compiles script → generates 4 bash scripts
3. Scripts uploaded to S3: `sls-scripts/{timestamp}/`
4. Orchestration begins

**Execution**:

```
Timeline:

T+0s: EC2 Splitter starts
  - Downloads s3://bucket/inputs/1M.txt to FIFO
  - r_split reads from FIFO, splits into 2 chunks
  - Sends chunk1 via pashlib send*uuid1*
  - Sends chunk2 via pashlib send*uuid2*

T+1s: Lambda 1 & Lambda 2 start (invoked asynchronously)
  - Lambda 1: pashlib recv*uuid1* receives chunk1 → sort → send back via send*uuid3*
  - Lambda 2: pashlib recv*uuid2* receives chunk2 → sort → send back via send*uuid4*

T+2s: Lambdas complete sorting (1.4s execution time typical)

T+3s: EC2 Merger starts
  - pashlib recv*uuid3* receives sorted chunk1
  - pashlib recv*uuid4* receives sorted chunk2
  - sort -m merges both streams
  - Uploads result to s3://bucket/outputs/result.txt

T+5s: Complete
```

### Network Connections

```
EC2 Splitter (port 50000) ←→ Lambda 1 (port 44000)
    UUID: send*49e60090...  ↔  recv*49e60090...

EC2 Splitter (port 50001) ←→ Lambda 2 (port 44001)
    UUID: send*e64098c5...  ↔  recv*e64098c5...

Lambda 1 (port 44100) ←→ EC2 Merger (port 50100)
    UUID: send*f0082947...  ↔  recv*f0082947...

Lambda 2 (port 44101) ←→ EC2 Merger (port 50101)
    UUID: send*61d430ca...  ↔  recv*61d430ca...
```

---

## Configuration & Environment

### Required Environment Variables

**For EC2 & Lambda**:
```bash
AWS_BUCKET="your-s3-bucket"           # S3 bucket for scripts and data
AWS_ACCOUNT_ID="123456789012"         # AWS account ID
AWS_REGION="us-east-1"                # AWS region (default)
```

**For PaSh**:
```bash
PASH_TOP="/home/ubuntu/pash"          # PaSh installation directory
PASH_TMP_PREFIX="/tmp/pash_xxx"       # Temp directory for execution
```

### AWS Resources Required

1. **S3 Bucket**: Store scripts, inputs, outputs
2. **Lambda Function**: Named "lambda", Python 3.9 runtime
   - Must have PaSh runtime layer attached
   - IAM role with S3 read/write permissions
   - Sufficient memory (1769 MB recommended)
3. **EC2 Instance**: Run coordinator and ec2-handler.py
   - IAM role with Lambda invoke + S3 permissions
   - Port 9999 open (for local socket communication)
4. **CloudWatch Logs**: `/aws/lambda/lambda` log group

---

## Debugging & Monitoring

### Log Collection

**Lambda Logs**:
```bash
# Fetch and delete logs
cd /home/ubuntu/pash
python3 scripts/serverless/utils.py output_dir
```
Saves to: `output_dir/logs/log_id_*_scriptid.log`

**Script Collection**:
Scripts automatically saved from S3 to: `output_dir/scripts/{folder_id}/{script_id}.sh`

### Log Analysis

**Parse logs**:
```python
# scripts/serverless/utils.py
analyze_logs("output_dir/logs")
```

**Metrics**:
- Total billed time (milliseconds)
- Number of Lambda invocations
- Cost estimate
- Error detection (Connection reset, Timeout, OOM)

### Common Issues

**"Connection reset by peer"**:
- Lambda cold start timing issues
- Network connectivity problems between Lambda and EC2
- Solution: Add retries, increase timeout

**"Port already in use"**:
- pashlib port conflicts
- Multiple executions interfering
- Solution: Use unique temp directories per run

**"Timeout"**:
- Lambda timeout too short for workload
- Increase Lambda timeout in configuration

**"Out of Memory"**:
- Lambda memory insufficient for data size
- Increase Lambda memory allocation

---

## Performance Characteristics

### Lambda Execution Time

**Typical (1M line sort, width=2)**:
- Initialization: 250-300ms (cold start)
- Execution: 1000-1400ms
- Total billed: 1400-1700ms

**Costs** (us-east-1, 1769MB):
- Per ms: $0.000000028849902
- Per execution: ~$0.000049
- 1000 executions: ~$0.049

### Network Transfer

**S3 Download**: Streams directly to FIFO, bandwidth-limited by S3
**Inter-node**: Uses AWS internal networking when possible
**Throughput**: Depends on data size and network conditions

### Scaling

**Width=2**: 2 Lambda workers, 2-3 EC2 coordinators (split + merge)
**Width=4**: 4 Lambda workers, 4-5 EC2 coordinators (split + 2 L1 merges + final merge)
**Width=N**: N Lambda workers, log₂(N) merge levels

**Bottlenecks**:
- Lambda concurrent execution limit (default: 1000)
- EC2 network bandwidth
- S3 request rate limits

---

## File Locations Reference

### Core Components
```
/home/ubuntu/pash/
├── aws/
│   ├── invoke-lambda.py          # Invoke Lambda/EC2 workers
│   ├── s3-get-object.py          # Download from S3
│   └── s3-put-object.py          # Upload to S3
├── compiler/serverless/
│   ├── serverless_executor.py    # Main orchestration logic
│   └── ir_helper.py              # IR to scripts converter
├── runtime/serverless/
│   ├── lambda-function.py        # Lambda entry point
│   └── aws/invoke-lambda.py      # Runtime version of invoker
├── scripts/serverless/
│   ├── ec2-handler.py            # EC2 socket server (port 9999)
│   └── utils.py                  # Log collection & analysis
└── runtime/
    ├── r_split                   # Round-robin splitter
    ├── dgsh-tee                  # Buffering tool
    └── wait_for_output_and_sigpipe_rest.sh
```

### Generated Artifacts
```
S3: s3://{bucket}/
├── sls-scripts/{folder_id}/
│   ├── {script_id_1}.sh          # Splitter script
│   ├── {script_id_2}.sh          # Lambda worker script
│   └── ...
├── oneliners/inputs/
│   └── 1M.txt                    # Input data
└── oneliners/outputs/
    └── result.txt                # Output data

Local: /tmp/
└── pash_{random}/
    └── {hash}/
        ├── #fifo1                # Named pipes
        └── script-{folder_id}-{script_id}.sh
```

---

## Advanced Topics

### Custom Commands

To add support for new commands:
1. Define aggregator in `compiler/aggregators.py`
2. Add parallelization spec in `compiler/annotations/`
3. Test with: `pa.sh --serverless --width N new_command.sh`

### Hybrid Execution

Mix local and serverless execution:
- Set node execution location in IR
- Some nodes run locally, others on Lambda
- Useful for debugging or partial offloading

### Recovery & Checkpointing

Enable with `--recovery` flag:
- Stores intermediate results in S3
- On failure, restarts from last checkpoint
- Implemented in `serverless_executor.py`

### Batch Invocation

For many small tasks:
- Group multiple script IDs in single invocation
- Reduces Lambda invocation overhead
- Implemented in `invoke_lambda()` batch handling

---

## Future Enhancements

**Potential Improvements**:
1. Dynamic width selection based on input size
2. Lambda-to-Lambda direct communication (bypass EC2)
3. Multi-region execution for geo-distributed data
4. Cost optimization: Lambda vs Fargate vs EC2 Spot
5. Speculative execution for fault tolerance

---

## Version History

- **Last Updated**: 2025-11-07
- **PaSh Version**: Commit eadfb21e (Nov 2024)
- **Leash Status**: Experimental, active development

---

## Related Documentation

- `LEASH_ANALYSIS.md`: Detailed analysis of sort execution with width=2 & 4
- `manual_orchestrate_width2.sh`: Manual orchestration script example
- PaSh main docs: `/home/ubuntu/pash/docs/`

---

**Remember**: Update this document whenever you learn new details about Leash's implementation, architecture, or behavior!
