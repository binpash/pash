# PaSh Holepunch Integration Guide

**Date:** 2025-11-18
**Purpose:** Complete guide to understanding and running PaSh serverless execution with holepunch NAT traversal

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Holepunch Library (pashlib)](#holepunch-library-pashlib)
3. [Script Execution Flow](#script-execution-flow)
4. [How to Run the 4 Scripts](#how-to-run-the-4-scripts)
5. [Step-by-Step Execution Guide](#step-by-step-execution-guide)
6. [Debugging and Troubleshooting](#debugging-and-troubleshooting)

---

## Architecture Overview

### Components

```
┌─────────────────────────────────────────────────────────┐
│              DynamoDB "rdv" Table                        │
│  Rendezvous point for NAT traversal coordination       │
│                                                          │
│  Entry structure:                                       │
│  {                                                      │
│    "key": "49e60090-...",  // rdv_key (UUID)          │
│    "N_0": "{...EC2 addr...}",  // Node 0 address      │
│    "N_1": "{...Lambda addr...}" // Node 1 address     │
│  }                                                      │
└─────────────────────────────────────────────────────────┘
           ▲                                  ▲
           │ Register                         │ Register
           │ IP:Port                          │ IP:Port
           │                                  │
    ┌──────┴──────────┐              ┌───────┴────────┐
    │  EC2 Instance   │              │ Lambda Worker  │
    │  (this machine) │              │                │
    │                 │              │                │
    │  1. Get STUN    │              │  1. Get STUN   │
    │     address     │              │     address    │
    │  2. Register in │              │  2. Register   │
    │     DynamoDB    │              │     in DDB     │
    │  3. Query peer  │              │  3. Query peer │
    │  4. TCP connect │◄──holepunch──│  4. TCP conn   │
    │  5. Send data   │     via      │  5. Recv data  │
    │     through     │   TCP with   │     through    │
    │     fifos       │  simultaneous│     fifos      │
    │                 │     open     │                │
    └─────────────────┘              └────────────────┘
```

### Key Infrastructure

- **EC2 Instance:** `3.133.71.249` (public), `172.31.12.152` (private)
- **DynamoDB Table:** `rdv` in `us-east-1`
- **Lambda Function:** `lambda` (python3.9)
- **Lambda Layer:** `pashlib:1` (contains `/opt/pashlib` binary)
- **S3 Bucket:** `$AWS_BUCKET` (for script storage)

---

## Holepunch Library (pashlib)

### Location and Source

- **Binary:** `/opt/pashlib` (on both EC2 and Lambda)
- **Source:** `/home/ubuntu/splash-stun-lib/bin/pashlib-oneproc/main.rs`
- **Language:** Rust with Tokio async runtime
- **AWS SDK:** Uses DynamoDB client for coordination

### Command Format

The pashlib binary accepts multiple send/recv commands in a single invocation:

```bash
/opt/pashlib \
  send*<rdv_key>*<me>*<peer>*<fifo_path> \
  recv*<rdv_key>*<me>*<peer>*<fifo_path> \
  ...
```

**Parameters:**
- `rdv_key`: UUID identifying this connection (e.g., `49e60090-bea3-43fe-b3a2-039491f8f37c`)
- `me`: My node ID (0 for EC2, 1 for Lambda)
- `peer`: Peer node ID (1 for Lambda from EC2's perspective, 0 from Lambda's)
- `fifo_path`: Path to named pipe for data I/O

**Example from splitter script:**
```bash
/opt/pashlib \
  send*49e60090-bea3-43fe-b3a2-039491f8f37c*0*1*/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo10 \
  send*e64098c5-470e-4972-9803-5b7b13c16ffd*0*1*/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo14
```

This means:
- **Connection 1:** Send data from fifo10 to peer (node 1) via rdv_key `49e60090...`
- **Connection 2:** Send data from fifo14 to peer (node 1) via rdv_key `e64098c5...`

### How Holepunch Works

#### Step 1: STUN Discovery
```rust
// Get local and external addresses via STUN
let (local_addr, external_addr) = stun_helper::get_addr().await;
// Example output:
// local_addr: 172.31.12.152:54744
// external_addr: 3.133.71.249:54744
```

#### Step 2: DynamoDB Registration
```rust
// Register my address in DynamoDB
let ctx = PashCtx::new(me, rdv_key, local_addr, external_addr).await;

// This creates/updates DynamoDB item:
// {
//   "key": "<rdv_key>",
//   "N_0": "{\"local_addr\":\"172.31.12.152:54744\",\"external_addr\":\"3.133.71.249:54744\"}"
// }
```

DynamoDB attribute naming: `N_<node_id>` (e.g., `N_0`, `N_1`)

#### Step 3: Peer Discovery
```rust
// Wait for peer to register (polls every 1 second)
loop {
    let peer_info = get_attr(&client, &rdv_key, &format!("N_{}", peer)).await;
    if let Some(info) = peer_info {
        let peer_ctx: PashCtx = serde_json::from_str(&info).unwrap();
        break peer_ctx.external_addr;
    }
    tokio::time::sleep(Duration::from_secs(1)).await;
}
```

#### Step 4: Simultaneous TCP Open (Holepunch)
```rust
// Both sides attempt to connect to each other simultaneously
let sock = TcpSocket::new_v4().unwrap();
sock.set_reuseport(true).unwrap();  // Critical for holepunch!
sock.bind(local_addr.parse().unwrap()).unwrap();
let stream = sock.connect(peer_external_addr.parse().unwrap()).await.unwrap();
stream.set_nodelay(true).unwrap();
```

**Key insight:** Both sides bind to their STUN-discovered port and connect to each other simultaneously. This creates a "hole" in the NAT that allows bidirectional communication.

#### Step 5: Data Transfer
```rust
// Send mode: read from fifo, write to TCP stream
let (mut rd, mut wr) = stream.into_split();
let mut buf = [0u8; 128];
let mut fifo = File::open(fifo_path).await.unwrap();
rd.read_exact(&mut buf).await.unwrap();  // Handshake
io::copy(&mut fifo, &mut wr).await.unwrap();

// Recv mode: read from TCP stream, write to fifo
let (mut rd, mut wr) = stream.into_split();
let buf = [0u8; 128];
let mut fifo = File::create(fifo_path).await.unwrap();
wr.write_all(&buf).await.unwrap();  // Handshake
io::copy(&mut rd, &mut fifo).await.unwrap();
```

---

## Script Execution Flow

### The 4 Scripts in Your Directory

Location: `/home/ubuntu/pash/evaluation/benchmarks/oneliners/logs/sort.sh:1M.txt:2/scripts/1762478811/`

1. **559c679d-49ed-49c3-9a1f-08c3ce9e5f02.sh** - Splitter (EC2)
2. **7422cb0d-628b-49ae-b29d-387640f949f6.sh** - Lambda Worker 1
3. **6b9eae70-3757-42a4-a3fb-eee1dc6a73ca.sh** - Lambda Worker 2
4. **839ba2b4-ed80-4864-a800-04694e745f76.sh** - Merger (EC2)

### Dataflow

```
[S3 Input]
    ↓
[EC2: Splitter Script]
    ↓ (via pashlib holepunch)
    ├─→ [Lambda Worker 1: sort] ─┐
    │                              │
    └─→ [Lambda Worker 2: sort] ─┤
                                  ↓ (via pashlib holepunch)
                          [EC2: Merger Script]
                                  ↓
                            [S3 Output]
```

### Script 1: Splitter (EC2)

**Purpose:** Download input from S3, split into 2 chunks, send to Lambda workers

**Key operations:**
```bash
# Line 98: Download from S3
python3.9 aws/s3-get-object.py "oneliners/inputs/1M.txt" "/tmp/pash_p4qiBDD/.../fifo29"

# Line 96: Split into 2 chunks (round-robin lines)
runtime/r_split -r "/tmp/pash_p4qiBDD/.../fifo2" 1000000 \
  "/tmp/pash_p4qiBDD/.../fifo10" \
  "/tmp/pash_p4qiBDD/.../fifo14"

# Line 100: Send chunks to Lambda workers via holepunch
/opt/pashlib \
  send*49e60090-bea3-43fe-b3a2-039491f8f37c*0*1*/tmp/.../fifo10 \
  send*e64098c5-470e-4972-9803-5b7b13c16ffd*0*1*/tmp/.../fifo14
```

**Connection mapping:**
- `49e60090...` → Lambda Worker 1 (fifo10)
- `e64098c5...` → Lambda Worker 2 (fifo14)

### Script 2: Lambda Worker 1

**Purpose:** Receive chunk from splitter, sort it, send to merger

**Key operations:**
```bash
# Line 95: Receive from splitter, sort, send to merger
/opt/pashlib \
  recv*49e60090-bea3-43fe-b3a2-039491f8f37c*1*0*/tmp/.../fifo12 \
  send*f0082947-1beb-4dda-b5f4-3e82e8a2ef0b*0*1*/tmp/.../fifo18

# Line 91: Sort the data
sort <"/tmp/.../fifo13" >"/tmp/.../fifo18"
```

**Connection mapping:**
- Receives from: `49e60090...` (splitter) → fifo12
- Sends to: `f0082947...` (merger) → fifo18

### Script 3: Lambda Worker 2

**Purpose:** Same as Worker 1, different connection IDs

**Connection mapping:**
- Receives from: `e64098c5...` (splitter) → fifo16
- Sends to: `61d430ca...` (merger) → fifo25

### Script 4: Merger (EC2)

**Purpose:** Receive sorted chunks from both Lambda workers, merge, upload to S3

**Key operations:**
```bash
# Line 101: Receive from both Lambda workers
/opt/pashlib \
  recv*f0082947-1beb-4dda-b5f4-3e82e8a2ef0b*1*0*/tmp/.../fifo20 \
  recv*61d430ca-72a0-4745-bf97-94441b9af0fc*1*0*/tmp/.../fifo27

# Line 95: Merge sorted chunks
sort -m "/tmp/.../fifo21" "/tmp/.../fifo28" >"/tmp/.../fifo22"

# Line 103: Upload to S3
python3.9 aws/s3-put-object.py \
  "oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt" \
  "/tmp/.../fifo22" $1
```

**Connection mapping:**
- Receives from: `f0082947...` (Lambda Worker 1) → fifo20
- Receives from: `61d430ca...` (Lambda Worker 2) → fifo27

---

## How to Run the 4 Scripts

### Prerequisites

**Environment Variables:**
```bash
export PASH_TOP=/home/ubuntu/pash
export AWS_BUCKET=inout741448956691  # Your S3 bucket
export AWS_REGION=us-east-1
```

**Verify pashlib:**
```bash
ls -la /opt/pashlib  # Should exist
aws dynamodb describe-table --table-name rdv --region us-east-1  # Should exist
```

**Verify Lambda function:**
```bash
aws lambda get-function --function-name lambda --region us-east-1
# Should show runtime: python3.9, handler: lambda-function.lambda_handler
```

### Lambda Function Handler

**Location:** `/home/ubuntu/pash/runtime/serverless/lambda-function.py`

**What it does:**
```python
def lambda_handler(event, context):
    # event = {"ids": ["7422cb0d-628b-49ae-b29d-387640f949f6"],
    #          "folder_ids": ["1762478811"]}

    for i, folder_id in enumerate(event['folder_ids']):
        id_ = event['ids'][i]

        # Download script from S3
        key = f"sls-scripts/{folder_id}/{id_}.sh"
        s3.get_object(Bucket=BUCKET, Key=key)

        # Execute script
        subprocess.run(["/bin/bash", f"/tmp/script-{folder_id}-{id_}.sh", folder_id])
```

### Critical Timing Requirements

**The Problem:** Holepunch requires both nodes to attempt connection within ±15 seconds

**Current behavior:**
1. EC2 scripts start immediately → Register in DynamoDB → Wait for Lambda
2. Lambda invoked async → Cold start (0.2-10s) → Register in DynamoDB
3. If timing window exceeded → Connection timeout

**Success requires:**
- Both nodes register in DynamoDB within ~15 seconds
- Both nodes query for peer and find each other
- Both nodes attempt simultaneous TCP connect

---

## Step-by-Step Execution Guide

### Option 1: Using Manual Orchestrator Script

The existing `manual_orchestrate_width2.sh` script does most of the work:

```bash
cd /home/ubuntu/pash
export AWS_BUCKET=inout741448956691

# This will:
# 1. Copy EC2 scripts to /tmp
# 2. Upload Lambda scripts to S3
# 3. Start all 4 nodes in parallel
./manual_orchestrate_width2.sh
```

**What it does:**
```bash
FOLDER_ID=$(date +%s)

# Upload Lambda scripts to S3
aws s3 cp "$SCRIPTS_DIR/${LAMBDA1}.sh" \
  "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA1}.sh"
aws s3 cp "$SCRIPTS_DIR/${LAMBDA2}.sh" \
  "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA2}.sh"

# Start all nodes in parallel (critical for timing!)
/bin/bash "/tmp/$SPLITTER" "$FOLDER_ID" &
python3 "$INVOKE_SCRIPT" "$LAMBDA1" "$FOLDER_ID" "lambda" &
python3 "$INVOKE_SCRIPT" "$LAMBDA2" "$FOLDER_ID" "lambda" &
/bin/bash "/tmp/$MERGER" "$FOLDER_ID" &

wait  # Wait for all to complete
```

### Option 2: Manual Step-by-Step

If you want to understand each step:

```bash
cd /home/ubuntu/pash
export AWS_BUCKET=inout741448956691
export PASH_TOP=/home/ubuntu/pash

SCRIPTS_DIR="/home/ubuntu/pash/evaluation/benchmarks/oneliners/logs/sort.sh:1M.txt:2/scripts/1762478811"
FOLDER_ID=$(date +%s)

# Script IDs
SPLITTER="559c679d-49ed-49c3-9a1f-08c3ce9e5f02"
LAMBDA1="7422cb0d-628b-49ae-b29d-387640f949f6"
LAMBDA2="6b9eae70-3757-42a4-a3fb-eee1dc6a73ca"
MERGER="839ba2b4-ed80-4864-a800-04694e745f76"

# Step 1: Prepare EC2 scripts
cp "$SCRIPTS_DIR/${SPLITTER}.sh" "/tmp/$SPLITTER"
cp "$SCRIPTS_DIR/${MERGER}.sh" "/tmp/$MERGER"
chmod +x "/tmp/$SPLITTER" "/tmp/$MERGER"

# Step 2: Upload Lambda scripts to S3
aws s3 cp "$SCRIPTS_DIR/${LAMBDA1}.sh" \
  "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA1}.sh"
aws s3 cp "$SCRIPTS_DIR/${LAMBDA2}.sh" \
  "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA2}.sh"

# Step 3: Start all nodes in parallel (MUST BE PARALLEL for holepunch timing!)
(cd "$PASH_TOP" && /bin/bash "/tmp/$SPLITTER" "$FOLDER_ID" 2>&1 | tee /tmp/splitter.log) &
SPLITTER_PID=$!

python3 "$PASH_TOP/aws/invoke-lambda.py" "$LAMBDA1" "$FOLDER_ID" "lambda" &
LAMBDA1_PID=$!

python3 "$PASH_TOP/aws/invoke-lambda.py" "$LAMBDA2" "$FOLDER_ID" "lambda" &
LAMBDA2_PID=$!

(cd "$PASH_TOP" && /bin/bash "/tmp/$MERGER" "$FOLDER_ID" 2>&1 | tee /tmp/merger.log) &
MERGER_PID=$!

# Wait for all
wait $SPLITTER_PID $LAMBDA1_PID $LAMBDA2_PID $MERGER_PID
echo "Exit codes: splitter=$?, lambda1=$?, lambda2=$?, merger=$?"
```

### Option 3: Test Individual Components

**Test pashlib connection between EC2 and Lambda:**

```bash
# On EC2: Start a receiver
/opt/pashlib recv*test-connection-1*0*1*/tmp/test-recv.fifo &
echo "Waiting for data..." &
cat /tmp/test-recv.fifo

# In another terminal: Invoke Lambda that sends
# (You'll need to create a test Lambda function that calls pashlib send)
```

---

## Debugging and Troubleshooting

### Check DynamoDB Registrations

```bash
# See recent holepunch attempts
aws dynamodb scan --table-name rdv --limit 10 --region us-east-1

# Check specific connection
RDV_KEY="49e60090-bea3-43fe-b3a2-039491f8f37c"
aws dynamodb get-item --table-name rdv \
  --key "{\"key\": {\"S\": \"$RDV_KEY\"}}" \
  --region us-east-1
```

Expected output:
```json
{
  "Item": {
    "key": {"S": "49e60090-bea3-43fe-b3a2-039491f8f37c"},
    "N_0": {"S": "{\"local_addr\":\"172.31.12.152:54744\",\"external_addr\":\"3.133.71.249:54744\"}"},
    "N_1": {"S": "{\"local_addr\":\"169.254.100.6:33471\",\"external_addr\":\"34.207.116.54:33471\"}"}
  }
}
```

### Check Lambda Logs

```bash
# Get recent Lambda invocations
aws logs tail /aws/lambda/lambda --follow --region us-east-1

# Search for specific execution
aws logs filter-log-events \
  --log-group-name /aws/lambda/lambda \
  --filter-pattern "lambda-function.py" \
  --region us-east-1 \
  --start-time $(date -d '5 minutes ago' +%s)000
```

### Monitor EC2 Script Execution

```bash
# Run with logging
(cd $PASH_TOP && /bin/bash /tmp/559c679d-49ed-49c3-9a1f-08c3ce9e5f02 \
  $(date +%s) 2>&1 | tee /tmp/splitter-debug.log)

# In another terminal, watch the log
tail -f /tmp/splitter-debug.log
```

### Common Issues

#### Issue 1: "Connection timed out"
**Cause:** Lambda cold start took too long, missed the 15-second window

**Solutions:**
- Pre-warm Lambda (invoke it first)
- Increase timeout in pashlib code
- Use provisioned concurrency

#### Issue 2: "Script not found in S3"
**Cause:** Lambda scripts weren't uploaded before invocation

**Solution:**
```bash
# Verify scripts exist in S3
aws s3 ls s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/
```

#### Issue 3: "Peer not found in DynamoDB"
**Cause:** Node didn't register properly

**Solutions:**
- Check STUN service is reachable
- Verify DynamoDB permissions
- Check network connectivity

#### Issue 4: Fifos blocking
**Cause:** Named pipes need both reader and writer

**Solution:** All processes using fifos must be started in parallel (backgrounded)

### Verifying Holepunch Success

Look for these log messages:

**EC2 side:**
```
Get local_addr: 172.31.12.152:54744 external_addr: 3.133.71.249:54744
connecting from 3.133.71.249:54744 to 34.207.116.54:33471
```

**Lambda side (in CloudWatch Logs):**
```
Get local_addr: 169.254.100.6:33471 external_addr: 34.207.116.54:33471
connecting from 34.207.116.54:33471 to 3.133.71.249:54744
```

If you see both, holepunch succeeded!

---

## Next Steps

To improve reliability, consider:

1. **Pre-warm Lambdas** - Invoke them with no-op before real execution
2. **Increase timeouts** - Modify pashlib to wait longer for peer
3. **Add retries** - Retry holepunch connection with exponential backoff
4. **Use S3 orchestration** - Replace holepunch with S3-based communication (see `S3_ORCHESTRATOR_README.md`)
5. **Monitor timing** - Log timestamps at each step to understand delays

---

## Summary

The key to running these 4 scripts is:

1. **All 4 must start in parallel** (within ~15 seconds)
2. **Lambda scripts must be in S3** before Lambda invocation
3. **DynamoDB table "rdv"** must exist and be accessible
4. **/opt/pashlib binary** must exist on both EC2 and Lambda
5. **Connection IDs (UUIDs)** must match between sender and receiver

The existing `manual_orchestrate_width2.sh` script handles most of this, but timing is critical and often fails due to Lambda cold starts.
