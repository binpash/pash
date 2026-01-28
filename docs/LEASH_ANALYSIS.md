# Leash Computation Graph Analysis - Sort with Width 2 & 4

## Overview

This document analyzes how Leash (PaSh serverless) splits and merges computation through EC2 and Lambda for the `sort` operation at different parallelization widths.

---

## Width=4 Computation Graph

**Location**: `/home/ubuntu/pash/evaluation/benchmarks/oneliners/logs/sort.sh:1M.txt:4/scripts/1762465916/`

### Graph Structure

```
                    [EC2: Input Splitter]
                       (253b2697)
                  Downloads 1M.txt from S3
                  Splits into 4 parts (250K lines each)
                           |
        +--------+--------+--------+--------+
        |        |        |        |        |
    [Lambda 1] [Lambda 2] [Lambda 3] [Lambda 4]
    (27a22ad8) (d21c33f5) (26c6892a) (6db40e5f)
      sort       sort       sort       sort
        |        |        |        |        |
        +----+---+        +----+---+        |
             |                 |            |
     [EC2: Merger L1-L]  [EC2: Merger L1-R]
        (80a60533)         (db76fd71)
         sort -m             sort -m
             |                 |
             +--------+--------+
                      |
            [EC2: Final Merger]
               (0de9f449)
                sort -m
                    |
              Upload to S3
```

### Key Scripts (Width=4)

#### Lambda 1 (27a22ad8) - Parallel Sort Worker
```bash
{ sort <"#fifo23" >"#fifo36" & }
{ runtime/dgsh-tee -i "#fifo22" -o "#fifo23" -I -b 5M & }
{ /opt/pashlib recv*eb04493e-ed53-4671-8fcd-5146987cf675*1*0*#fifo22 \
                send*b6f1a3ac-dd6a-4e46-8bee-1711ecd91f45*0*1*#fifo36 & }
```
**Flow**: Receive data via network → buffer → sort → send sorted data back

#### EC2 Merger L1-L (80a60533) - First Level Merge
```bash
{ sort -m "#fifo51" "#fifo58" >"#fifo52" & }
{ runtime/dgsh-tee -i "#fifo50" -o "#fifo51" -I -b 5M & }
{ runtime/dgsh-tee -i "#fifo57" -o "#fifo58" -I -b 5M & }
{ /opt/pashlib recv*1213986b-cd5a-4d3a-b6e1-94cad699a1a8*1*0*#fifo50 \
                recv*e3a079b5-8265-4d08-bfc0-606c885d6745*1*0*#fifo57 & }
```
**Flow**: Receive from Lambda 1 & 2 → buffer → `sort -m` merges

#### EC2 Final Merger (0de9f449) - Output to S3
```bash
{ sort -m "#fifo51" "#fifo58" >"#fifo52" & }
{ /opt/pashlib recv*... }
{ python3.9 aws/s3-put-object.py oneliners/outputs/sort.sh:1M.txt:4:leashstdout.txt \
                                 "#fifo52" $1 & }
```
**Flow**: Receive from both L1 mergers → merge → stream to S3

---

## Width=2 Computation Graph

**Location**: `/home/ubuntu/pash/evaluation/benchmarks/oneliners/logs/sort.sh:1M.txt:2/scripts/1762478811/`

### Graph Structure

```
              [EC2: Input Splitter]
                 (559c679d)
            Downloads 1M.txt from S3
            Splits into 2 parts (500K lines each)
                     |
              +------+------+
              |             |
          [Lambda 1]    [Lambda 2]
          (7422cb0d)    (6b9eae70)
            sort          sort
              |             |
              +------+------+
                     |
            [EC2: Final Merger]
                (839ba2b4)
                 sort -m
                     |
              Upload to S3
```

**Key Difference from Width=4**: No intermediate merge layer! Just split → parallel sort → final merge.

### Script IDs (Width=2)
- **Splitter**: `559c679d-49ed-49c3-9a1f-08c3ce9e5f02`
- **Lambda 1**: `7422cb0d-628b-49ae-b29d-387640f949f6`
- **Lambda 2**: `6b9eae70-3757-42a4-a3fb-eee1dc6a73ca`
- **Merger**: `839ba2b4-ed80-4864-a800-04694e745f76`

### Key Scripts (Width=2)

#### Lambda 1 (7422cb0d) - Parallel Sort Worker
```bash
{ sort <"#fifo13" >"#fifo18" & }
{ runtime/dgsh-tee -i "#fifo12" -o "#fifo13" -I -b 5M & }
{ /opt/pashlib recv*49e60090-bea3-43fe-b3a2-039491f8f37c*1*0*#fifo12 \
                send*f0082947-1beb-4dda-b5f4-3e82e8a2ef0b*0*1*#fifo18 & }
```

#### EC2 Splitter (559c679d) - Download & Split
```bash
{ cat "#fifo29" >"#fifo2" & }
{ runtime/r_split -r "#fifo2" 1000000 "#fifo10" "#fifo14" & }
{ python3.9 aws/s3-get-object.py "oneliners/inputs/1M.txt" "#fifo29" & }
{ /opt/pashlib send*49e60090-bea3-43fe-b3a2-039491f8f37c*0*1*#fifo10 \
                send*e64098c5-470e-4972-9803-5b7b13c16ffd*0*1*#fifo14 & }
```
**Flow**: S3 download → r_split divides at 1M lines into 2 chunks → send to Lambdas

#### EC2 Merger (839ba2b4) - Merge & Upload
```bash
{ sort -m "#fifo21" "#fifo28" >"#fifo22" & }
{ runtime/dgsh-tee -i "#fifo20" -o "#fifo21" -I -b 5M & }
{ runtime/dgsh-tee -i "#fifo27" -o "#fifo28" -I -b 5M & }
{ /opt/pashlib recv*f0082947-1beb-4dda-b5f4-3e82e8a2ef0b*1*0*#fifo20 \
                recv*61d430ca-72a0-4745-bf97-94441b9af0fc*1*0*#fifo27 & }
{ python3.9 aws/s3-put-object.py oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt \
                                 "#fifo22" $1 & }
```

---

## Network Communication

### pashlib Protocol

The `/opt/pashlib` binary creates point-to-point TCP connections using unique IDs:

```
Splitter → Lambda 1:  send*49e60090... → recv*49e60090...
Lambda 1 → Merger:    send*f0082947... → recv*f0082947...

Splitter → Lambda 2:  send*e64098c5... → recv*e64098c5...
Lambda 2 → Merger:    send*61d430ca... → recv*61d430ca...
```

### Connection Format
- **Send**: `send*{uuid}*{port}*{mode}*{fifo_path}`
- **Recv**: `recv*{uuid}*{port}*{mode}*{fifo_path}`

Matching UUIDs establish connections between nodes.

---

## S3 Download/Upload

### s3-get-object.py
- **Location**: `/home/ubuntu/pash/aws/s3-get-object.py`
- **Usage**: `python3.9 aws/s3-get-object.py "oneliners/inputs/1M.txt" "/tmp/fifo"`
- **Mechanism**: Streams S3 object directly to a FIFO pipe (no disk storage)

### s3-put-object.py
- **Location**: `/home/ubuntu/pash/aws/s3-put-object.py`
- **Usage**: `python3.9 aws/s3-put-object.py "output/key" "/tmp/fifo"`
- **Mechanism**: Reads from FIFO pipe and streams to S3 (no intermediate files)

**This is pure streaming I/O** - data flows through memory (FIFOs) without hitting disk!

---

## Orchestration Infrastructure

### EC2 Handler
- **Location**: `/home/ubuntu/pash/scripts/serverless/ec2-handler.py`
- **Port**: 9999
- **Function**: Listens for script execution requests via TCP socket
- **Message Format**: `{"ids": [script_id], "folder_ids": [folder_id]}`
- **Process**:
  1. Receives JSON request
  2. Downloads script from S3: `s3://{bucket}/sls-scripts/{folder_id}/{script_id}.sh`
  3. Executes: `/bin/bash /tmp/script-{folder_id}-{script_id}.sh {folder_id}`

### Lambda Function
- **Location**: `/home/ubuntu/pash/runtime/serverless/lambda-function.py`
- **Function Name**: `lambda`
- **Invocation**: Event (async)
- **Payload**: `{"ids": [script_id], "folder_ids": [folder_id]}`
- **Process**: Same as EC2 handler but runs in Lambda environment

### Invoke Script
- **Location**: `/home/ubuntu/pash/aws/invoke-lambda.py`
- **Usage**: `python3 invoke-lambda.py {script_id} {folder_id} {type}`
  - `type`: "lambda" or "ec2"
- **EC2_IP**: localhost (for local EC2 handler)
- **EC2_PORT**: 9999

**Fixed Issues**:
- ✅ Changed EC2 message format to use arrays: `{"ids": [...], "folder_ids": [...]}`
- ✅ Set EC2_IP to `localhost`

---

## Manual Orchestration Script

### File: `/home/ubuntu/pash/manual_orchestrate_width2.sh`

This script runs the **exact 4 pre-generated scripts** without recompiling through PaSh.

```bash
#!/bin/bash
# Manual orchestration for width=2 sort using pre-generated scripts

SCRIPTS_DIR="/home/ubuntu/pash/evaluation/benchmarks/oneliners/logs/sort.sh:1M.txt:2/scripts/1762478811"
PASH_TOP="/home/ubuntu/pash"
INVOKE_SCRIPT="$PASH_TOP/aws/invoke-lambda.py"

SPLITTER="559c679d-49ed-49c3-9a1f-08c3ce9e5f02"
LAMBDA1="7422cb0d-628b-49ae-b29d-387640f949f6"
LAMBDA2="6b9eae70-3757-42a4-a3fb-eee1dc6a73ca"
MERGER="839ba2b4-ed80-4864-a800-04694e745f76"

FOLDER_ID=$(date +%s)

# Step 1: Upload scripts to S3
for script_id in $SPLITTER $LAMBDA1 $LAMBDA2 $MERGER; do
    aws s3 cp "$SCRIPTS_DIR/${script_id}.sh" \
        "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${script_id}.sh"
done

# Step 2: Start splitter on EC2
python3 "$INVOKE_SCRIPT" "$SPLITTER" "$FOLDER_ID" "ec2"
sleep 2

# Step 3: Start Lambda workers (parallel)
python3 "$INVOKE_SCRIPT" "$LAMBDA1" "$FOLDER_ID" "lambda" &
python3 "$INVOKE_SCRIPT" "$LAMBDA2" "$FOLDER_ID" "lambda" &
wait
sleep 2

# Step 4: Start merger on EC2
python3 "$INVOKE_SCRIPT" "$MERGER" "$FOLDER_ID" "ec2"
```

### How to Use

**Terminal 1** - Start EC2 Handler:
```bash
cd /home/ubuntu/pash
python3 scripts/serverless/ec2-handler.py
```

**Terminal 2** - Run Orchestration:
```bash
cd /home/ubuntu/pash
./manual_orchestrate_width2.sh
```

### Prerequisites
- `AWS_BUCKET` environment variable set
- `AWS_ACCOUNT_ID` environment variable set (for Lambda)
- EC2 handler running on port 9999
- Lambda function deployed and named "lambda"
- Input file at: `s3://{bucket}/oneliners/inputs/1M.txt`

### Output
- **S3 Path**: `s3://{bucket}/oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt`
- **Lambda Logs**: CloudWatch `/aws/lambda/lambda`
- **EC2 Logs**: Console output from ec2-handler.py

---

## Key Insights

1. **Binary Merge Tree**: Width=4 uses a 2-level merge tree; Width=2 is flat (no intermediate merges)

2. **Streaming Architecture**: No intermediate disk files - all data flows through FIFOs and network

3. **Network IDs**: UUID-based send/recv pairs establish point-to-point connections via pashlib

4. **S3 Integration**: Direct streaming to/from S3 via Python boto3 scripts

5. **Orchestration**: Scripts are stored in S3, workers pull and execute them

6. **Execution Model**:
   - **EC2**: Long-running handler listens on port 9999
   - **Lambda**: Event-driven invocation via boto3
   - Both execute bash scripts from S3

---

## Files Modified

### `/home/ubuntu/pash/aws/invoke-lambda.py`
- Changed `EC2_IP = "localhost"`
- Fixed EC2 message format: `{"ids": [id_], "folder_ids": [folder_id]}`

### Created Files
- `/home/ubuntu/pash/manual_orchestrate_width2.sh` - Manual orchestration script
- `/home/ubuntu/pash/docs/LEASH_ANALYSIS.md` - This document

---

## Next Steps

To extend this analysis:

1. **Width=8**: Analyze 3-level merge tree structure
2. **Different Commands**: Compare computation graphs for `wf.sh`, `spell.sh`, etc.
3. **Performance Analysis**: Measure Lambda execution times, network transfer speeds
4. **Cost Optimization**: Analyze Lambda billed duration vs EC2 compute
5. **Failure Handling**: Investigate retry mechanisms and error recovery

---

## References

- **Original Script**: `/home/ubuntu/pash/evaluation/benchmarks/oneliners/scripts/sort.sh`
- **Run Script**: `/home/ubuntu/pash/evaluation/benchmarks/oneliners/run-leash.sh`
- **Serverless Executor**: `/home/ubuntu/pash/compiler/serverless/serverless_executor.py`
- **Analysis Utils**: `/home/ubuntu/pash/scripts/serverless/utils.py`
