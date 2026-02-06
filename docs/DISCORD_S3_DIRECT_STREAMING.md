# Discord Message: S3 Direct Streaming Findings

Hey everyone, wanted to share what I found trying to optimize the S3 -> Lambda data path. Spent a while on this and hit some walls, so looking for ideas.
Good new is that it's correct now (on all the scripts I've tested - 10 - and for inputs up to 20GB) (it encounters issues very frequently with the approximate boundaries method outlined below whereas the smart boundaries and original leash don't) but performance is suffering.

---

## The Problem I Was Trying to Solve

Right now when we process a big file (say 20GB), the data flow looks like this:

```
S3 --> EC2 (download entire file) --> r_split --> pashlib send --> 16 Lambdas
```

That EC2 hop felt wasteful. The file is already in S3, why download it to EC2 just to split and re-upload chunks to Lambda? Why not have each Lambda grab its slice directly from S3?

```
S3 --> Lambda 0 (bytes 0-1.25GB)
   --> Lambda 1 (bytes 1.25-2.5GB)
   --> ...
   --> Lambda 15 (bytes 18.75-20GB)
```

---

## The Line Boundary Problem

You can't just divide a text file by bytes. If you tell Lambda 1 "read bytes 1250000000-2500000000", you're probably cutting a line in half. Lambda 0 has the beginning of that line, Lambda 1 has the end. Your output is garbage.

So I tried two approaches to fix this:

---

## Approach 1: Lambda-to-Lambda Communication (smart distributed reader based on DFR in DISH)

The idea: each Lambda downloads its approximate byte range, then talks to its neighbors to figure out where lines actually start/end.

Code is in `aws/s3-chunk-reader-approx-tail-coordination.py`

```
Lambda 0: downloads bytes 0-1.25GB, processes, tells Lambda 1 "my last line ended at byte X"
Lambda 1: downloads bytes 1.25-2.5GB, waits to hear from Lambda 0, skips partial first line
Lambda 2: same thing, waits for Lambda 1
...
```

**The problem:** This creates a sequential dependency chain. Lambda 1 can't finalize until Lambda 0 tells it where the boundary is. Lambda 2 waits for Lambda 1. And so on.

Each handshake involves S3 marker polling (~200ms per check) (so we don't open servers we don't need to) plus pashlib TCP. With 16 Lambdas on a 20GB file, these delays cascade and we hit FIFO timeouts a lot of the time or even lambda timeouts if it's too slow.

For small files it works, but it's still slower than the original because of all the coordination overhead.

There might be another issue which is causing this to timeout but I'm not sure.
Also this fails in many cases but I'm not sure why. I still have to debug it more.
---

## Approach 2: Smart Boundaries (EC2 Pre-Calculates)

What if EC2 does a quick scan upfront to find the exact line boundaries, then tells each Lambda "read bytes X to Y" where X and Y are guaranteed to be line-aligned?

Code is in `compiler/serverless/ir_helper.py:815` - the `find_line_boundaries_smart()` function.

```python
# For each boundary between shards
approx_boundary = i * (file_size // num_shards)
# Download a 1MB window around that point
response = s3.get_object(Range=f'bytes={start}-{end}')
# Find the first newline after our approximate boundary
newline_pos = chunk.index(b'\n', offset)
actual_boundary = start + newline_pos + 1
```

For a 20GB file with 16 shards, we need 15 boundaries (1 MB each). That's about 15MB of S3 reads - 0.075% of the file. Takes 2-3 seconds. Pretty cheap.

**This works!** No inter-Lambda communication needed. Each Lambda gets its exact byte range and just processes it.

But it's still slower than the original in many cases.

### Smart Boundary Scanning with Adaptive Window Sizing

The smart boundaries approach now uses **multi-point adaptive window sizing** to minimize data transfer while handling variable line lengths:

#### How It Works

1. **Multi-Point Sampling (automatic):**
   - Takes 5 samples of 256KB each from different file positions
   - Positions: beginning (0%), quarter (25%), middle (50%), three-quarter (75%), end (100%)
   - Total sampling: 1.25MB across 5 S3 requests (~600ms)
   - Calculates average line length across all samples

2. **Adaptive Window Calculation:**
   - Sets window to hold ~500 lines based on detected average
   - Clamped to 4KB minimum, 1MB maximum
   - Provides safety margin for line length variance

3. **Benefits:**
   - **10-20 byte lines:** Uses 5-10KB windows (100× less data than 1MB)
   - **100 byte lines:** Uses 50KB windows (20× less data)
   - **Variable lines:** Multi-point sampling captures variance across file
   - **Very long lines:** Falls back to 1MB window (current behavior)

#### Performance Impact

**For 10-20 byte lines (16 shards example):**
- Sampling: 1.25MB, ~600ms (5 requests)
- Boundary scan: 120KB, ~1.8s (15 requests)
- **Total: 1.37MB, ~2.4s** (vs previous 15MB, ~3s)
- **Improvement: 11× less data, 20% faster**

**For variable line files:**
- Multi-point sampling provides more accurate estimate
- Reduces risk of under-sized windows missing long lines

#### Manual Override (optional)

**Environment Variable:**
- `PASH_BOUNDARY_WINDOW_KB`: Force specific window size in KB
  - Example: `export PASH_BOUNDARY_WINDOW_KB=16` (use 16KB windows)
  - Skips multi-point sampling when set
  - Use when line length is known or for testing

**When to Override:**
- Known uniform line length → set explicit window
- Debugging → force large window to rule out boundary issues
- Very small files → skip sampling overhead

---

## The Numbers

For the covid benchmark it seems quicker but for unix50/6.sh we have:
Script	Input	Mode 1 (No Opt)	Mode 2 (S3 Opt - Smart Boundaries)	Mode 2 (S3 Opt - Non-Smart Boundaries)
6.sh	3_1G.txt	18.312s	         31.949s	                           28.384s
9.sh	4_1G.txt	14.860s	         10.799s	                           11.652s
10.sh	4_1G.txt	17.651s	         9.478s	                           9.253s
19.sh	8_1G.txt	17.451s	         24.415s	                           21.376s


See the below results for 6.sh:

| Input Size | Original (EC2 splitter)| Smart Boundaries | Lambda-to-Lambda* |
|------------|------------------------|------------------|------------------|
| 1MB        | 4-7s                    |  2.5-8s             | 30s              |
| 1GB        | 20-47s                  |   ~30s             |  28s             |
| 20GB       | 3-4 min                | 5.5-6 mins          | a lot of TIMEOUTs, sometimes 6 mins     |

*: a lot of overhead for starting requester opt/pashlib services and polling

So depends on the script a lot but for 20GB it's almost always slower (50-60%) to do s3 direct opt or it even times out.
But for smaller to medium input sizes it is faster at least for smart boundaries. 

---

## What's Going Wrong?

### Issue 1: FIFO Backpressure

I kept seeing this error in the Lambda logs:

```
[STREAM] Shard range: bytes=15728460011-17039165010 (1249.99 MB)
[STREAM] Chunk 1: 256.0 MB / 1250.0 MB (20.5%) at 0.7 MB/s
[STREAM] ERROR: 273256276 read, but total bytes expected is 1310705000
botocore.exceptions.IncompleteReadError
```

What's happening: Lambda is streaming S3 data into a FIFO. The downstream shell command (sort, cut, whatever) reads from that FIFO. But if the shell is slow, the FIFO buffer (64KB by default) fills up. `fifo.write()` blocks. Meanwhile, boto3's S3 connection has a read timeout (~60s). If we're blocked too long, the S3 connection dies and throws IncompleteReadError.

So basically: slow downstream consumer -> FIFO fills -> S3 times out -> job fails.


---

## What I Already Tried

- **Streaming directly to FIFO** using s3.getobj and shutil.copy into the output fifo instead instead of buffering in /tmp - helped a lot saved around 1-2 minutes
- **Larger chunks** (256MB instead of 8MB) - fewer boto3 calls, didn't change much
- **Rewrote in Go** - same performance, so it's not Python being slow, it's network-bound
- **markers in a DynamoDB table for servers** instead of creating many opt/pashlib servers constantly waiting for a possible connection, lambda i instead polls this DDB table for markers from any of the lambdas < i (they put a marker when they want to connect to lambda i). This is very similar to the original implementation of the DFR in DISH. This leads to faster connections between lambdas and less overhead, but doesn't fix the fundamental issue

---

## Questions

1. **Is there a heuristic we should use?** Like "use S3 direct for files under X GB, use EC2 splitter for bigger files"? Or is the 50-60% slowdown just not worth it ever?

2. **Why is smart boundaries slower?** My best guesses:
   - S3 rate limiting when 16 Lambdas download simultaneously
   - Network bandwidth saturated
   - Something in the merge phase on EC2?

   Anyone have ideas?

3. **The backpressure problem** - is there a better way to handle this? Bigger FIFO buffer? Pre-buffer some data to Lambda /tmp? Accept that some jobs will fail?

4. Maybe it's because these lambdas only have 3GB of memory vs a large EC2 instance with 32GB. I tried extending to 10GB but I can't get it to work.

5. What numbers should I collect/ what experiments should I run? What are the next steps here?

---

## Files for Reference

- `compiler/serverless/ir_helper.py`
   - l. 815: def find_line_boundaries_smart() - smart boundary calculation
   - l. 968:  main optimization entry point calls def optimize_s3_lambda_direct_streaming() which does the s3 optimization
- `aws/s3-chunk-reader-approx-tail-coordination.py` - lambda-to-lambda approx. boundaries approach
- `aws/s3-chunk-reader-smart-prealigned.py` - smart boundaries approach
- `compiler/definitions/ir/nodes/serverless_remote_pipe.py:44-76` - chooses whether to do 'no opt', 's3 opt with smart boundaries', 's3 opt with non-smart (approx.) boundaries'

Here is the git branch I'm working on (just committed): 
s3_direct_lambda
[text](https://github.com/binpash/pash/tree/s3_direct_lambda)
