# Plan: Fix STUN Timeout Error in Holepunching Library

## Executive Summary

**What's broken:** The unix50 benchmarks crash with `ErrTransactionTimeOut` when STUN server is unreachable
**Why it happens:** Poor error handling - uses `.unwrap()` instead of graceful fallback
**The fix:** 4 targeted Rust code changes across 3 files (42 LOC changed, ~100 LOC added)
**Outcome:** Holepunch continues working when STUN succeeds, falls back gracefully when it fails

## Problem Summary

The unix50 benchmark is failing with a STUN timeout error when running with `--serverless_exec`:

```
thread 'main' panicked at /code/src/stun_helper.rs:30:32:
called `Result::unwrap()` on an `Err` value: ErrTransactionTimeOut
```

**Root Cause:** The holepunching library (`splash-stun-lib`) attempts to query a STUN server at `stun.l.google.com:19302` to discover the external IP address for NAT traversal. When the STUN server is unreachable (network issues, firewall, Lambda NAT restrictions), the query times out after ~8 seconds. The code uses `.unwrap()` on the result, causing a panic instead of handling the error gracefully.

**Impact:** Blocks all serverless benchmarks that rely on pashlib for EC2-Lambda communication.

## Critical Files

1. `/home/ubuntu/splash-stun-lib/src/stun_helper.rs` - STUN query implementation (line 30 panics)
2. `/home/ubuntu/splash-stun-lib/bin/pashlib-oneproc/main.rs` - Main entry point that calls STUN helper (line 94)
3. `/home/ubuntu/splash-stun-lib/src/holepunch.rs` - Holepunch context and environment variable checks (lines 172-197)
4. `/home/ubuntu/pash/evaluation/benchmarks/unix50/run-leash-compare-generic.sh` - Benchmark script that triggers the error

## Recommended Solution: Graceful Error Handling with Fallback

### Approach

Implement a robust fix to the STUN library that handles timeouts gracefully while preserving holepunch functionality for environments where STUN works (EC2-to-EC2 communication):

1. **Handle STUN errors without panicking** - Replace `.unwrap()` with proper error handling
2. **Implement intelligent fallback strategies** - Use AWS metadata service, environment variables, or local addresses when STUN fails
3. **Preserve holepunch for working cases** - Keep STUN-based NAT traversal functional for EC2 environments where it succeeds
4. **Add configuration options** - Allow timeout tuning and fallback behavior via environment variables

### Implementation Steps

#### Step 1: Fix STUN Helper Error Handling

**File:** `splash-stun-lib/src/stun_helper.rs`

**Current code (lines 9-40):**
```rust
pub async fn get_addr() -> (String, String) {
    let server = "stun.l.google.com:19302";
    // ... setup code ...
    let event = handler_rx.recv().await.unwrap();
    let msg = event.event_body.unwrap();  // LINE 30 - PANICS HERE
    assert_eq!(msg.typ, BINDING_SUCCESS);
    // ... process response ...
    (local_addr, addr)
}
```

**Changes needed:**

1. **Change function signature** (line 9):
   ```rust
   pub async fn get_addr() -> Result<(String, String), String>
   ```

2. **Replace panic with error handling** (lines 29-40):
   ```rust
   let event = handler_rx.recv().await.unwrap();

   // Handle timeout or other STUN errors
   let msg = match event.event_body {
       Ok(msg) => msg,
       Err(err) => {
           eprintln!("STUN query failed: {:?}", err);
           return Err(format!("STUN timeout: {:?}", err));
       }
   };

   // Validate response
   if msg.typ != BINDING_SUCCESS {
       return Err(format!("STUN response type mismatch: expected BINDING_SUCCESS, got {:?}", msg.typ));
   }
   if msg.transaction_id != txn_id {
       return Err("STUN transaction ID mismatch".to_string());
   }

   let mut xor_addr = XorMappedAddress::default();
   if let Err(e) = xor_addr.get_from(&msg) {
       return Err(format!("Failed to parse XOR mapped address: {:?}", e));
   }

   let addr = xor_addr.to_string();
   client.close().await.unwrap();

   Ok((local_addr, addr))
   ```

3. **Replace assertions with proper validation** - Remove `assert_eq!` statements and return errors instead.

#### Step 2: Update Caller to Handle STUN Failures

**File:** `splash-stun-lib/bin/pashlib-oneproc/main.rs`

**Add import at top of file:**
```rust
use pashlib::holepunch::resolve_public_addr_fallback;  // New function we'll create
```

**Current code (lines 91-100):**
```rust
let (local_addr_global, external_addr_global) = if no_holepunch_enabled() {
    (String::new(), String::new())
} else {
    let (local_addr_global, external_addr_global) = stun_helper::get_addr().await;
    eprintln!("Get local_addr: {} external_addr: {}", local_addr_global, external_addr_global);
    (local_addr_global, external_addr_global)
};
```

**Replace with:**
```rust
let (local_addr_global, external_addr_global) = if no_holepunch_enabled() {
    (String::new(), String::new())
} else {
    match stun_helper::get_addr().await {
        Ok((local_addr, external_addr)) => {
            eprintln!("✓ STUN discovery succeeded - local_addr: {} external_addr: {}", local_addr, external_addr);
            (local_addr, external_addr)
        },
        Err(e) => {
            eprintln!("⚠ STUN discovery failed: {}", e);
            eprintln!("  Using fallback address resolution...");

            // Fallback: bind to get local address, then resolve public address
            let (local_addr, external_addr) = resolve_public_addr_fallback().await;
            eprintln!("✓ Fallback succeeded - local_addr: {} external_addr: {}", local_addr, external_addr);

            (local_addr, external_addr)
        }
    }
};
```

#### Step 3: Implement Fallback Address Resolution

**File:** `splash-stun-lib/src/holepunch.rs`

**Add new public function** (insert after line 201, after the `connect_no_holepunch` function):

```rust
/// Fallback address resolution when STUN fails
/// Binds a UDP socket to get local address, then uses resolve_public_addr for external address
pub async fn resolve_public_addr_fallback() -> (String, String) {
    use tokio::net::UdpSocket;

    // Bind UDP socket to get local address (same as STUN helper does)
    let local_addr = match UdpSocket::bind("0.0.0.0:0").await {
        Ok(socket) => {
            match socket.local_addr() {
                Ok(addr) => {
                    eprintln!("  Bound local address: {}", addr);
                    addr
                },
                Err(e) => {
                    eprintln!("  Warning: Failed to get local_addr from socket: {}", e);
                    "0.0.0.0:0".parse().unwrap()
                }
            }
        },
        Err(e) => {
            eprintln!("  Warning: Failed to bind UDP socket: {}", e);
            "0.0.0.0:0".parse().unwrap()
        }
    };

    // Use existing resolve_public_addr function to determine external address
    // This function already checks PASH_PUBLIC_ADDR and PASH_PUBLIC_IP environment variables
    let external_addr = resolve_public_addr(local_addr);

    // If still using local address after resolve_public_addr, try AWS metadata
    let external_addr = if external_addr == local_addr.to_string() {
        match query_aws_metadata_public_ip().await {
            Ok(public_ip) => {
                eprintln!("  AWS metadata service returned public IP: {}", public_ip);
                format!("{}:{}", public_ip, local_addr.port())
            },
            Err(e) => {
                eprintln!("  AWS metadata query failed: {}", e);
                external_addr
            }
        }
    } else {
        external_addr
    };

    (local_addr.to_string(), external_addr)
}

/// Query AWS EC2 metadata service for public IPv4 address
/// Returns the public IP if successful, error otherwise
async fn query_aws_metadata_public_ip() -> Result<String, String> {
    use tokio::net::TcpStream;
    use tokio::io::{AsyncReadExt, AsyncWriteExt};
    use std::time::Duration;

    // AWS metadata service endpoint
    let metadata_endpoint = "169.254.169.254:80";

    // Short timeout - metadata service should respond quickly
    let socket = tokio::time::timeout(
        Duration::from_secs(2),
        TcpStream::connect(metadata_endpoint)
    )
    .await
    .map_err(|_| "Timeout connecting to AWS metadata service".to_string())?
    .map_err(|e| format!("Failed to connect to AWS metadata service: {}", e))?;

    let mut stream = socket;

    // HTTP request for public IPv4
    let request = b"GET /latest/meta-data/public-ipv4 HTTP/1.0\r\nHost: 169.254.169.254\r\n\r\n";

    stream
        .write_all(request)
        .await
        .map_err(|e| format!("Failed to write request: {}", e))?;

    // Read response
    let mut response = String::new();
    stream
        .read_to_string(&mut response)
        .await
        .map_err(|e| format!("Failed to read response: {}", e))?;

    // Parse HTTP response - body comes after "\r\n\r\n"
    if let Some(body) = response.split("\r\n\r\n").nth(1) {
        let ip = body.trim();
        // Validate it looks like an IP address
        if ip.split('.').count() == 4 && ip.len() >= 7 {
            return Ok(ip.to_string());
        }
    }

    Err("Failed to parse public IP from metadata response".to_string())
}
```

**Notes:**
- The existing `resolve_public_addr` function (line 221-236) already checks `PASH_PUBLIC_ADDR` and `PASH_PUBLIC_IP` environment variables
- We're reusing that logic and adding AWS metadata as an additional fallback
- The function is async to match the Tokio runtime used by pashlib

#### Step 4: (Optional) Make Pashlib Export Public Function

**File:** `splash-stun-lib/src/lib.rs`

Ensure the new `resolve_public_addr_fallback` function is exported in the public API:

```rust
pub mod holepunch;
pub mod stun_helper;

// Re-export commonly used functions
pub use holepunch::{resolve_public_addr_fallback, no_holepunch_enabled};
```

This allows the function to be imported in `main.rs` with:
```rust
use pashlib::holepunch::resolve_public_addr_fallback;
```

## Testing During Development (Using Environment Variables)

While implementing the Rust fixes, you can test the benchmarks using these environment variable workarounds:

### Option A: Disable Holepunch Entirely
```bash
export PASH_NO_HOLEPUNCH=1
./run-leash-compare-generic.sh --small
```

### Option B: Manually Specify Public Address
```bash
# Get your EC2 public IP
PUBLIC_IP=$(curl -s http://169.254.169.254/latest/meta-data/public-ipv4)
export PASH_PUBLIC_IP="$PUBLIC_IP"
./run-leash-compare-generic.sh --small
```

### Option C: Reduce STUN Timeout
```bash
export PASH_STUN_TIMEOUT_MS=3000  # 3 seconds instead of 8
./run-leash-compare-generic.sh --small
```

## Testing & Verification

After implementing the fix:

1. **Test STUN failure handling:**
   ```bash
   # Block STUN server to force timeout
   sudo iptables -A OUTPUT -d 64.233.160.0/19 -j DROP  # Block Google IPs
   ./run-leash-compare-generic.sh --small
   sudo iptables -D OUTPUT -d 64.233.160.0/19 -j DROP  # Cleanup
   ```

2. **Test with environment variable fallback:**
   ```bash
   export PASH_PUBLIC_IP="3.133.71.249"
   ./run-leash-compare-generic.sh --small
   ```

3. **Test with holepunch disabled:**
   ```bash
   export PASH_NO_HOLEPUNCH=1
   ./run-leash-compare-generic.sh --small
   ```

4. **Verify logs show fallback messages:**
   - Check for "STUN discovery failed" message
   - Check for "Fallback - local_addr: ... external_addr: ..." message
   - Confirm no panic/crash occurs

5. **Validate benchmark results:**
   - Check that noopt and s3opt outputs match
   - Verify logs are generated correctly
   - Confirm final summary shows "✓ MATCH"

## Build & Deployment

After implementing the Rust code changes:

### 1. Build the Modified Pashlib Binary

```bash
cd /home/ubuntu/splash-stun-lib

# Clean previous builds
cargo clean

# Build release version
cargo build --release --bin pashlib-oneproc

# Verify the binary was built
ls -lh target/release/pashlib-oneproc
```

### 2. Deploy to Local Runtime

```bash
# Copy to pash runtime directory (for local EC2 execution)
cp target/release/pashlib-oneproc /home/ubuntu/pash/runtime/serverless/runtime/pashlib

# Make executable
chmod +x /home/ubuntu/pash/runtime/serverless/runtime/pashlib

# Verify
file /home/ubuntu/pash/runtime/serverless/runtime/pashlib
```

### 3. Update Lambda Layer (If Using Lambda)

If you're running benchmarks that use AWS Lambda:

```bash
cd /home/ubuntu/pash/runtime/serverless

# Copy new binary to Lambda runtime
cp /home/ubuntu/splash-stun-lib/target/release/pashlib-oneproc runtime/pashlib

# Rebuild Lambda layer zip
cd runtime
zip -r ../lambda-layer.zip .

# Upload to AWS Lambda (update layer name as needed)
aws lambda publish-layer-version \
    --layer-name pash-runtime \
    --zip-file fileb://../lambda-layer.zip \
    --compatible-runtimes python3.9 python3.10

# Update Lambda function to use new layer version
# (This step depends on your Lambda function name)
```

### 4. Verify Installation

```bash
# Test that pashlib runs without crashing
/home/ubuntu/pash/runtime/serverless/runtime/pashlib --help

# Or run a simple send/recv test with holepunch disabled
export PASH_NO_HOLEPUNCH=1
/home/ubuntu/pash/runtime/serverless/runtime/pashlib send*test*0*1*/tmp/test.fifo
```

## Benefits of This Approach

1. **Non-breaking:** Existing working setups continue to work with STUN
2. **Graceful degradation:** Falls back to alternative methods when STUN fails
3. **Configurable:** Environment variables allow easy testing and debugging
4. **AWS-aware:** Automatically queries EC2 metadata for public IP
5. **Better error messages:** Clear indication when STUN fails and why

## Notes

- The STUN timeout is caused by network restrictions in AWS Lambda/EC2 environments
- The architecture documentation (PASH_SERVERLESS_ARCHITECTURE.md) notes that holepunch has 0% success rate in Lambda
- The S3 direct streaming optimization (`--enable_s3_direct`) is designed to bypass holepunch entirely
- Consider longer-term: deprecate holepunch in favor of S3-based orchestration
