# STUN Timeout During Holepunch

Canonical guide for diagnosing and fixing `ErrTransactionTimeOut` in holepunch (`pashlib`).

## Incident Summary

The failure appears as a panic during STUN discovery when running serverless LEASH workloads.

Typical panic:

```text
thread 'main' panicked at /code/src/stun_helper.rs:30:32:
called `Result::unwrap()` on an `Err` value: ErrTransactionTimeOut
```

Impact:

- Holepunch transport crashes before TCP data flow begins.
- Serverless benchmark runs fail when they depend on STUN reachability.

## Symptoms

Example log sequence:

```text
Get local_addr: 172.31.12.152:46604 external_addr: 3.133.71.249:46604
thread 'main' panicked at /code/src/stun_helper.rs:30:32:
called `Result::unwrap()` on an `Err` value: ErrTransactionTimeOut
```

The panic happens in STUN discovery, before send/recv holepunch connections are established.

## Root Cause

`stun_helper.rs` unwraps STUN event results directly. On timeout, `event.event_body` is an error, so `unwrap()` panics.

Problematic pattern:

```rust
let event = handler_rx.recv().await.unwrap();
let msg = event.event_body.unwrap();
```

Common underlying causes:

- No internet/NAT path from EC2 private subnet.
- Blocked UDP egress to STUN (`19302`), or restrictive NACL/SG rules.
- DNS/packet-loss/transient STUN server issues.

## Affected Code Paths

Primary files:

1. `/home/ubuntu/splash-stun-lib/src/stun_helper.rs`
2. `/home/ubuntu/splash-stun-lib/bin/pashlib-oneproc/main.rs`
3. `/home/ubuntu/splash-stun-lib/src/holepunch.rs`
4. `/home/ubuntu/pash/runtime/serverless/runtime/pashlib` (deployed binary)

Useful benchmark entrypoint for repro:

- `evaluation/benchmarks/unix50/run-leash-benchmark-matrix.sh`

## Quick Triage

### 1) Check STUN reachability

```bash
cd /home/ubuntu/splash-stun-lib
cargo run --bin stun-test
```

If this fails, fix network egress first.

### 2) Validate EC2 egress requirements

- Route to internet (IGW for public subnet, NAT gateway for private subnet).
- Security group egress allows UDP to port `19302`.
- NACL allows outbound UDP and return traffic.

### 3) Temporary workarounds

Disable holepunch:

```bash
export PASH_NO_HOLEPUNCH=1
```

Provide known public IP (if your build supports it):

```bash
PUBLIC_IP=$(curl -s http://169.254.169.254/latest/meta-data/public-ipv4)
export PASH_PUBLIC_IP="$PUBLIC_IP"
```

Set custom STUN server (if supported):

```bash
export PASH_STUN_SERVER=stun.l.google.com:19302
```

## Permanent Fix Plan

### Step 1: Make STUN helper non-fatal

File: `/home/ubuntu/splash-stun-lib/src/stun_helper.rs`

Changes:

- Change `get_addr()` return type to `Result<(String, String), String>`.
- Replace panic-based unwrapping with error-aware handling.
- Replace assertions with explicit validation + error returns.
- Add retry/backoff around STUN request.

Sketch:

```rust
pub async fn get_addr() -> Result<(String, String), String> {
    for attempt in 1..=5 {
        // send STUN request
        let event = handler_rx.recv().await.ok_or("stun channel closed")?;
        match event.event_body {
            Ok(msg) => {
                // validate response + parse xor mapped address
                return Ok((local_addr, external_addr));
            }
            Err(e) => {
                eprintln!("stun attempt {} failed: {:?}", attempt, e);
                tokio::time::sleep(Duration::from_millis(200 * attempt)).await;
            }
        }
    }
    Err("STUN failed after retries".to_string())
}
```

### Step 2: Handle STUN failure at caller

File: `/home/ubuntu/splash-stun-lib/bin/pashlib-oneproc/main.rs`

Changes:

- Replace direct `get_addr().await` usage with `match`.
- On STUN failure, call fallback address resolver instead of panicking.

Sketch:

```rust
let (local_addr_global, external_addr_global) = if no_holepunch_enabled() {
    (String::new(), String::new())
} else {
    match stun_helper::get_addr().await {
        Ok((local_addr, external_addr)) => (local_addr, external_addr),
        Err(e) => {
            eprintln!("STUN discovery failed: {}", e);
            resolve_public_addr_fallback().await
        }
    }
};
```

### Step 3: Add robust fallback address resolution

File: `/home/ubuntu/splash-stun-lib/src/holepunch.rs`

Add:

- `resolve_public_addr_fallback()` that:
  - binds UDP socket to obtain local address,
  - uses existing env-based public-address resolution,
  - optionally queries AWS metadata (`169.254.169.254`) for public IPv4.

Optional helper:

- `query_aws_metadata_public_ip()` with short timeout and defensive parsing.

### Step 4: Optional public API export

File: `/home/ubuntu/splash-stun-lib/src/lib.rs`

If needed by your import style, re-export fallback helpers.

## Testing and Verification

### Development checks

```bash
cd /home/ubuntu/splash-stun-lib
cargo run --bin stun-test
```

```bash
cd /home/ubuntu/pash/evaluation/benchmarks/unix50
./run-leash-benchmark-matrix.sh --small --noopt
```

### Verify fallback behavior

- Confirm logs show STUN failure message but no panic.
- Confirm fallback address selection log appears.
- Confirm benchmark completes and outputs are produced.

### Optional failure-injection test

If you intentionally block STUN traffic, ensure pashlib no longer panics and falls back.

## Build and Deployment

### 1) Build updated binary

```bash
cd /home/ubuntu/splash-stun-lib
cargo clean
cargo build --release --bin pashlib
ls -lh target/release/pashlib
```

### 2) Deploy locally for PaSh runtime

```bash
cp /home/ubuntu/splash-stun-lib/target/release/pashlib /home/ubuntu/pash/runtime/serverless/runtime/pashlib
chmod +x /home/ubuntu/pash/runtime/serverless/runtime/pashlib
/home/ubuntu/pash/runtime/serverless/runtime/pashlib --help
```

### 3) Lambda layer update (if Lambda path is used)

```bash
cd /home/ubuntu/pash/runtime/serverless
cp /home/ubuntu/splash-stun-lib/target/release/pashlib runtime/pashlib
cd runtime
zip -r ../lambda-layer.zip .
```

Then publish/update your Lambda layer and attach it to the function version used for LEASH.

## Why This Fix Is Worth Keeping

- Non-breaking for environments where STUN works.
- Graceful degradation when STUN fails.
- Better diagnostics in logs.
- Cleaner ops posture for mixed EC2/Lambda networking conditions.

## Related Docs

- `LEASH.md` for holepunch architecture context.
- `S3_DIRECT_STREAMING_METHODS.md` for non-holepunch S3 execution paths.
