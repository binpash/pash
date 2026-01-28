# STUN Timeout Panic During Holepunch

This doc explains the `ErrTransactionTimeOut` panic seen when running
`run-leash-compare-generic.sh`, and the steps to fix it.

## Symptoms

You see output like:

```text
Get local_addr: 172.31.12.152:46604 external_addr: 3.133.71.249:46604
thread 'main' panicked at /code/src/stun_helper.rs:30:32:
called `Result::unwrap()` on an `Err` value: ErrTransactionTimeOut
```

The panic happens inside the STUN discovery step, before any holepunch TCP
connection attempts are made.

## Root Cause Analysis

The STUN helper calls `stun.l.google.com:19302` over UDP and waits for a
binding response. If no response arrives before the STUN client times out,
the event handler yields `ErrTransactionTimeOut`.

In `/home/ubuntu/splash-stun-lib/src/stun_helper.rs`, we do:

```rust
let event = handler_rx.recv().await.unwrap();
let msg = event.event_body.unwrap();
```

When the STUN client times out, `event.event_body` is `Err(...)`, so the
`unwrap()` panics. The earlier "Get local_addr..." lines indicate some STUN
transactions succeeded; a later one timed out and crashed the process.

Common underlying causes:

- The instance cannot reach the public STUN server (no NAT gateway, blocked
  UDP egress, or restrictive NACL/security group).
- DNS resolution or UDP packet loss causes dropped STUN responses.
- STUN server rate limiting or transient network issues.

## Fix Steps

### 1) Confirm STUN reachability

Run the bundled STUN test:

```console
cd /home/ubuntu/splash-stun-lib
cargo run --bin stun-test
```

If this prints a local and external address, STUN is reachable. If it hangs
or errors, STUN is failing at the network layer.

### 2) Validate network egress

On EC2, confirm:

- The instance has a route to the Internet (IGW for public subnet or NAT
  gateway for private subnet).
- Security group egress allows UDP to port 19302.
- NACLs allow outbound UDP and the return path.

If you are in a private subnet without NAT, STUN will always fail.

### 3) Make STUN failures non-fatal (recommended)

Change the STUN helper to handle timeouts and retry instead of panicking.
Suggested approach:

- Add a retry loop around the STUN request.
- Increase the STUN client RTO (retransmission timeout).
- Convert the function to return `Result<...>` so callers can decide whether
  to fall back or fail gracefully.

Pseudo-code sketch:

```rust
for attempt in 1..=5 {
    let mut client = ClientBuilder::new()
        .with_conn(cc.clone())
        .with_rto(Duration::from_secs(1))
        .build()?;

    client.send(&msg, Some(Arc::new(handler_tx))).await?;
    let event = handler_rx.recv().await.ok_or("stun channel closed")?;

    match event.event_body {
        Ok(msg) => { /* success */ }
        Err(e) => {
            eprintln!("stun attempt {} failed: {}", attempt, e);
            tokio::time::sleep(Duration::from_millis(200 * attempt)).await;
        }
    }
}
```

### 4) Add configurability for the STUN server (optional)

Allow overriding the STUN server with an env var, for example:

```
PASH_STUN_SERVER=stun.l.google.com:19302
```

This helps when the default server is blocked or rate limited.

### 5) Short-term workaround if you can avoid NAT traversal

If holepunching is not required for your setup, disable it:

```console
export PASH_NO_HOLEPUNCH=1
```

This bypasses STUN entirely. Use only if direct connectivity is acceptable.

## Verification

After applying the changes:

1. Run `cargo run --bin stun-test` again to confirm STUN succeeds or retries.
2. Re-run the benchmark script and confirm no STUN-related panics.
3. If STUN still fails, fall back to `PASH_NO_HOLEPUNCH=1` or fix networking.
