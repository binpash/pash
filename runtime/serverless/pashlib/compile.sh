#!/bin/bash

cd "$(dirname "$0")"

cargo build --release

cp target/release/pashlib ../runtime/pashlib # Copy the compiled pashlib to the runtime directory for packaging.
cp target/release/pashlib $PASH_TOP/runtime/pashlib # Copy the compiled pashlib to the ec2 runtime directory