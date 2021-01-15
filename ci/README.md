This directory holds a draft implementation of a thin CI controller
used to programmatically orchastrate tests/reporting on EC2 instances.

`main.js` is the entry point. Other modules are supplemental.

The existing implementation uses `scripts/webhooks.js` to react to
requests by running Bash subprocesses. This version uses similar
dispatching rules, but it does not use the host system to run tests.
EC2 instances are responsible for the heavy lifting.

Performance tests are treated differently than correctness tests, in
that the former has higher compute requirements and are therefore more
expensive to run. All the nuance in how to treat tests will
evolve. The short-term focus is to reproduce correctness tests using a
long-running EC2 instance.



## Changelog

### Jan 15 2021

- Created `ci` directory
- Wrote `install-aws-cli-v2.sh`, a helper script to install AWS CLI v2
- Wrote draft HMAC comparison to protect the server.
- Refactored code to control more functionality solely in terms of arguments.
