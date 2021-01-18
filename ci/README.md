This is a CI controller used to run tests using AWS SSM.

To start the server, run `node main.js`.


## Changelog

### Jan 18 2021

- Move over active code for route implementations.
- Call aws ec2 instead of local Bash scripts


### Jan 15 2021

- Created `ci` directory
- Wrote `install-aws-cli-v2.sh`, a helper script to install AWS CLI v2
- Wrote draft HMAC comparison to protect the server.
- Refactored code to control more functionality solely in terms of arguments.
