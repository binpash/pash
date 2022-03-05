#! /usr/bin/env bash

# Fire hook script in response to an event. Hooks are format-specific,
# so Debian packages can respond differently than Arch packages, for
# example.

set -eu
event_name="$1"
output_format="$2"

here="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
event="$here/hooks/${event_name}"

if [ -f "$event/hook" ]; then
    # Pass directory name so the script don't need to figure out where
    # it is.  -f implies the script does not have to be Bash, but use
    # of non-executables must raise an error.
    "$event/hook" "$event" "$output_format"
fi
