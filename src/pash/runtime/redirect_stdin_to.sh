#!/usr/bin/env bash

## TODO: Is this a hack?
{ cat > "${1?No file to redirect to}" <&3 3<&- & } 3<&0
