#!/usr/bin/env dgsh
#
# SYNOPSIS C/C++ symbols that should be static
# DESCRIPTION
# Given as an argument a directory containing object files, show which
# symbols are declared with global visibility, but should have been
# declared with file-local (static) visibility instead.
# Demonstrates the use of dgsh-capable comm (1) to combine data from
# two sources.
#
#  Copyright 2014 Diomidis Spinellis
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

# Find object files
find "$1" -name \*.o |

# Print defined symbols
xargs nm |

tee |
{{

  # List all defined (exported) symbols
  awk 'NF == 3 && $2 ~ /[A-Z]/ {print $3}' | sort

  # List all undefined (imported) symbols
  awk '$1 == "U" {print $2}' | sort

}} |
# Print exports that are not imported
comm -23
