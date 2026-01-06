# Ensure directories exist
mkdir -p /tmp/pash_8Z86Ixy/ba360f5a7b9b4067b3943eb10f5ec080
mkdir -p /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0

# Remove any stale files
rm -f \
  /tmp/pash_8Z86Ixy/ba360f5a7b9b4067b3943eb10f5ec080/#fifo2 \
  /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/#fifo10 \
  /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/#fifo14 \
  /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/#fifo29

# Recreate FIFOs
mkfifo /tmp/pash_8Z86Ixy/ba360f5a7b9b4067b3943eb10f5ec080/#fifo2
mkfifo /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/#fifo10
mkfifo /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/#fifo14
mkfifo /tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/#fifo29

SRC_FILE="/tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/fifo29.out"
MID_FILE="/tmp/pash_8Z86Ixy/ba360f5a7b9b4067b3943eb10f5ec080/fifo2.out"
OUT1_FILE="/tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/fifo10.out"
OUT2_FILE="/tmp/pash_8Z86Ixy/6e9e138c3e834eb78417be4b5ac285a0/fifo14.out"

# 1) Download
python3.9 aws/s3-get-object.py \
  "oneliners/inputs/1M.txt" \
  "$SRC_FILE"

# 2) Copy
cat "$SRC_FILE" > "$MID_FILE"

# 3) Split
runtime/r_split -r \
  "$MID_FILE" \
  1000000 \
  "$OUT1_FILE" \
  "$OUT2_FILE"




