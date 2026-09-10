#!/bin/bash
set -u
cd "$(dirname "$0")" || exit 1
: "${PASH_TOP:?PASH_TOP not set}"
: "${AWS_BUCKET:?AWS_BUCKET not set}"

NBASE=${NBASE:-50}
NTOTAL=${NTOTAL:-1000}
SIZE_MB=${SIZE_MB:-100}
S3_PREFIX="agent/inputs/logs"
AVG_LINE=64                                   # ~bytes/line, for the line count
N=$(( SIZE_MB * 1024 * 1024 / AVG_LINE ))

STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

# Deterministic generator (self-contained LCG; identical bytes on any awk).
gen_awk='
BEGIN {
    split("ERROR WARNING INFO INFO WARNING INFO ERROR DEBUG INFO WARNING", LV, " "); nlv = 10;
    msg[0]="Database connection established"; msg[1]="Deadlock detected in transaction ID";
    msg[2]="Disk space low: remaining"; msg[3]="Scheduled backup completed successfully";
    msg[4]="API response time exceeded threshold"; msg[5]="User login successful for user";
    msg[6]="Unhandled exception: TimeoutError"; msg[7]="Cache cleared successfully";
    msg[8]="High memory usage detected"; msg[9]="Slow query detected: execution time"; nmsg = 10;
    state = SEED;
    for (i = 0; i < N; i++) {
        state=(state*1103515245+12345)%2147483648; hh=int(state/60)%24;
        state=(state*1103515245+12345)%2147483648; mm=state%60;
        state=(state*1103515245+12345)%2147483648; ss=state%60;
        state=(state*1103515245+12345)%2147483648; lv=LV[(state%nlv)+1];
        state=(state*1103515245+12345)%2147483648; m=msg[state%nmsg];
        state=(state*1103515245+12345)%2147483648; nnum=state%10000;
        printf "%s %02d:%02d:%02d [%s] %s %d\n", DATE, hh, mm, ss, lv, m, nnum;
    }
}'

# Build the NTOTAL target names (valid date prefix so filedate=${f%%_*} works).
srcs=(api app auth db)
names=()
for i in $(seq 1 "$NTOTAL"); do
    s=${srcs[$(( (i-1) % 4 ))]}
    printf -v nm "2025-08-12_%s_%04d.log" "$s" "$i"
    names+=("$nm")
done
printf "%s\n" "${names[@]}" > logs.txt
echo "wrote logs.txt with ${#names[@]} names"

# 1+2. Generate + upload the NBASE distinct base files.
echo "generating + uploading $NBASE base files (~${SIZE_MB}MB each, N=$N lines)..."
for i in $(seq 0 $((NBASE-1))); do
    nm="${names[$i]}"
    awk -v SEED="$((i+1))" -v N="$N" -v DATE="2025-08-12" "$gen_awk" > "$STAGE/$nm"
    aws s3 cp "$STAGE/$nm" "s3://$AWS_BUCKET/$S3_PREFIX/$nm" --no-progress >/dev/null
    rm -f "$STAGE/$nm"
    echo "  [$((i+1))/$NBASE] uploaded $nm"
done

# 3. Fan out the rest with server-side copies (no data leaves S3).
echo "creating $((NTOTAL-NBASE)) server-side copies..."
for i in $(seq "$NBASE" $((NTOTAL-1))); do
    nm="${names[$i]}"
    base="${names[$(( i % NBASE ))]}"
    aws s3 cp "s3://$AWS_BUCKET/$S3_PREFIX/$base" "s3://$AWS_BUCKET/$S3_PREFIX/$nm" --no-progress >/dev/null
    if [ $(( (i+1) % 100 )) -eq 0 ]; then echo "  copied up to $((i+1))/$NTOTAL"; fi
done

echo "done: $NTOTAL objects ($NBASE distinct, ${SIZE_MB}MB each) under s3://$AWS_BUCKET/$S3_PREFIX/"
