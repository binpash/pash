cd $(dirname "$0")

output_dir="outputs-old"

SCRIPTS=(
    1_1.sh
    2_1.sh
    2_2.sh
    3_1.sh
    3_2.sh
    3_3.sh
    4_3.sh
    6_1_1.sh
    6_1_2.sh
    6_2.sh
    6_3.sh
    6_4.sh
    6_5.sh
    6_7.sh
    7_1.sh
    7_2.sh
    8.2_1.sh
    6_1.sh
    8_1.sh
    8.3_2.sh
    8.3_3.sh
  )

for script in "${SCRIPTS[@]}"
do
    outfile_analysis="${output_dir}/leash-${script}-pash-w1-log-analysis.log"
    outfile_time="${output_dir}/leash-${script}-pash-w1-time.log"

    # --- parse real time (e.g., 4m46.421s) ---
    real_line=$(grep '^real' "$outfile_time")
    mins=$(echo "$real_line" | awk -F'm' '{print $1}' | awk '{print $2}')
    secs=$(echo "$real_line" | awk -F'm' '{print $2}' | sed 's/s//')
    real_seconds=$(awk -v m="$mins" -v s="$secs" 'BEGIN {print m*60+s}')

    # --- parse billed time ---
    billed_ms=$(grep -oP 'Total billed time:\s*\K[0-9]+' "$outfile_analysis")

    # --- output tuple ---
    echo "[$script $real_seconds, $billed_ms]"
done
